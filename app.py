import streamlit as st
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

st.set_page_config(page_title="自動多重高斯擬合 (AIC 最佳化)", layout="centered")
st.title("CSV 波形自動分析與 AIC 最佳化擬合")

# 定義多重高斯函數
def multi_gaussian(x, *params):
    y = np.zeros_like(x, dtype=float)
    for i in range(0, len(params), 3):
        A = params[i]
        mu = params[i+1]
        sigma = params[i+2]
        y += A * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))
    return y

# 計算 AIC 與 RSS
def calculate_aic(y_true, y_pred, num_params):
    n = len(y_true)
    rss = np.sum((y_true - y_pred) ** 2)
    # 避免極端情況下 rss 為 0 導致 log(0) 錯誤
    rss = max(rss, 1e-12)
    aic = n * np.log(rss / n) + 2 * num_params
    return aic, rss

uploaded_file = st.file_uploader("請上傳包含 X, Y 座標的 CSV 檔案", type=["csv"])

if uploaded_file is not None:
    df = pd.read_csv(uploaded_file)
    cols = df.columns.tolist()
    
    col1, col2 = st.columns(2)
    x_col = col1.selectbox("選擇 X 軸欄位", cols, index=0)
    y_col = col2.selectbox("選擇 Y 軸欄位", cols, index=1 if len(cols) > 1 else 0)
    
    # 確保依 X 軸排序並去除空值
    clean_df = df[[x_col, y_col]].dropna().sort_values(by=x_col)
    x_data = clean_df[x_col].values
    y_data = clean_df[y_col].values
    
    max_peaks = st.slider("設定最大搜尋波峰上限 (N_max)", min_value=1, max_value=12, value=6)
    
    if st.button("開始自動偵測與 AIC 擬合"):
        with st.spinner("正在評估不同波峰數量的模型並計算 AIC..."):
            # 1. 利用 prominence (波峰顯著性) 預先找出候選波峰並依高度排序
            y_range = np.max(y_data) - np.min(y_data)
            peaks, properties = find_peaks(y_data, prominence=y_range * 0.03)
            
            # 依據波峰顯著性由大到小排序，確保優先擬合主要波峰
            if len(peaks) > 0:
                sorted_idx = np.argsort(properties["prominences"])[::-1]
                peaks = peaks[sorted_idx]
            
            models_record = []
            best_aic = np.inf
            best_popt = None
            best_n = 1
            
            # 2. 遍歷 1 到 max_peaks 進行擬合與 AIC 評估
            for n in range(1, max_peaks + 1):
                p0 = []
                sigma_init = (np.max(x_data) - np.min(x_data)) / (n * 5)
                
                for i in range(n):
                    if i < len(peaks):
                        p0.extend([y_data[peaks[i]], x_data[peaks[i]], sigma_init])
                    else:
                        # 若候選峰不足，從目前殘差最大的位置補上初始猜測
                        if len(p0) > 0:
                            current_y = multi_gaussian(x_data, *p0)
                            residual = y_data - current_y
                            max_res_idx = np.argmax(residual)
                            p0.extend([max(residual[max_res_idx], 1e-3), x_data[max_res_idx], sigma_init])
                        else:
                            p0.extend([np.max(y_data), np.median(x_data), sigma_init])
                
                lower_bounds = [0, np.min(x_data), 1e-6] * n
                upper_bounds = [np.inf, np.max(x_data), np.inf] * n
                
                try:
                    popt, _ = curve_fit(
                        multi_gaussian, x_data, y_data, 
                        p0=p0, bounds=(lower_bounds, upper_bounds), maxfev=15000
                    )
                    y_pred = multi_gaussian(x_data, *popt)
                    # 參數數量 k = 3*n + 1 (包含誤差變異數)
                    aic, rss = calculate_aic(y_data, y_pred, num_params=3*n + 1)
                    
                    models_record.append({"波峰數量 (N)": n, "AIC 分數": round(aic, 2), "殘差平方和 (RSS)": round(rss, 4)})
                    
                    if aic < best_aic:
                        best_aic = aic
                        best_popt = popt
                        best_n = n
                except Exception:
                    # 若該數量無法收斂則跳過
                    continue
            
            if best_popt is not None:
                st.success(f"分析完成！AIC 評估最佳波峰數量為：**{best_n} 個**（最低 AIC = {best_aic:.2f}）")
                
                # 整理最佳模型的參數（依 mu 中心位置由小到大排序，方便閱讀）
                peaks_params = []
                for i in range(best_n):
                    peaks_params.append({
                        "A": best_popt[i*3],
                        "mu": best_popt[i*3+1],
                        "sigma": best_popt[i*3+2]
                    })
                peaks_params.sort(key=lambda item: item["mu"])
                
                res_df = pd.DataFrame([
                    {
                        "高斯分量": f"Peak {i+1}",
                        "振幅 (A)": round(p["A"], 4),
                        "中心位置 (mu)": round(p["mu"], 4),
                        "標準差 (sigma)": round(p["sigma"], 4)
                    }
                    for i, p in enumerate(peaks_params)
                ])
                
                # 繪製雙子圖：上圖為波形擬合結果，下圖為 AIC 模型選擇曲線
                fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7.5), gridspec_kw={'height_ratios': [2, 1]})
                
                # 子圖 1：波形與各高斯分量
                ax1.scatter(x_data, y_data, label="原始數據", color="gray", alpha=0.5, s=15)
                x_smooth = np.linspace(np.min(x_data), np.max(x_data), 600)
                y_fit = multi_gaussian(x_smooth, *best_popt)
                ax1.plot(x_smooth, y_fit, label=f"最佳總擬合曲線 (N={best_n})", color="red", linewidth=2)
                
                for i, p in enumerate(peaks_params):
                    y_sub = multi_gaussian(x_smooth, p["A"], p["mu"], p["sigma"])
                    ax1.plot(x_smooth, y_sub, linestyle="--", label=f"Peak {i+1} (μ={p['mu']:.2f})")
                
                ax1.set_title("多重高斯最佳擬合波形圖")
                ax1.set_xlabel(x_col)
                ax1.set_ylabel(y_col)
                ax1.legend()
                ax1.grid(True, linestyle=":", alpha=0.6)
                
                # 子圖 2：AIC 評估曲線
                aic_df = pd.DataFrame(models_record)
                ax2.plot(aic_df["波峰數量 (N)"], aic_df["AIC 分數"], marker="o", color="steelblue")
                ax2.scatter([best_n], [best_aic], color="red", s=80, zorder=5, label=f"最佳模型 (N={best_n})")
                ax2.set_title("赤池資訊量準則 (AIC) 評估曲線 — 越低越佳")
                ax2.set_xlabel("高斯分量數量 (N)")
                ax2.set_ylabel("AIC 分數")
                ax2.set_xticks(aic_df["波峰數量 (N)"])
                ax2.legend()
                ax2.grid(True, linestyle=":", alpha=0.6)
                
                plt.tight_layout()
                st.pyplot(fig)
                
                # 顯示參數表與 AIC 比較表
                st.subheader(f"最佳模型參數輸出（共 {best_n} 個高斯分佈）")
                st.dataframe(res_df, use_container_width=True)
                
                with st.expander("查看所有波峰數量的 AIC 比較數據"):
                    st.dataframe(aic_df, use_container_width=True)
            else:
                st.error("所有模型皆無法收斂，請確認數據是否包含過多雜訊或負值。")