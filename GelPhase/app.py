import os
import io
import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import streamlit as st

from gelphase.deconvolution import fit_gel_lane
from gelphase.threshold import fit_saturation_threshold
from gelphase.visualizer import plot_spectral_decomposition, plot_threshold_fit

# Page configuration
st.set_page_config(
    page_title="GelPhase: Biomolecular Gel Deconvolution & Threshold Analyzer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom Styling
st.markdown("""
<style>
    .main-header {
        font-size: 2.2rem;
        font-weight: 700;
        color: #1e3d59;
        margin-bottom: 0.2rem;
    }
    .sub-header {
        font-size: 1.1rem;
        color: #17b978;
        font-weight: 500;
        margin-bottom: 1.5rem;
    }
    .metric-card {
        background-color: #f8f9fa;
        border-radius: 8px;
        padding: 15px;
        border-left: 5px solid #17b978;
    }
</style>
""", unsafe_allow_html=True)

st.markdown('<div class="main-header">🧬 GelPhase: Biomolecular Gel Deconvolution & Threshold Analyzer</div>', unsafe_allow_html=True)
st.markdown('<div class="sub-header">Automated Gel Electrophoresis Lane Deconvolution & Saturation Phase-Transition Identification</div>', unsafe_allow_html=True)

# Sidebar
st.sidebar.image("https://img.shields.io/badge/BIOMOD-Open--Source-orange?style=for-the-badge", use_container_width=True)
st.sidebar.title("Navigation & Settings")
app_mode = st.sidebar.radio(
    "Choose Mode:",
    ["📊 Benchmark Dashboard (Pre-loaded Data)", "🧪 Custom Gel Analysis (Upload CSV)", "📐 Mathematical Theory & Methods", "📝 Citation & Wiki Attribution"]
)

st.sidebar.markdown("---")
st.sidebar.info("""
**GelPhase Pipeline**:
1. Multi-peak Gaussian deconvolution
2. AICc-based model arbitration
3. Continuous broken-line threshold regression
""")

# =========================================================================
# Mode 1: Benchmark Dashboard
# =========================================================================
if app_mode == "📊 Benchmark Dashboard (Pre-loaded Data)":
    st.subheader("📊 BIOMOD Serum Stability Benchmark Results (FBS Assay)")
    st.markdown("""
    This dashboard presents the experimental validation across **24 gel lane profiles** evaluated under two reaction times (2hr vs 4hr)
    and biological replicates (C2, E2-1, E2-2).
    """)
    
    results_dir = os.path.join(os.path.dirname(__file__), "results")
    thresh_csv = os.path.join(results_dir, "threshold_summary.csv")
    master_plot = os.path.join(results_dir, "plots", "master_threshold_comparison.png")
    
    col1, col2 = st.columns([1, 1])
    
    if os.path.exists(thresh_csv):
        df_th = pd.read_csv(thresh_csv)
        with col1:
            st.markdown("### 🎯 Critical Saturation Threshold Summary")
            st.dataframe(df_th, use_container_width=True)
            
            st.markdown("""
            > **Key Biological Finding**:
            > * Four independent datasets (**2hr C2, 2hr E2-1, 4hr C2, 4hr E2-1**) tightly converge on a **~4.7% - 5.3% FBS breakdown threshold**!
            > * Post-transition kinetic slopes ($s_2$) surge 20–80x higher than pre-transition slopes ($s_1$).
            """)
    else:
        st.warning("threshold_summary.csv not found in results directory.")

    with col2:
        st.markdown("### 📈 Master Threshold Comparison (2hr vs 4hr)")
        if os.path.exists(master_plot):
            st.image(master_plot, use_container_width=True, caption="Multi-Condition Threshold Transition Kinetics")
        else:
            st.info("Master plot image not found.")

    st.markdown("---")
    st.markdown("### 🔍 Interactive Condition & Group Inspector")
    
    cond_choice = st.selectbox("Select Condition:", ["G2 2hr", "G2 4hr"])
    grp_choice = st.selectbox("Select Sample Group:", ["C2", "E2-1", "E2-2"])
    
    plot_folder = os.path.join(results_dir, "plots", f"{cond_choice}_{grp_choice}")
    thresh_fig = os.path.join(results_dir, "plots", f"threshold_{cond_choice}_{grp_choice}.png")
    
    tab_th, tab_spec = st.tabs(["📉 Saturation Threshold Fit", "🔬 Deconvoluted Spectral Bands (0%, 1%, 5%, 10%)"])
    
    with tab_th:
        if os.path.exists(thresh_fig):
            st.image(thresh_fig, use_container_width=True, caption=f"{cond_choice} {grp_choice} Threshold Fit")
        else:
            st.info(f"Threshold figure for {cond_choice} {grp_choice} not found.")
            
    with tab_spec:
        if os.path.exists(plot_folder):
            c_cols = st.columns(4)
            concentrations = ["0.0", "1.0", "5.0", "10.0"]
            for i, conc in enumerate(concentrations):
                p_file = os.path.join(plot_folder, f"fit_{conc}pct.png")
                with c_cols[i]:
                    if os.path.exists(p_file):
                        st.image(p_file, use_container_width=True, caption=f"Conc: {conc}%")
                    else:
                        st.caption(f"No plot for {conc}%")
        else:
            st.info(f"Deconvolution plots folder {plot_folder} not found.")

# =========================================================================
# Mode 2: Custom Gel Analysis (Upload CSV)
# =========================================================================
elif app_mode == "🧪 Custom Gel Analysis (Upload CSV)":
    st.subheader("🧪 Upload & Analyze Your Own Gel Profiles")
    st.markdown("""
    Upload your ImageJ line-scan CSV files below. 
    Each CSV file must contain two columns: **Column 1 = pixel distance ($x$)**, **Column 2 = intensity ($y$)**.
    """)
    
    uploaded_files = st.file_uploader(
        "Upload ImageJ CSV files (Multiple files allowed)",
        type=["csv", "txt"],
        accept_multiple_files=True
    )
    
    if uploaded_files:
        st.success(f"Uploaded {len(uploaded_files)} file(s).")
        
        # Parse or ask for concentrations
        st.markdown("#### Specify Concentrations for Uploaded Files")
        file_concs = []
        c_cols = st.columns(min(4, len(uploaded_files)))
        
        for idx, file in enumerate(uploaded_files):
            # Try to auto-detect concentration
            m = re.search(r'(\d+(?:\.\d+)?)%', file.name)
            default_val = float(m.group(1)) if m else float(idx)
            
            with c_cols[idx % 4]:
                conc_val = st.number_input(
                    f"{file.name[:20]}...",
                    value=default_val,
                    step=0.5,
                    key=f"conc_{idx}"
                )
                file_concs.append((conc_val, file))
                
        if st.button("🚀 Run GelPhase Deconvolution & Threshold Fitting", type="primary"):
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            extracted_results = []
            spectral_figs = []
            
            # Sort by concentration
            file_concs.sort(key=lambda x: x[0])
            
            for i, (conc, file_obj) in enumerate(file_concs):
                status_text.text(f"Deconvolving: {file_obj.name} (Conc: {conc}%)...")
                file_obj.seek(0)
                df = pd.read_csv(file_obj, header=None)
                x = df[0].values
                y = df[1].values
                
                # Run deconvolution
                fit_res = fit_gel_lane(x, y)
                fig_spec = plot_spectral_decomposition(fit_res, title=f"{file_obj.name} ({conc}%)")
                spectral_figs.append((conc, fig_spec))
                
                extracted_results.append({
                    "Filename": file_obj.name,
                    "Concentration": conc,
                    "Model": fit_res["best_model"],
                    "AICc": round(fit_res["best_aic"], 1),
                    "State1_Center_mu": round(fit_res["main_mu"], 2),
                    "State1_Width_sigma": round(fit_res["main_sigma"], 2),
                    "Width_StdErr": round(fit_res["main_sigma_err"], 2),
                })
                progress_bar.progress((i + 1) / len(file_concs))
                
            status_text.text("Fitting continuous saturation threshold...")
            df_out = pd.DataFrame(extracted_results)
            
            # Fit threshold
            concs = df_out["Concentration"].values
            sigmas = df_out["State1_Width_sigma"].values
            sig_errs = df_out["Width_StdErr"].values
            
            st.markdown("### 🎯 Analysis Results")
            st.dataframe(df_out, use_container_width=True)
            
            if len(concs) >= 3:
                th_res = fit_saturation_threshold(concs, sigmas, sig_errs)
                fig_th = plot_threshold_fit(th_res, title="Fitted Saturation Threshold")
                
                m1, m2, m3 = st.columns(3)
                m1.metric("Critical Breakpoint (C_th)", f"{th_res['threshold']:.2f} %")
                m2.metric("Low-Conc Slope (s1)", f"{th_res['slope1']:.3f}")
                m3.metric("Post-Break Slope (s2)", f"{th_res['slope2']:.3f}")
                
                st.pyplot(fig_th)
            else:
                st.warning("At least 3 concentration points are required to compute a saturation threshold.")
                
            st.markdown("### 🔬 Deconvoluted Profiles")
            spec_cols = st.columns(2)
            for j, (conc, fig) in enumerate(spectral_figs):
                with spec_cols[j % 2]:
                    st.pyplot(fig)

# =========================================================================
# Mode 3: Mathematical Theory & Methods
# =========================================================================
elif app_mode == "📐 Mathematical Theory & Methods":
    st.subheader("📐 Mathematical Modeling & Algorithm Principles")
    st.markdown("""
    GelPhase integrates multi-scale biophysical modeling to ensure objective, reproducible analysis:
    """)
    
    st.markdown("#### 1. Multi-Component Gaussian Deconvolution")
    st.latex(r"y(x) = (\beta_0 + \beta_1 x) + \sum_{i=1}^{M} A_i \exp\left(-\frac{(x - \mu_i)^2}{2\sigma_i^2}\right)")
    st.markdown("""
    * $\mu_i$: Migration distance along lane.
    * $\sigma_i$: Band standard deviation (width), reflecting structural dispersion / degradation.
    * Linear background $\beta_0 + \beta_1 x$ eliminates gel staining gradients.
    """)
    
    st.markdown("#### 2. Model Selection via Corrected AIC (AICc)")
    st.latex(r"\text{AICc} = n \ln\left(\frac{\text{RSS}}{n}\right) + 2k + \frac{2k(k + 1)}{n - k - 1}")
    st.markdown("""
    Implements Occam's Razor: a third (shoulder) peak is adopted if and only if $\Delta\text{AICc} > 2.0$.
    """)
    
    st.markdown("#### 3. Continuous Piecewise Linear (Broken-Line) NLLS Regression")
    st.latex(r"\sigma(C) = y_0 + s_1 C + \Delta s \max(0, C - C_{\text{th}})")
    st.markdown("""
    * **$C_{\text{th}}$**: Phase-transition saturation breakpoint threshold.
    * Solved via **Trust Region Reflective (TRF)** non-linear least squares.
    """)

# =========================================================================
# Mode 4: Citation & Wiki Attribution
# =========================================================================
elif app_mode == "📝 Citation & Wiki Attribution":
    st.subheader("📝 Citation & Attribution Guidelines for BIOMOD / iGEM")
    st.markdown("""
    If you use GelPhase in your competition Wiki, presentation, or scientific publication, please cite as follows:
    """)
    
    st.code("""
@software{gelphase2026,
  author = {BIOMOD Research Team and GelPhase Contributors},
  title = {GelPhase: Automated Biomolecular Gel Lane Profile Deconvolution and Threshold Analyzer},
  year = {2026},
  url = {https://github.com/your-username/GelPhase},
  version = {1.0.0}
}
    """, language="latex")
    
    st.markdown("#### Suggested Wiki Statement (Attributions Page):")
    st.info("""
    "Gel lane profile deconvolution, AICc model selection, and saturation threshold broken-line regressions 
    were performed using the open-source tool **GelPhase** (developed with assistance from Google DeepMind Antigravity). 
    The software is publicly archived on GitHub under the MIT License."
    """)
