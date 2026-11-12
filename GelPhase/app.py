import os
import sys
import re
import io

# Ensure the directory containing app.py is in sys.path for Streamlit Cloud
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import streamlit as st

# Import core analytical functions from gelphase package
from gelphase.deconvolution import fit_gel_lane
from gelphase.threshold import fit_saturation_threshold
from gelphase.visualizer import plot_spectral_decomposition, plot_threshold_fit

# =========================================================================
# Page Setup & Configuration
# =========================================================================
st.set_page_config(
    page_title="GelPhase - Automated Gel Profile Analyzer",
    page_icon="🧬",
    layout="wide"
)

st.title("🧬 GelPhase: Automated Gel Profile & Threshold Analyzer")
st.markdown("""
**Welcome to GelPhase!** An open-source computational tool for automated gel electrophoresis lane deconvolution 
and phase-transition saturation threshold identification.
""")

st.info("💡 **Instructions**: Upload your ImageJ line-profile CSV files below. The tool will automatically deconvolve the bands using multi-Gaussian modeling (with AICc arbitration) and fit a continuous piecewise linear saturation threshold.")

# =========================================================================
# File Upload Section (Drag & Drop)
# =========================================================================
uploaded_files = st.file_uploader(
    "📁 Upload ImageJ Lane Profile CSV Files (Select multiple files):",
    type=["csv", "txt"],
    accept_multiple_files=True
)

if uploaded_files:
    st.success(f"Successfully selected {len(uploaded_files)} file(s).")
    
    st.markdown("### 1. Confirm Concentrations for Each File")
    st.caption("GelPhase attempts to auto-detect concentrations from filenames (e.g. '...5%...'). Please verify or adjust below:")
    
    file_concs = []
    cols = st.columns(min(4, len(uploaded_files)))
    
    for idx, f in enumerate(uploaded_files):
        # Auto-detect concentration percentage from filename
        m = re.search(r'(\d+(?:\.\d+)?)%', f.name)
        detected_conc = float(m.group(1)) if m else float(idx)
        
        with cols[idx % 4]:
            val = st.number_input(
                f"Conc for {f.name[:20]}...",
                value=detected_conc,
                step=0.5,
                key=f"file_conc_{idx}"
            )
            file_concs.append((val, f))
            
    st.markdown("---")
    
    # Run Button
    if st.button("🚀 Run GelPhase Deconvolution & Threshold Fitting", type="primary"):
        # Sort files by concentration
        file_concs.sort(key=lambda x: x[0])
        
        progress = st.progress(0)
        status = st.empty()
        
        extracted_data = []
        spectral_figures = []
        
        for i, (conc, file_obj) in enumerate(file_concs):
            status.text(f"Processing ({i+1}/{len(file_concs)}): {file_obj.name} @ {conc}%...")
            file_obj.seek(0)
            
            try:
                # Read 2-column CSV (Column 0: pixel/distance x, Column 1: intensity y)
                df = pd.read_csv(file_obj, header=None)
                x = df[0].values
                y = df[1].values
                
                # Run Multi-Gaussian Deconvolution
                fit_res = fit_gel_lane(x, y)
                fig_spec = plot_spectral_decomposition(
                    fit_res,
                    title=f"{file_obj.name} (Conc: {conc}%)"
                )
                spectral_figures.append((conc, fig_spec))
                
                extracted_data.append({
                    "Filename": file_obj.name,
                    "Concentration (%)": conc,
                    "Best Model (AICc)": fit_res["best_model"],
                    "AICc Score": round(fit_res["best_aic"], 1),
                    "State 1 Center (mu)": round(fit_res["main_mu"], 2),
                    "State 1 Width (sigma)": round(fit_res["main_sigma"], 2),
                    "Width StdErr": round(fit_res["main_sigma_err"], 2),
                })
            except Exception as e:
                st.error(f"Error processing {file_obj.name}: {e}")
                
            progress.progress((i + 1) / len(file_concs))
            
        status.empty()
        progress.empty()
        
        if extracted_data:
            df_results = pd.DataFrame(extracted_data)
            
            st.markdown("### 2. Extracted Band Parameters")
            st.dataframe(df_results, use_container_width=True)
            
            # CSV Download Button
            csv_buffer = io.StringIO()
            df_results.to_csv(csv_buffer, index=False)
            st.download_button(
                label="📥 Download Extracted Parameters (CSV)",
                data=csv_buffer.getvalue(),
                file_name="gelphase_extracted_parameters.csv",
                mime="text/csv"
            )
            
            st.markdown("---")
            st.markdown("### 3. Saturation Threshold Analysis (Piecewise Broken-Line Fit)")
            
            concs = df_results["Concentration (%)"].values
            sigmas = df_results["State 1 Width (sigma)"].values
            sig_errs = df_results["Width StdErr"].values
            
            if len(concs) >= 3:
                try:
                    th_res = fit_saturation_threshold(concs, sigmas, sig_errs)
                    fig_th = plot_threshold_fit(th_res, title="Fitted Saturation Threshold Curve")
                    
                    m1, m2, m3, m4 = st.columns(4)
                    m1.metric("🎯 Saturation Threshold (C_th)", f"{th_res['threshold']:.2f} %")
                    m2.metric("📉 Pre-Transition Slope (s1)", f"{th_res['slope1']:.3f}")
                    m3.metric("📈 Post-Transition Slope (s2)", f"{th_res['slope2']:.3f}")
                    m4.metric("📊 Baseline Intercept (y0)", f"{th_res['y_intercept']:.2f}")
                    
                    st.pyplot(fig_th)
                    
                    # Image Download Button
                    img_buffer = io.BytesIO()
                    fig_th.savefig(img_buffer, format="png", dpi=300)
                    st.download_button(
                        label="📥 Download Threshold Plot (PNG, 300 DPI)",
                        data=img_buffer.getvalue(),
                        file_name="gelphase_saturation_threshold.png",
                        mime="image/png"
                    )
                except Exception as e:
                    st.error(f"Failed to fit saturation threshold: {e}")
            else:
                st.warning("⚠️ At least 3 concentration points are required to compute a saturation threshold.")
                
            st.markdown("---")
            st.markdown("### 4. Deconvoluted Spectral Bands (Individual Fits)")
            spec_cols = st.columns(2)
            for j, (conc, fig) in enumerate(spectral_figures):
                with spec_cols[j % 2]:
                    st.pyplot(fig)
else:
    st.write("👆 Please upload CSV files above to begin analysis.")

st.markdown("---")
st.caption("GelPhase | Developed for BIOMOD Nanotechnology Research | Licensed under MIT")
