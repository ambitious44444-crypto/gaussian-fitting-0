"""
GelPhase Visualizer:
High-resolution, publication-ready plotting modules for spectral
deconvolution, piecewise threshold identification, and multi-condition comparisons.
"""

import os
import numpy as np
import matplotlib
# Headless safety backend
if not os.environ.get("MPLCONFIGDIR"):
    os.environ["MPLCONFIGDIR"] = os.path.abspath(".mplconfig")
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .models import gaussian, model_2peaks, model_3peaks, piecewise_linear

def plot_spectral_decomposition(fit_result, title="Gel Lane Deconvolution", output_path=None, dpi=300):
    """
    Plot single-spectrum decomposed Gaussian components and baseline.
    """
    x = fit_result["x"]
    y = fit_result["y"]
    popt = fit_result["best_popt"]
    model_name = fit_result["best_model"]
    
    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=dpi)
    ax.scatter(x, y, color="#2c3e50", s=8, alpha=0.5, label="Raw Data (ImageJ Profile)")
    
    bg = popt[0] * x + popt[1]
    ax.plot(x, bg, "k--", alpha=0.6, label="Linear Baseline")
    
    if model_name == "2-Peak":
        p1 = gaussian(x, popt[2], popt[3], popt[4])
        p2 = gaussian(x, popt[5], popt[6], popt[7])
        y_fit = bg + p1 + p2
        ax.plot(x, y_fit, color="#e74c3c", lw=2.0, label=f"Total Fit (AIC={fit_result['best_aic']:.0f})")
        ax.plot(x, bg + p1, color="#2980b9", lw=1.5, label=f"State 1 ($\mu$={popt[3]:.1f}, $\sigma$={popt[4]:.2f})")
        ax.plot(x, bg + p2, color="#27ae60", lw=1.5, label=f"State 2 ($\mu$={popt[6]:.1f}, $\sigma$={popt[7]:.2f})")
    else:
        p1 = gaussian(x, popt[2], popt[3], popt[4])
        p2 = gaussian(x, popt[5], popt[6], popt[7])
        p3 = gaussian(x, popt[8], popt[9], popt[10])
        y_fit = bg + p1 + p2 + p3
        ax.plot(x, y_fit, color="#e74c3c", lw=2.0, label=f"Total Fit (3-Peak, AIC={fit_result['best_aic']:.0f})")
        ax.plot(x, bg + p1, color="#2980b9", lw=1.5, label=f"State 1 ($\mu$={popt[3]:.1f}, $\sigma$={popt[4]:.2f})")
        ax.plot(x, bg + p2, color="#27ae60", lw=1.5, label=f"State 2 ($\mu$={popt[6]:.1f}, $\sigma$={popt[7]:.2f})")
        ax.plot(x, bg + p3, color="#8e44ad", lw=1.2, linestyle=":", label=f"Shoulder ($\mu$={popt[9]:.1f}, $\sigma$={popt[10]:.2f})")
        
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("Distance / Channel (pixels)", fontsize=11)
    ax.set_ylabel("Intensity (a.u.)", fontsize=11)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, linestyle=":", alpha=0.5)
    fig.tight_layout()
    
    if output_path:
        fig.savefig(output_path, dpi=dpi)
        plt.close(fig)
    return fig

def plot_threshold_fit(th_result, title="Saturation Threshold Analysis", output_path=None, dpi=300):
    """
    Plot concentration vs peak width with continuous piecewise linear fit.
    """
    c = th_result["concentrations"]
    s = th_result["sigmas"]
    s_err = th_result["sigma_errors"]
    th = th_result["threshold"]
    popt = th_result["popt"]
    s1 = th_result["slope1"]
    s2 = th_result["slope2"]
    
    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=dpi)
    if s_err is not None:
        ax.errorbar(c, s, yerr=s_err, fmt="o", color="#1f77b4", markersize=6, capsize=4, label=r"Extracted $\sigma_{main}$")
    else:
        ax.scatter(c, s, color="#1f77b4", s=50, label=r"Extracted $\sigma_{main}$")
        
    x_dense = np.linspace(c.min(), c.max(), 300)
    y_dense = piecewise_linear(x_dense, *popt)
    ax.plot(x_dense, y_dense, color="#d62728", lw=2.2,
            label=f"Piecewise Fit\nThreshold: {th:.2f}%\nSlope 1: {s1:.2f}, Slope 2: {s2:.2f}")
    ax.axvline(th, color="gray", linestyle="--", alpha=0.8, label=f"Breakpoint = {th:.2f}%")
    
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xlabel("Concentration (%)", fontsize=11)
    ax.set_ylabel(r"Band Width ($\sigma_{main}$)", fontsize=11)
    ax.legend(fontsize=9, loc="upper left")
    ax.grid(True, linestyle=":", alpha=0.6)
    fig.tight_layout()
    
    if output_path:
        fig.savefig(output_path, dpi=dpi)
        plt.close(fig)
    return fig

def plot_master_comparison(master_data, output_path=None, dpi=300):
    """
    Plot multi-condition side-by-side threshold comparisons.
    master_data: dict of {group_name: {'cond1': th_result1, 'cond2': th_result2}}
    """
    n_groups = len(master_data)
    fig, axes = plt.subplots(1, n_groups, figsize=(6 * n_groups, 5), dpi=dpi)
    if n_groups == 1:
        axes = [axes]
        
    colors = ["#27ae60", "#2980b9", "#8e44ad", "#e67e22"]
    
    for idx, (grp, cond_dict) in enumerate(master_data.items()):
        ax = axes[idx]
        for c_idx, (cond_name, res) in enumerate(cond_dict.items()):
            c = res["concentrations"]
            s = res["sigmas"]
            th = res["threshold"]
            popt = res["popt"]
            col = colors[c_idx % len(colors)]
            
            x_dense = np.linspace(c.min(), c.max(), 250)
            ax.scatter(c, s, color=col, s=40, label=f"{cond_name} ($C_{{th}}={th:.2f}\%$)")
            ax.plot(x_dense, piecewise_linear(x_dense, *popt), color=col, lw=1.8, linestyle="--")
            
        ax.set_title(f"{grp} Threshold Comparison", fontsize=12, fontweight="bold")
        ax.set_xlabel("Concentration (%)", fontsize=11)
        ax.set_ylabel(r"Band Width ($\sigma_{main}$)", fontsize=11)
        ax.legend(fontsize=9, loc="upper left")
        ax.grid(True, linestyle=":", alpha=0.6)
        
    fig.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=dpi)
        plt.close(fig)
    return fig
