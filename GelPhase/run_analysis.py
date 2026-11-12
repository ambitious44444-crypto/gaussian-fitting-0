#!/usr/bin/env python3
"""
=============================================================================
GelPhase: Automated Gel Lane Profile Deconvolution & Threshold Pipeline
Main Batch Execution Script
=============================================================================
"""

import os
import re
import glob
import argparse
import pandas as pd
import numpy as np

from gelphase.deconvolution import fit_gel_lane
from gelphase.threshold import fit_saturation_threshold
from gelphase.visualizer import (
    plot_spectral_decomposition,
    plot_threshold_fit,
    plot_master_comparison,
)

def parse_concentration_from_filename(filename):
    """Extract numeric concentration (e.g., 0, 1, 5, 10) from filename."""
    m = re.search(r'(\d+)%', filename)
    if m:
        return float(m.group(1))
    return -1.0

def run_batch_pipeline(data_root, output_dir):
    """
    Execute full batch analysis on data directory structure.
    Expected structure:
    data_root/
        ├── <Condition1> (e.g. G2 2hr) /
        │     ├── <Group1> (e.g. C2) / *.csv
        │     └── ...
        └── <Condition2> (e.g. G2 4hr) / ...
    """
    os.makedirs(output_dir, exist_ok=True)
    plots_dir = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    print("\n" + "="*70)
    print(" 🧬 GelPhase: Automated Gel Profile & Phase-Transition Analyzer")
    print("="*70)
    print(f"[*] Scanning data directory: {os.path.abspath(data_root)}")
    print(f"[*] Output directory       : {os.path.abspath(output_dir)}\n")
    
    # 1. Discover Conditions and Groups
    condition_folders = [d for d in glob.glob(os.path.join(data_root, "*")) if os.path.isdir(d)]
    if not condition_folders:
        condition_folders = [data_root]
        
    all_spec_records = []
    all_threshold_records = []
    master_comp_data = {}
    
    for cond_path in condition_folders:
        cond_name = os.path.basename(cond_path)
        group_folders = [d for d in glob.glob(os.path.join(cond_path, "*")) if os.path.isdir(d)]
        if not group_folders:
            group_folders = [cond_path]
            
        for grp_path in group_folders:
            grp_name = os.path.basename(grp_path)
            csv_files = sorted(glob.glob(os.path.join(grp_path, "*.csv")))
            if not csv_files:
                continue
                
            print(f">> Processing: [{cond_name}] -> Group [{grp_name}] ({len(csv_files)} files)")
            
            # Sort files by concentration
            file_entries = []
            for f in csv_files:
                conc = parse_concentration_from_filename(os.path.basename(f))
                file_entries.append((conc, f))
            file_entries.sort(key=lambda item: item[0])
            
            grp_concs = []
            grp_sigmas = []
            grp_sig_errs = []
            
            # Subplot directory for this group
            grp_plot_dir = os.path.join(plots_dir, f"{cond_name}_{grp_name}")
            os.makedirs(grp_plot_dir, exist_ok=True)
            
            for conc, fpath in file_entries:
                fname = os.path.basename(fpath)
                try:
                    df = pd.read_csv(fpath, header=None)
                    x = df[0].values
                    y = df[1].values
                    
                    # Deconvolution
                    fit_res = fit_gel_lane(x, y)
                    
                    # Save individual spectral fit plot
                    spec_plot_path = os.path.join(grp_plot_dir, f"fit_{conc}pct.png")
                    plot_spectral_decomposition(
                        fit_res,
                        title=f"{grp_name} ({cond_name}) - {conc}% Conc",
                        output_path=spec_plot_path
                    )
                    
                    grp_concs.append(conc)
                    grp_sigmas.append(fit_res["main_sigma"])
                    grp_sig_errs.append(fit_res["main_sigma_err"])
                    
                    all_spec_records.append({
                        "Condition": cond_name,
                        "Group": grp_name,
                        "Concentration": conc,
                        "Filename": fname,
                        "Best_Model": fit_res["best_model"],
                        "AIC": round(fit_res["best_aic"], 1),
                        "Main_Peak_Mu": round(fit_res["main_mu"], 2),
                        "Main_Peak_Sigma": round(fit_res["main_sigma"], 2),
                        "Main_Peak_Sigma_Err": round(fit_res["main_sigma_err"], 2),
                    })
                except Exception as e:
                    print(f"   [!] Error processing {fname}: {e}")
                    
            # 2. Threshold Broken-Line Regression for this group
            if len(grp_concs) >= 3:
                try:
                    th_res = fit_saturation_threshold(grp_concs, grp_sigmas, grp_sig_errs)
                    
                    th_plot_path = os.path.join(plots_dir, f"threshold_{cond_name}_{grp_name}.png")
                    plot_threshold_fit(
                        th_res,
                        title=f"{grp_name} ({cond_name}) Saturation Threshold",
                        output_path=th_plot_path
                    )
                    
                    all_threshold_records.append({
                        "Condition": cond_name,
                        "Group": grp_name,
                        "Threshold_Breakpoint (%)": round(th_res["threshold"], 2),
                        "Low_Conc_Slope1": round(th_res["slope1"], 3),
                        "High_Conc_Slope2": round(th_res["slope2"], 3),
                        "Baseline_Intercept": round(th_res["y_intercept"], 2),
                    })
                    
                    if grp_name not in master_comp_data:
                        master_comp_data[grp_name] = {}
                    master_comp_data[grp_name][cond_name] = th_res
                    
                    print(f"   --> Found Breakpoint: {th_res['threshold']:.2f}% | Slope1: {th_res['slope1']:.3f} | Slope2: {th_res['slope2']:.3f}")
                except Exception as e:
                    print(f"   [!] Failed threshold fitting for {grp_name}: {e}")

    # 3. Save Summary CSVs
    df_specs = pd.DataFrame(all_spec_records)
    specs_csv = os.path.join(output_dir, "spectral_fit_summary.csv")
    df_specs.to_csv(specs_csv, index=False, encoding="utf-8-sig")
    
    df_thresh = pd.DataFrame(all_threshold_records)
    thresh_csv = os.path.join(output_dir, "threshold_summary.csv")
    df_thresh.to_csv(thresh_csv, index=False, encoding="utf-8-sig")
    
    # 4. Generate Master Comparison Plot
    if master_comp_data:
        master_plot_path = os.path.join(plots_dir, "master_threshold_comparison.png")
        plot_master_comparison(master_comp_data, output_path=master_plot_path)
        print(f"\n[+] Master Comparison Plot saved: {master_plot_path}")

    # 5. Display Final Table
    print("\n" + "="*70)
    print(" 📊 GELPHASE ANALYSIS SUMMARY")
    print("="*70)
    if not df_thresh.empty:
        print(df_thresh.to_string(index=False))
    print("\n[+] All tasks successfully completed!")
    print(f"[+] Output CSVs and Plots generated in: {os.path.abspath(output_dir)}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GelPhase: Automated Gel Lane Profile Deconvolution & Threshold Pipeline")
    parser.add_argument(
        "--data",
        type=str,
        default=os.path.join("..", "非線性擬和"),
        help="Path to root data directory (default: ../非線性擬和)"
    )
    parser.add_argument(
        "--output",
        type=str,
        default="results",
        help="Path to output results directory (default: results)"
    )
    args = parser.parse_args()
    
    run_batch_pipeline(args.data, args.output)
