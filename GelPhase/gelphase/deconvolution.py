"""
GelPhase Deconvolution Engine:
Automated peak identification, multi-component Gaussian fitting,
and AIC-based model arbitration.
"""

import numpy as np
from scipy.optimize import curve_fit
from .models import model_2peaks, model_3peaks, calculate_aic

def fit_gel_lane(x_data, y_data, band1_range=(50, 160), band2_range=(180, 340), enable_3peaks=True):
    """
    Perform automated deconvolution of a single gel lane intensity profile.
    
    Parameters:
    -----------
    x_data : np.ndarray
        Distance / pixel coordinates along the lane.
    y_data : np.ndarray
        Intensity profile.
    band1_range : tuple
        Approximate search range for State 1 / Main Band (min_x, max_x).
    band2_range : tuple
        Approximate search range for State 2 Band (min_x, max_x).
    enable_3peaks : bool
        Whether to attempt and evaluate a 3-peak model with a shoulder peak.
        
    Returns:
    --------
    result : dict
        Fitting parameters, AIC scores, best model name, extracted sigmas, and errors.
    """
    x = np.asarray(x_data, dtype=float)
    y = np.asarray(y_data, dtype=float)
    
    # 1. Dynamic peak localization
    m1 = (x >= band1_range[0]) & (x <= band1_range[1])
    if np.any(m1):
        p1_x = x[m1][np.argmax(y[m1])]
        p1_y = np.max(y[m1])
    else:
        p1_x = (band1_range[0] + band1_range[1]) / 2.0
        p1_y = np.max(y)
        
    m2 = (x >= band2_range[0]) & (x <= band2_range[1])
    if np.any(m2):
        p2_x = x[m2][np.argmax(y[m2])]
        p2_y = np.max(y[m2])
    else:
        p2_x = (band2_range[0] + band2_range[1]) / 2.0
        p2_y = np.mean(y)
        
    bg_init = np.mean(y[:20]) if len(y) >= 20 else np.min(y)
    
    # 2. Fit Two-Peak Model
    p0_2 = [-0.01, bg_init, p1_y, p1_x, 12.0, p2_y, p2_x, 18.0]
    bounds_2 = (
        [-1.0, -30.0, 0.0, band1_range[0] * 0.8, 1.0, 0.0, band2_range[0] * 0.8, 1.0],
        [1.0, 100.0, max(500.0, p1_y * 3.0), band1_range[1] * 1.2, 80.0, max(500.0, p2_y * 3.0), band2_range[1] * 1.2, 80.0]
    )
    
    popt2, pcov2 = None, None
    aic2, rss2 = np.inf, np.inf
    try:
        popt2, pcov2 = curve_fit(model_2peaks, x, y, p0=p0_2, bounds=bounds_2, maxfev=15000)
        aic2, rss2 = calculate_aic(y, model_2peaks(x, *popt2), len(popt2))
    except Exception as e:
        pass
        
    # 3. Fit Three-Peak Model (Main + Shoulder + Band 2)
    popt3, pcov3 = None, None
    aic3, rss3 = np.inf, np.inf
    if enable_3peaks and popt2 is not None:
        p0_3 = [-0.01, bg_init, p1_y, p1_x, 10.0, p2_y, p2_x, 18.0, p1_y * 0.3, p1_x - 15.0, 8.0]
        bounds_3 = (
            [-1.0, -30.0, 0.0, band1_range[0] * 0.8, 1.0, 0.0, band2_range[0] * 0.8, 1.0, 0.0, 40.0, 1.0],
            [1.0, 100.0, max(500.0, p1_y * 3.0), band1_range[1] * 1.2, 80.0, max(500.0, p2_y * 3.0), band2_range[1] * 1.2, 80.0, max(200.0, p1_y), band1_range[1], 40.0]
        )
        try:
            popt3, pcov3 = curve_fit(model_3peaks, x, y, p0=p0_3, bounds=bounds_3, maxfev=15000)
            aic3, rss3 = calculate_aic(y, model_3peaks(x, *popt3), len(popt3))
        except Exception:
            pass

    # 4. Objective AIC Arbitration (Delta AIC > 2 threshold)
    if (aic2 - aic3) > 2.0 and popt3 is not None:
        best_model = "3-Peak"
        best_popt = popt3
        best_pcov = pcov3
        best_aic = aic3
        best_rss = rss3
        main_mu = popt3[3]
        main_sigma = abs(popt3[4])
        main_sigma_err = np.sqrt(pcov3[4, 4]) if pcov3 is not None else 0.0
    elif popt2 is not None:
        best_model = "2-Peak"
        best_popt = popt2
        best_pcov = pcov2
        best_aic = aic2
        best_rss = rss2
        main_mu = popt2[3]
        main_sigma = abs(popt2[4])
        main_sigma_err = np.sqrt(pcov2[4, 4]) if pcov2 is not None else 0.0
    else:
        raise RuntimeError("Curve fitting failed for both 2-peak and 3-peak models.")

    return {
        "x": x,
        "y": y,
        "best_model": best_model,
        "best_popt": best_popt,
        "best_pcov": best_pcov,
        "best_aic": best_aic,
        "best_rss": best_rss,
        "aic_2p": aic2,
        "aic_3p": aic3,
        "delta_aic": aic2 - aic3,
        "main_mu": main_mu,
        "main_sigma": main_sigma,
        "main_sigma_err": main_sigma_err,
        "popt2": popt2,
        "pcov2": pcov2,
        "popt3": popt3,
        "pcov3": pcov3,
    }
