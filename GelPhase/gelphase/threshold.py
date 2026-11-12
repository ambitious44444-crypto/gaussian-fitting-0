"""
GelPhase Threshold Determination Engine:
Continuous piecewise linear (broken-line) regression for saturation
breakpoint and transition slope extraction.
"""

import numpy as np
from scipy.optimize import curve_fit
from .models import piecewise_linear

def fit_saturation_threshold(concentrations, sigmas, sigma_errors=None, max_conc_bound=15.0):
    """
    Fit continuous piecewise linear broken-line regression to identify
    the saturation threshold (C_th) and transition kinetics.
    
    Parameters:
    -----------
    concentrations : array-like
        List of tested concentrations (e.g. [0, 1, 5, 10]).
    sigmas : array-like
        Extracted band width metric (sigma) for each concentration.
    sigma_errors : array-like, optional
        Standard errors associated with sigma values.
    max_conc_bound : float
        Upper search boundary for concentration threshold.
        
    Returns:
    --------
    threshold_results : dict
        Breakpoint (threshold), baseline intercept, low-concentration slope,
        high-concentration slope, and fitted parameter array.
    """
    c = np.asarray(concentrations, dtype=float)
    s = np.asarray(sigmas, dtype=float)
    
    # Sort by concentration ascending
    sort_idx = np.argsort(c)
    c = c[sort_idx]
    s = s[sort_idx]
    
    if sigma_errors is not None:
        s_err = np.asarray(sigma_errors, dtype=float)[sort_idx]
    else:
        s_err = None

    # Initial parameter guess [x_break, y_intercept, slope1, delta_slope]
    mid_conc = (c.min() + c.max()) / 2.0
    p0 = [mid_conc, s[0], 0.2, 5.0]
    
    bounds = (
        [c.min(), 0.0, -1.0, 0.0],
        [max(max_conc_bound, c.max()), max(50.0, s.max()), 5.0, 30.0]
    )
    
    sigma_weight = s_err if (s_err is not None and np.all(s_err > 0)) else None
    
    popt, pcov = curve_fit(
        piecewise_linear,
        c,
        s,
        p0=p0,
        bounds=bounds,
        sigma=sigma_weight,
        absolute_sigma=True if sigma_weight is not None else False,
        maxfev=10000
    )
    
    x_break = popt[0]
    y_int = popt[1]
    slope1 = popt[2]
    delta_slope = popt[3]
    slope2 = slope1 + delta_slope
    
    # Attempt to extract standard error if degrees of freedom > 0
    try:
        perr = np.sqrt(np.diag(pcov))
        x_break_err = perr[0]
    except Exception:
        x_break_err = np.nan

    return {
        "threshold": x_break,
        "threshold_err": x_break_err,
        "y_intercept": y_int,
        "slope1": slope1,
        "delta_slope": delta_slope,
        "slope2": slope2,
        "popt": popt,
        "pcov": pcov,
        "concentrations": c,
        "sigmas": s,
        "sigma_errors": s_err,
    }
