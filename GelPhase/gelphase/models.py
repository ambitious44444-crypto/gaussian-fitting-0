"""
GelPhase Core Mathematical Models:
- Gaussian peak representation
- Multi-component spectral models with linear baseline
- Akaike Information Criterion (AICc) with small-sample correction
- Continuous Piecewise Linear (Broken-Line) regression model
"""

import numpy as np

def gaussian(x, amp, mu, sigma):
    """
    Standard Gaussian (Normal distribution profile) for spectral peaks.
    
    Parameters:
    -----------
    x : array-like
        Channel, pixel, or distance coordinate along the gel lane.
    amp : float
        Peak amplitude (maximum intensity above local baseline).
    mu : float
        Center position of the peak.
    sigma : float
        Standard deviation (peak width metric).
    """
    return amp * np.exp(-((x - mu) ** 2) / (2.0 * (sigma ** 2)))

def model_2peaks(x, bg_a, bg_b, A1, mu1, sig1, A2, mu2, sig2):
    """
    Two-peak Gaussian model with an underlying linear background.
    Total parameters = 8.
    """
    background = bg_a * x + bg_b
    peak1 = gaussian(x, A1, mu1, sig1)
    peak2 = gaussian(x, A2, mu2, sig2)
    return background + peak1 + peak2

def model_3peaks(x, bg_a, bg_b, A1, mu1, sig1, A2, mu2, sig2, A3, mu3, sig3):
    """
    Three-peak Gaussian model (e.g., Main peak + Shoulder peak + State 2 peak)
    with an underlying linear background. Total parameters = 11.
    """
    background = bg_a * x + bg_b
    peak1 = gaussian(x, A1, mu1, sig1)
    peak2 = gaussian(x, A2, mu2, sig2)
    peak3 = gaussian(x, A3, mu3, sig3)
    return background + peak1 + peak2 + peak3

def calculate_aic(y_true, y_pred, k):
    """
    Calculate Akaike Information Criterion with small-sample bias correction (AICc).
    
    AICc = n * ln(RSS / n) + 2k + [2k(k + 1) / (n - k - 1)]
    
    Parameters:
    -----------
    y_true : array-like
        Observed experimental intensity values.
    y_pred : array-like
        Fitted model intensity values.
    k : int
        Number of estimable parameters in the model.
        
    Returns:
    --------
    aic : float
        Calculated AICc score (lower is better).
    rss : float
        Residual sum of squares.
    """
    n = len(y_true)
    rss = np.sum((y_true - y_pred) ** 2)
    if rss <= 0:
        return -np.inf, 0.0
    
    aic = n * np.log(rss / n) + 2.0 * k
    if n - k - 1 > 0:
        aic += (2.0 * k * (k + 1.0)) / (n - k - 1.0)
    return aic, rss

def piecewise_linear(x, x_break, y_intercept, slope1, delta_slope):
    """
    Continuous Two-Phase Piecewise Linear (Broken-Line) Regression Model:
    
    For x <= x_break:
        y = y_intercept + slope1 * x
    For x > x_break:
        y = y_intercept + slope1 * x + delta_slope * (x - x_break)
        Total slope = slope1 + delta_slope
        
    Parameters:
    -----------
    x : array-like
        Concentration values (e.g., % FBS, % ligand).
    x_break : float
        The breakpoint / saturation threshold concentration (C_th).
    y_intercept : float
        Baseline parameter value at zero concentration.
    slope1 : float
        Slope in the low-concentration pre-transition regime.
    delta_slope : float
        Slope increase beyond the threshold transition.
    """
    return y_intercept + slope1 * x + delta_slope * np.maximum(0.0, x - x_break)
