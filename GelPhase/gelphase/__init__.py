"""
GelPhase: Open-Source Biomolecular Gel Lane Profile Deconvolution
and Phase-Transition Saturation Threshold Analyzer.
"""

from .models import (
    gaussian,
    model_2peaks,
    model_3peaks,
    calculate_aic,
    piecewise_linear,
)
from .deconvolution import fit_gel_lane
from .threshold import fit_saturation_threshold
from .visualizer import (
    plot_spectral_decomposition,
    plot_threshold_fit,
    plot_master_comparison,
)

__version__ = "1.0.0"
__author__ = "BIOMOD Research Team & GelPhase Contributors"
__all__ = [
    "gaussian",
    "model_2peaks",
    "model_3peaks",
    "calculate_aic",
    "piecewise_linear",
    "fit_gel_lane",
    "fit_saturation_threshold",
    "plot_spectral_decomposition",
    "plot_threshold_fit",
    "plot_master_comparison",
]
