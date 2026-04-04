"""Plotting style configuration for CFRAM diagnostics.

Matches the visual style of Wu et al. (2025) J. Climate figures.
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


# Nonlinear contour levels matching Fig.3 (partial temperature changes)
LEVELS_DT = [-15, -7, -2, -0.5, 0, 0.5, 2, 7, 15]

# Contour levels for aerosol subtypes (Fig.5)
LEVELS_AER = [-3, -1.5, -0.5, 0, 0.5, 1.5, 3]

# Default colormap: diverging red-blue
CMAP_DIVERGE = 'RdBu_r'

# Font sizes
FONTSIZE_TITLE = 12
FONTSIZE_LABEL = 10
FONTSIZE_TICK = 8
FONTSIZE_PANEL = 11


def setup_style():
    """Set up matplotlib defaults for publication-quality figures."""
    plt.rcParams.update({
        'font.size': FONTSIZE_LABEL,
        'axes.titlesize': FONTSIZE_TITLE,
        'axes.labelsize': FONTSIZE_LABEL,
        'xtick.labelsize': FONTSIZE_TICK,
        'ytick.labelsize': FONTSIZE_TICK,
        'legend.fontsize': FONTSIZE_TICK,
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight',
        'font.family': 'sans-serif',
    })


def make_norm(levels):
    """Create a BoundaryNorm for nonlinear contour levels."""
    return mpl.colors.BoundaryNorm(levels, ncolors=256, clip=True)
