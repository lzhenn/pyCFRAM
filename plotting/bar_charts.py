"""PAP coefficient bar charts for CFRAM decomposition.

Reproduces Fig.4 of Wu et al. (2025): bar chart showing the relative
contribution (PAP coefficient) of each feedback process to the mean
amplitude of surface temperature anomalies.
"""

import numpy as np
import matplotlib.pyplot as plt
from .style import FONTSIZE_PANEL, FONTSIZE_LABEL, setup_style


def compute_pap(dT_total, dT_components, area_weights=None):
    """Compute Pattern-Amplitude Projection (PAP) coefficients.

    PAP_i = <ΔT> × <ΔT × ΔT_i> / <ΔT²>

    where <·> denotes area-weighted spatial average over the key region.

    Args:
        dT_total: (nlat, nlon) or scalar, total surface temperature change
        dT_components: dict of {name: (nlat, nlon) or scalar}
        area_weights: (nlat, nlon) weights (default: uniform)

    Returns:
        pap: dict of {name: PAP coefficient value}
    """
    if np.isscalar(dT_total) or dT_total.ndim == 0:
        # Single column: PAP simplifies to dT_i
        return {name: float(val) for name, val in dT_components.items()}

    if area_weights is None:
        area_weights = np.ones_like(dT_total)

    w = area_weights / area_weights.sum()

    dt_mean = np.sum(w * dT_total)
    dt2_mean = np.sum(w * dT_total ** 2)

    if dt2_mean == 0:
        return {name: 0.0 for name in dT_components}

    pap = {}
    for name, dT_i in dT_components.items():
        dt_dti_mean = np.sum(w * dT_total * dT_i)
        pap[name] = float(dt_mean * dt_dti_mean / dt2_mean)

    return pap


def compute_cos_weights(lats):
    """Compute cosine latitude area weights."""
    cos_lat = np.cos(np.deg2rad(lats))
    return cos_lat[:, np.newaxis] if cos_lat.ndim == 1 else cos_lat


def plot_pap_bars(pap_values, title='', obs_value=None,
                  aerosol_subtypes=None, figsize=(8, 4), save_path=None):
    """Plot PAP coefficient bar chart (Fig.4 style).

    Args:
        pap_values: dict of {process_name: PAP value}
        title: subplot title (e.g., "EH13" or "EH22")
        obs_value: observed total PAP for hatched bar overlay
        aerosol_subtypes: dict of {subtype: PAP} for inset plot
        figsize: figure size
        save_path: output file path
    """
    setup_style()

    # Standard order matching Fig.4
    order = ['SR', 'AL', 'WV', 'O3', 'CO2', 'CH4', 'CLD', 'AER', 'SURF', 'ATMD', 'SUM']

    # Map from internal names to display names
    name_map = {
        'solar': 'SR', 'albedo': 'AL', 'q': 'WV', 'o3': 'O3',
        'co2': 'CO2', 'ch4': 'CH4', 'cloud': 'CLD', 'aerosol': 'AER',
        'ts': 'SURF', 'atd': 'ATMD', 'warm': 'SUM',
    }

    # Reorder
    labels = []
    values = []
    for display in order:
        # Find matching key
        found = False
        for key, disp in name_map.items():
            if disp == display and key in pap_values:
                labels.append(display)
                values.append(pap_values[key])
                found = True
                break
        if not found:
            labels.append(display)
            values.append(0.0)

    values = np.array(values)

    fig, ax = plt.subplots(figsize=figsize)

    colors = ['#4393c3' if v < 0 else '#d6604d' for v in values]
    bars = ax.bar(labels, values, color=colors, edgecolor='black', linewidth=0.5)

    # Hatched bar for observed SUM
    if obs_value is not None:
        sum_idx = labels.index('SUM') if 'SUM' in labels else -1
        if sum_idx >= 0:
            ax.bar([labels[sum_idx]], [obs_value], color='none',
                   edgecolor='black', linewidth=1, hatch='///')

    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.set_ylabel('PAP (K)', fontsize=FONTSIZE_LABEL)
    ax.set_title(title, fontsize=FONTSIZE_PANEL, fontweight='bold')
    ax.tick_params(axis='x', rotation=0)

    # Inset for aerosol subtypes
    if aerosol_subtypes:
        inset = ax.inset_axes([0.02, 0.65, 0.25, 0.3])
        aer_labels = list(aerosol_subtypes.keys())
        aer_vals = list(aerosol_subtypes.values())
        aer_colors = ['#d6604d' if v > 0 else '#4393c3' for v in aer_vals]
        inset.bar(aer_labels, aer_vals, color=aer_colors,
                  edgecolor='black', linewidth=0.3)
        inset.axhline(y=0, color='black', linewidth=0.3)
        inset.tick_params(labelsize=6)
        inset.set_ylabel('PAP', fontsize=6)

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")

    return fig


def plot_pap_comparison(pap_eh13, pap_eh22, obs_eh13=None, obs_eh22=None,
                        aer_eh13=None, aer_eh22=None, save_path=None):
    """Plot side-by-side PAP bar charts for two events (Fig.4 style)."""
    setup_style()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 7))

    for ax, pap, title, obs, aer in [
        (ax1, pap_eh13, '(a) EH13', obs_eh13, aer_eh13),
        (ax2, pap_eh22, '(b) EH22', obs_eh22, aer_eh22),
    ]:
        plt.sca(ax)
        plot_pap_bars(pap, title=title, obs_value=obs,
                      aerosol_subtypes=aer)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")
    return fig
