#!/usr/bin/env python3
"""Reproduce Fig.4: PAP (Pattern-Amplitude Projection) bar charts.

PAP_i = <dT_total> * <dT_total * dT_i> / <dT_total^2>
where <> denotes area-weighted average over key region.

Run on hqlx204:
    python3 scripts/plot_fig4.py
"""
import os
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER_BASE = os.path.join(PROJECT_ROOT, "paper_data", "cfram_out")
CASES = {
    'EH13': ('case_eh13_c20250102', [110, 122, 28, 36]),
    'EH22': ('case_eh22_c20250118', [103, 122, 27, 35]),
}
OUTDIR = os.path.join(PROJECT_ROOT, "figures")

# Terms for PAP (order matching paper Fig.4)
PAP_TERMS = [
    ('wv',      'WV',       '#2166AC'),
    ('cld',     'Cloud',    '#D6604D'),
    ('aer',     'Aerosol',  '#F4A582'),
    ('sfcdyn',  'SfcDyn',   '#92C5DE'),
    ('dyn',     'AtmDyn',   '#4393C3'),
    ('co2',     'CO2',      '#FDDBC7'),
    ('ozone',   'O3',       '#B2182B'),
    ('solar',   'Solar',    '#FEE08B'),
    ('albedo',  'Albedo',   '#D9EF8B'),
    ('ch4',     'CH4',      '#E0E0E0'),
    ('lhflx',   'LH',       '#66BD63'),
    ('shflx',   'SH',       '#1A9850'),
    ('atmdyn',  'DynRes',   '#878787'),
]

# Aerosol sub-terms for inset
AER_SUB = [
    ('bc',   'BC',    '#333333'),
    ('oc',   'OC',    '#A6611A'),
    ('sulf', 'Sulf',  '#DFC27D'),
    ('dust', 'Dust',  '#BF812D'),
    ('seas', 'SS',    '#80CDC1'),
]


def load_surface(case_dir):
    """Load surface partial_t and total_t."""
    files = os.listdir(case_dir)
    nc = Dataset(os.path.join(case_dir, [f for f in files if 'partial_t' in f][0]))
    lats = np.array(nc.variables['lat'][:])
    lons = np.array(nc.variables['lon'][:])
    sfc_idx = 0  # surface = first level (1013 hPa)
    data = {}
    skip = {'time', 'lat', 'lon', 'lev', 'bounds_time', 'bounds_latitude',
            'bounds_longitude', 'bounds_level'}
    for v in nc.variables:
        if v in skip: continue
        arr = np.array(nc.variables[v][0, sfc_idx, :, :], dtype=np.float64)
        data[v] = np.where(np.abs(arr) > 900, np.nan, arr)
    nc.close()
    nc = Dataset(os.path.join(case_dir, [f for f in files if 'total_t' in f][0]))
    tv = [v for v in nc.variables if v not in skip][0]
    data['total'] = np.where(
        np.abs(nc.variables[tv][0, sfc_idx, :, :]) > 900, np.nan,
        np.array(nc.variables[tv][0, sfc_idx, :, :], dtype=np.float64))
    nc.close()
    # Aerosol sum
    data['aer'] = sum(data.get(v, 0) for v in ['bc', 'oc', 'sulf', 'dust', 'seas'])
    return lats, lons, data


def compute_pap(dT_total, dT_i, weights):
    """Compute PAP coefficient with area weights.
    PAP_i = <dT> * <dT * dT_i> / <dT^2>
    """
    mask = np.isfinite(dT_total) & np.isfinite(dT_i) & np.isfinite(weights)
    if mask.sum() == 0:
        return 0.0
    w = weights[mask]
    dt = dT_total[mask]
    di = dT_i[mask]
    mean_dt = np.average(dt, weights=w)
    mean_dt_di = np.average(dt * di, weights=w)
    mean_dt2 = np.average(dt ** 2, weights=w)
    if mean_dt2 == 0:
        return 0.0
    return mean_dt * mean_dt_di / mean_dt2


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for col, (case_label, (case_dir_name, region)) in enumerate(CASES.items()):
        case_path = os.path.join(PAPER_BASE, case_dir_name)
        lats, lons, data = load_surface(case_path)

        # Key region mask
        lon0, lon1, lat0, lat1 = region
        lon2d, lat2d = np.meshgrid(lons, lats)
        region_mask = (lon2d >= lon0) & (lon2d <= lon1) & (lat2d >= lat0) & (lat2d <= lat1)
        # Area weights (cos latitude)
        cos_weights = np.cos(np.deg2rad(lat2d))
        weights = np.where(region_mask, cos_weights, np.nan)

        dT_total = data['total']

        # Compute PAP for each term
        pap_vals = []
        pap_labels = []
        pap_colors = []
        for term_key, term_label, color in PAP_TERMS:
            if term_key in data:
                pap = compute_pap(dT_total, data[term_key], weights)
            else:
                pap = 0.0
            pap_vals.append(pap)
            pap_labels.append(term_label)
            pap_colors.append(color)

        # Main bar chart
        ax = axes[col]
        x = np.arange(len(pap_vals))
        bars = ax.bar(x, pap_vals, color=pap_colors, edgecolor='k', linewidth=0.5)
        ax.set_xticks(x)
        ax.set_xticklabels(pap_labels, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel('PAP (K)', fontsize=10)
        ax.set_title('%s Key Region %s' % (case_label, region), fontsize=11, fontweight='bold')
        ax.axhline(y=0, color='k', linewidth=0.5)
        ax.set_ylim(-3, 6)

        # PAP sum annotation
        pap_sum = sum(pap_vals)
        mean_dt = np.nanmean(dT_total[region_mask])
        ax.text(0.02, 0.95, 'Sum=%.2f K\n<ΔT>=%.2f K' % (pap_sum, mean_dt),
                transform=ax.transAxes, fontsize=8, va='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

        # Inset: aerosol sub-terms
        ax_inset = ax.inset_axes([0.62, 0.55, 0.35, 0.40])
        aer_vals = []
        aer_labels = []
        aer_colors = []
        for sv, sl, sc in AER_SUB:
            if sv in data:
                pap = compute_pap(dT_total, data[sv], weights)
            else:
                pap = 0.0
            aer_vals.append(pap)
            aer_labels.append(sl)
            aer_colors.append(sc)
        ax_inset.bar(range(len(aer_vals)), aer_vals, color=aer_colors,
                     edgecolor='k', linewidth=0.4)
        ax_inset.set_xticks(range(len(aer_vals)))
        ax_inset.set_xticklabels(aer_labels, fontsize=6, rotation=30)
        ax_inset.set_ylabel('PAP (K)', fontsize=7)
        ax_inset.set_title('Aerosol', fontsize=7)
        ax_inset.axhline(y=0, color='k', linewidth=0.3)
        ax_inset.tick_params(labelsize=6)

        # Print values
        print("\n%s PAP values:" % case_label)
        for i, (lbl, val) in enumerate(zip(pap_labels, pap_vals)):
            print("  %-10s: %7.3f K" % (lbl, val))
        print("  Sum:        %7.3f K" % pap_sum)
        print("  <dT>:       %7.3f K" % mean_dt)

    fig.tight_layout()
    outpath = os.path.join(OUTDIR, 'fig4_pap.png')
    fig.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print("\nSaved: %s" % outpath)


if __name__ == '__main__':
    main()
