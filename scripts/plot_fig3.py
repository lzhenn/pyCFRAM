#!/usr/bin/env python3
"""Reproduce Fig.3: Spatial decomposition maps for EH13 and EH22.

Surface-level partial_t from paper_data.
Layout: N rows x 2 columns (EH13 left, EH22 right).

Run on hqlx204:
    python3 scripts/plot_fig3.py
"""
import os
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

# ---- Config ----
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER_BASE = os.path.join(PROJECT_ROOT, "paper_data", "cfram_out")
CASES = {
    'EH13': 'case_eh13_c20250102',
    'EH22': 'case_eh22_c20250118',
}
OUTDIR = os.path.join(PROJECT_ROOT, "figures")

# Wu et al. 2025 style
LEVELS_DT = [-15, -7, -2, -0.5, 0, 0.5, 2, 7, 15]
CMAP = plt.cm.RdBu_r

# Key regions (EH13 and EH22 boxes from paper)
KEY_REGIONS = {
    'EH13': [110, 122, 28, 36],  # lon_min, lon_max, lat_min, lat_max
    'EH22': [103, 122, 27, 35],
}

# Terms to plot (row order matching paper Fig.3)
PLOT_TERMS = [
    ('wv',     'Water Vapor'),
    ('cld',    'Cloud'),
    ('aer',    'Aerosol'),        # sum of bc+oc+sulf+dust+seas
    ('sfc',    'Surface Process'), # sfcdyn + lhflx + shflx
    ('atmdyn', 'Atm. Dynamics'),   # atmdyn only
    ('total',  'Total'),
]


def load_case(case_dir):
    """Load partial_t surface values and total_t for one case."""
    files = os.listdir(case_dir)
    # partial_t
    nc = Dataset(os.path.join(case_dir, [f for f in files if 'partial_t' in f][0]))
    lats = np.array(nc.variables['lat'][:])
    lons = np.array(nc.variables['lon'][:])
    levs = np.array(nc.variables['lev'][:])

    # Surface is lev index 0 (1013 hPa, surface->TOA order)
    sfc_idx = 0
    data = {}
    skip = {'time', 'lat', 'lon', 'lev', 'bounds_time', 'bounds_latitude',
            'bounds_longitude', 'bounds_level'}
    for v in nc.variables:
        if v in skip:
            continue
        arr = np.array(nc.variables[v][0, sfc_idx, :, :], dtype=np.float64)
        # Mask fill values
        arr = np.where(np.abs(arr) > 900, np.nan, arr)
        data[v] = arr
    nc.close()

    # total_t
    nc = Dataset(os.path.join(case_dir, [f for f in files if 'total_t' in f][0]))
    total_var = [v for v in nc.variables if v not in skip][0]
    data['total'] = np.where(
        np.abs(nc.variables[total_var][0, sfc_idx, :, :]) > 900,
        np.nan,
        np.array(nc.variables[total_var][0, sfc_idx, :, :], dtype=np.float64))
    nc.close()

    # Compute aerosol sum
    aer_vars = ['bc', 'oc', 'sulf', 'dust', 'seas']
    data['aer'] = sum(data.get(v, np.zeros_like(lats[:, None] * lons[None, :]))
                      for v in aer_vars if v in data)

    # Surface process = sfcdyn + lhflx + shflx
    zero = np.zeros_like(lats[:, None] * lons[None, :])
    data['sfc'] = (data.get('sfcdyn', zero) +
                   data.get('lhflx', zero) +
                   data.get('shflx', zero))

    return lats, lons, data


def plot_panel(ax, lons, lats, data, title, norm, cmap, key_region=None):
    """Plot one panel with filled contours on a map."""
    ax.set_extent([95, 125, 18, 42], crs=ccrs.PlateCarree())
    # Filled contours
    cf = ax.contourf(lons, lats, data, levels=LEVELS_DT, cmap=cmap, norm=norm,
                     transform=ccrs.PlateCarree(), extend='both')
    # Coastlines and borders
    ax.coastlines(resolution='50m', linewidth=0.5, color='k')
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle='-', edgecolor='gray')

    # Key region box
    if key_region is not None:
        lon0, lon1, lat0, lat1 = key_region
        ax.plot([lon0, lon1, lon1, lon0, lon0],
                [lat0, lat0, lat1, lat1, lat0],
                'b-', linewidth=1.5, transform=ccrs.PlateCarree())

    # Gridlines
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='gray',
                      alpha=0.5, linestyle='--')
    gl.top_labels = False
    gl.right_labels = False
    gl.xlocator = plt.FixedLocator([95, 100, 105, 110, 115, 120, 125])
    gl.ylocator = plt.FixedLocator([20, 25, 30, 35, 40])
    gl.xformatter = LongitudeFormatter()
    gl.yformatter = LatitudeFormatter()
    gl.xlabel_style = {'size': 7}
    gl.ylabel_style = {'size': 7}

    ax.set_title(title, fontsize=9, fontweight='bold')
    return cf


def main():
    os.makedirs(OUTDIR, exist_ok=True)

    # Load data
    case_data = {}
    for label, case_dir_name in CASES.items():
        case_path = os.path.join(PAPER_BASE, case_dir_name)
        lats, lons, data = load_case(case_path)
        case_data[label] = (lats, lons, data)

    # Color normalization
    norm = mcolors.BoundaryNorm(LEVELS_DT, CMAP.N, clip=True)

    nrows = len(PLOT_TERMS)
    fig, axes = plt.subplots(nrows, 2, figsize=(10, 3.0 * nrows),
                             subplot_kw={'projection': ccrs.PlateCarree()})

    labels = 'abcdefghijklmnopqrstuvwxyz'
    for row, (term_key, term_label) in enumerate(PLOT_TERMS):
        for col, case_label in enumerate(['EH13', 'EH22']):
            ax = axes[row, col]
            lats, lons, data = case_data[case_label]
            field = data.get(term_key, np.full((len(lats), len(lons)), np.nan))

            panel_label = labels[row * 2 + col]
            title = "(%s) %s %s" % (panel_label, term_label, case_label)
            key_region = KEY_REGIONS[case_label]

            cf = plot_panel(ax, lons, lats, field, title, norm, CMAP,
                           key_region=key_region)

    # Colorbar
    fig.subplots_adjust(bottom=0.06, top=0.96, left=0.06, right=0.94,
                        hspace=0.25, wspace=0.15)
    cbar_ax = fig.add_axes([0.15, 0.02, 0.7, 0.015])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP),
                      cax=cbar_ax, orientation='horizontal',
                      ticks=LEVELS_DT)
    cb.set_label('Partial Temperature (K)', fontsize=10)
    cb.ax.tick_params(labelsize=8)

    outpath = os.path.join(OUTDIR, 'fig3_decomposition.png')
    fig.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print("Saved: %s" % outpath)

    # Also print domain-mean values for verification
    print("\n--- Domain-mean surface dT (K) ---")
    print("%-15s %10s %10s" % ("Term", "EH13", "EH22"))
    for term_key, term_label in PLOT_TERMS:
        vals = []
        for cl in ['EH13', 'EH22']:
            _, _, d = case_data[cl]
            f = d.get(term_key, np.full(1, np.nan))
            vals.append(np.nanmean(f))
        print("%-15s %10.3f %10.3f" % (term_label, vals[0], vals[1]))


if __name__ == '__main__':
    main()
