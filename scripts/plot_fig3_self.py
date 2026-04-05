#!/usr/bin/env python3
"""Plot Fig.3 from self-computed CFRAM results.

All terms self-computed:
- WV, Cloud, Aerosol: RRTMG radiative perturbation
- Surface Process: sfcdyn+lhflx+shflx via Planck matrix × paper forcing
- Atm. Dynamics: residual = observed - radiative - surface
- Total: observed dT = T_warm - T_base

Run on hqlx204:
    python3 scripts/plot_fig3_self.py
"""
import os, sys
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

import argparse

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case, defaults, get_plev, PROJECT_ROOT

_defaults = defaults()
LEVELS_DT = _defaults['plotting']['levels_dt']
CMAP_NAME = _defaults['plotting']['cmap']
CMAP = plt.cm.get_cmap(CMAP_NAME)
NLEV = len(get_plev())

PLOT_ROWS = [
    ('wv',     'Water Vapor'),
    ('cld',    'Cloud'),
    ('aer',    'Aerosol'),
    ('sfc',    'Surface Process'),
    ('atmdyn', 'Atm. Dynamics'),
    ('total',  'Total'),
]


def load_case_data(case_name):
    """Load self-computed results from NetCDF."""
    case_cfg = load_case(case_name)
    self_file = os.path.join(case_cfg['_output_dir'], 'cfram_result.nc')
    nc = Dataset(self_file)
    lats = np.array(nc.variables['lat'][:])
    lons = np.array(nc.variables['lon'][:])
    sfc = -1  # surface = last index

    def get(vname):
        arr = np.array(nc.variables[vname][sfc, :, :], dtype=np.float64)
        return np.where(np.abs(arr) > 900, np.nan, arr)

    data = {}
    data['wv'] = get('dT_q')
    data['cld'] = get('dT_cloud')
    data['aer'] = get('dT_aerosol')
    data['total'] = get('dT_observed')
    data['atmdyn'] = get('dT_atmdyn')

    # Surface process = sfcdyn + lhflx + shflx (all via Planck matrix)
    sfc_sum = np.zeros_like(data['total'])
    for t in ['sfcdyn', 'lhflx', 'shflx']:
        if 'dT_' + t in nc.variables:
            arr = get('dT_' + t)
            arr = np.nan_to_num(arr, nan=0.0)
            sfc_sum += arr
    data['sfc'] = sfc_sum

    nc.close()
    return lats, lons, data


def plot_panel(ax, lons, lats, field, title, norm, cmap, key_region=None):
    ax.set_extent([95, 125, 18, 42], crs=ccrs.PlateCarree())
    cf = ax.contourf(lons, lats, field, levels=LEVELS_DT, cmap=cmap, norm=norm,
                     transform=ccrs.PlateCarree(), extend='both')
    ax.coastlines(resolution='50m', linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle='-', edgecolor='gray')
    if key_region:
        lon0, lon1, lat0, lat1 = key_region
        ax.plot([lon0,lon1,lon1,lon0,lon0], [lat0,lat0,lat1,lat1,lat0],
                'b-', linewidth=1.5, transform=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False; gl.right_labels = False
    gl.xlocator = plt.FixedLocator([95,100,105,110,115,120,125])
    gl.ylocator = plt.FixedLocator([20,25,30,35,40])
    gl.xformatter = LongitudeFormatter(); gl.yformatter = LatitudeFormatter()
    gl.xlabel_style = {'size': 7}; gl.ylabel_style = {'size': 7}
    ax.set_title(title, fontsize=9, fontweight='bold')
    return cf


def main():
    parser = argparse.ArgumentParser(description='Plot CFRAM decomposition maps')
    parser.add_argument('--case', nargs='+', required=True,
                        help='Case name(s), e.g., --case eh13 eh22')
    args = parser.parse_args()

    case_names = args.case
    ncols = len(case_names)
    norm = mcolors.BoundaryNorm(LEVELS_DT, CMAP.N, clip=True)
    nrows = len(PLOT_ROWS)

    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.0 * nrows),
                             subplot_kw={'projection': ccrs.PlateCarree()},
                             squeeze=False)

    labels_abc = 'abcdefghijklmnopqrstuvwxyz'
    for col, case_name in enumerate(case_names):
        case_cfg = load_case(case_name)
        lats, lons, data = load_case_data(case_name)
        key_region = case_cfg.get('plot', {}).get('key_region')
        kr = None
        if key_region:
            kr = key_region['lon'] + key_region['lat']  # [lon0,lon1,lat0,lat1]
        case_label = case_cfg.get('case_name', case_name.upper())

        for row, (key, label) in enumerate(PLOT_ROWS):
            ax = axes[row, col]
            field = data.get(key, np.full((len(lats), len(lons)), np.nan))
            field = np.clip(np.nan_to_num(field, nan=0), -20, 20)
            panel_label = labels_abc[row * ncols + col]
            title = "(%s) %s %s" % (panel_label, label, case_label)
            plot_panel(ax, lons, lats, field, title, norm, CMAP, kr)

    fig.subplots_adjust(bottom=0.06, top=0.96, left=0.06, right=0.94, hspace=0.25, wspace=0.15)
    cbar_ax = fig.add_axes([0.15, 0.02, 0.7, 0.015])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cbar_ax,
                      orientation='horizontal', ticks=LEVELS_DT)
    cb.set_label('Partial Temperature (K)', fontsize=10)

    # Save to first case's figures dir
    first_cfg = load_case(case_names[0])
    outdir = first_cfg['_figures_dir']
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, 'fig3_decomposition.png')
    fig.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print("Saved: %s" % outpath)


if __name__ == '__main__':
    main()
