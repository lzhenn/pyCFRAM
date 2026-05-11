#!/usr/bin/env python3
"""Side-by-side 13-panel comparison: cesm2_4xco2 (collaborator, 37 lev,
GOCART aerosols, interpolated cl) vs cesm2_4xco2_official (CMIP6 raw,
19 lev plev, hybrid→plev cl, no aerosols).

Quantifies whether the collaborator's preprocessing introduced bias by
interpolating ta/hus from 19→37 plev and constructing aerosol fields.

Usage:
    python3 scripts/compare_cesm2_official_vs_collaborator.py
"""
import os
import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

LEVELS = [-15.0, -7.0, -2.0, -0.5, 0.0, 0.5, 2.0, 7.0, 15.0]
CMAP = plt.cm.RdBu_r

# Same panel order as plot_13panel_global.py
PANELS = [
    ('AL',   'dT_albedo',   'a'),
    ('WV',   'dT_q',        'b'),
    ('CLD',  'dT_cloud',    'c'),
    ('CLDS', 'dT_cloud_sw', 'd'),
    ('CLDL', 'dT_cloud_lw', 'e'),
    ('CO2',  'dT_co2',      'f'),
    ('O3',   'dT_o3',       'g'),
    ('SR',   'dT_solar',    'h'),
    ('DYN',  'dT_dry',      'i'),
    ('ATM',  'dT_atmdyn',   'j'),
    ('OCH',  'dT_ocndyn',   'k'),
    ('SH',   'dT_shflx',    'l'),
    ('LH',   'dT_lhflx',    'm'),
]


def load_surface(nc_path):
    """Load surface row from cfram_result.nc, return dict {var: (lat, lon)}."""
    nc = Dataset(nc_path)
    lats = np.array(nc.variables['lat'][:])
    lons = np.array(nc.variables['lon'][:])
    out = {'lat': lats, 'lon': lons}
    for _, var, _ in PANELS:
        if var in nc.variables:
            # Surface row is index NLEV (last in CFRAM array). Post lev-fix,
            # NetCDF lev order is sfc→TOA + ground (1013) at end.
            # NetCDF index = -1 = ground, but NLEV row corresponds to surface temperature
            # which in the saved layout (atm[NLEV-1::-1] + sfc) is at index NLEV (last).
            # Equivalently: arr[-1] always returns the surface row.
            arr = np.array(nc.variables[var][-1, :, :], dtype=np.float64)
            arr = np.where(np.abs(arr) > 900, np.nan, arr)
        else:
            arr = np.full((len(lats), len(lons)), np.nan)
        out[var] = arr
    nc.close()
    return out


def plot_compare(off, col, out_path):
    """Side-by-side: official (left col), collaborator (right col)."""
    fig = plt.figure(figsize=(28, 28))
    norm = mcolors.BoundaryNorm(LEVELS, ncolors=CMAP.N, extend='both')
    proj = ccrs.PlateCarree(central_longitude=0)
    nrows = (len(PANELS) + 1) // 2  # 7

    for idx, (label, var, letter) in enumerate(PANELS):
        for col_idx, (data, src_label) in enumerate([(off, 'OFFICIAL (CMIP6 19lev)'),
                                                     (col, 'COLLAB (37lev interp)')]):
            ax = fig.add_subplot(nrows, 4, idx * 2 + col_idx + 1, projection=proj)
            field = data[var]
            if not np.all(np.isnan(field)):
                cf = ax.contourf(data['lon'], data['lat'], field,
                                 levels=LEVELS, cmap=CMAP, norm=norm,
                                 extend='both', transform=ccrs.PlateCarree())
            ax.coastlines(resolution='110m', linewidth=0.4)
            ax.set_global()
            ax.set_title(f'({letter}) {label}  [{src_label.split()[0]}]', fontsize=9, loc='left')
            dm = np.nanmean(field)
            ax.text(0.02, 0.02, f'mean={dm:+.2f}K', transform=ax.transAxes, fontsize=7,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

    cbar_ax = fig.add_axes([0.15, 0.04, 0.7, 0.012])
    fig.colorbar(cf, cax=cbar_ax, orientation='horizontal',
                 ticks=LEVELS, extend='both',
                 label='Partial Temperature dT (K)')
    fig.suptitle('CESM2 4×CO2 — OFFICIAL (CMIP6 raw, 19 plev) vs COLLABORATOR (37 plev interp)',
                 fontsize=13, y=0.96)
    plt.tight_layout(rect=[0, 0.06, 1, 0.94])
    plt.savefig(out_path, dpi=130, bbox_inches='tight')
    print('Saved:', out_path)


def main():
    off_nc = os.path.join(ROOT, 'cases/cesm2_4xco2_official/output/cfram_result.nc')
    col_nc = os.path.join(ROOT, 'cases/cesm2_4xco2/output/cfram_result.nc')
    if not os.path.exists(off_nc):
        sys.exit('Missing %s — run cesm2_4xco2_official first' % off_nc)
    if not os.path.exists(col_nc):
        sys.exit('Missing %s' % col_nc)

    print('Loading official (19lev)...')
    off = load_surface(off_nc)
    print('Loading collaborator (37lev)...')
    col = load_surface(col_nc)

    # Stats summary
    print('\n--- Surface global mean dT (K) ---')
    print('  %-8s | %-12s | %-12s | %-12s' % ('var', 'official', 'collab', 'diff'))
    for label, var, _ in PANELS:
        ov = np.nanmean(off[var])
        cv = np.nanmean(col[var])
        print('  %-8s | %+11.4f | %+11.4f | %+11.4f' % (label, ov, cv, ov - cv))

    out_dir = os.path.join(ROOT, 'cases/cesm2_4xco2_official/figures')
    os.makedirs(out_dir, exist_ok=True)
    plot_compare(off, col, os.path.join(out_dir, 'fig_13panel_official_vs_collab.png'))


if __name__ == '__main__':
    main()
