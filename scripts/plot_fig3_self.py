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

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SELF_DIR = os.path.join(PROJECT_ROOT, "cfram_output")
PAPER_BASE = os.path.join(PROJECT_ROOT, "paper_data", "cfram_out")
OUTDIR = os.path.join(PROJECT_ROOT, "figures")
CASES_PAPER = {'EH13': 'case_eh13_c20250102', 'EH22': 'case_eh22_c20250118'}

LEVELS_DT = [-15, -7, -2, -0.5, 0, 0.5, 2, 7, 15]
CMAP = plt.cm.RdBu_r
KEY_REGIONS = {'EH13': [110, 122, 28, 36], 'EH22': [103, 122, 27, 35]}
NLEV = 37

PLOT_ROWS = [
    ('wv',     'Water Vapor'),
    ('cld',    'Cloud'),
    ('aer',    'Aerosol'),
    ('sfc',    'Surface Process'),
    ('atmdyn', 'Atm. Dynamics'),
    ('total',  'Total'),
]


def load_case(case_label):
    """Load self-computed results from NetCDF."""
    cl = case_label.lower()
    self_file = os.path.join(SELF_DIR, 'cfram_%s_python.nc' % cl)
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
    os.makedirs(OUTDIR, exist_ok=True)
    norm = mcolors.BoundaryNorm(LEVELS_DT, CMAP.N, clip=True)
    nrows = len(PLOT_ROWS)

    fig, axes = plt.subplots(nrows, 2, figsize=(10, 3.0 * nrows),
                             subplot_kw={'projection': ccrs.PlateCarree()})

    labels_abc = 'abcdefghijklmnopqrstuvwxyz'
    for col, cl in enumerate(['EH13', 'EH22']):
        lats, lons, data = load_case(cl)
        for row, (key, label) in enumerate(PLOT_ROWS):
            ax = axes[row, col]
            field = data.get(key, np.full((len(lats), len(lons)), np.nan))
            field = np.clip(np.nan_to_num(field, nan=0), -20, 20)
            panel_label = labels_abc[row * 2 + col]
            title = "(%s) %s %s" % (panel_label, label, cl)
            plot_panel(ax, lons, lats, field, title, norm, CMAP, KEY_REGIONS[cl])

    fig.subplots_adjust(bottom=0.06, top=0.96, left=0.06, right=0.94, hspace=0.25, wspace=0.15)
    cbar_ax = fig.add_axes([0.15, 0.02, 0.7, 0.015])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cbar_ax,
                      orientation='horizontal', ticks=LEVELS_DT)
    cb.set_label('Partial Temperature (K)', fontsize=10)
    outpath = os.path.join(OUTDIR, 'fig3_self_computed.png')
    fig.savefig(outpath, dpi=200, bbox_inches='tight')
    plt.close()
    print("Saved: %s" % outpath)

    # Validation vs paper
    print("\n=== Validation vs paper_data ===")
    for cl in ['EH13', 'EH22']:
        lats, lons, ds = load_case(cl)
        paper_dir = os.path.join(PAPER_BASE, CASES_PAPER[cl])
        files = os.listdir(paper_dir)
        nc = Dataset(os.path.join(paper_dir, [f for f in files if 'partial_t' in f][0]))
        nc2 = Dataset(os.path.join(paper_dir, [f for f in files if 'total_t' in f][0]))
        skip = {'time','lat','lon','lev','bounds_time','bounds_latitude','bounds_longitude','bounds_level'}

        paper = {}
        paper['wv'] = np.array(nc.variables['wv'][0, 0, :, :], dtype=np.float64)
        paper['cld'] = np.array(nc.variables['cld'][0, 0, :, :], dtype=np.float64)
        paper['aer'] = sum(np.array(nc.variables[v][0, 0, :, :], dtype=np.float64)
                          for v in ['bc','oc','sulf','dust','seas'])
        paper['sfc'] = sum(np.array(nc.variables[v][0, 0, :, :], dtype=np.float64)
                          for v in ['sfcdyn','lhflx','shflx'])
        paper['atmdyn'] = np.array(nc.variables['atmdyn'][0, 0, :, :], dtype=np.float64)
        tv = [v for v in nc2.variables if v not in skip][0]
        paper['total'] = np.array(nc2.variables[tv][0, 0, :, :], dtype=np.float64)
        nc.close(); nc2.close()

        print("\n%s:" % cl)
        for key in ['wv', 'cld', 'aer', 'sfc', 'atmdyn', 'total']:
            s = ds[key]
            p = paper[key]
            mask = np.isfinite(s) & np.isfinite(p) & (np.abs(p) < 900) & (np.abs(s) < 900)
            if mask.sum() == 0: continue
            diff = s[mask] - p[mask]
            corr = np.corrcoef(s[mask], p[mask])[0, 1]
            print("  %-8s: RMSE=%.3f K, bias=%.3f K, corr=%.4f" %
                  (key, np.sqrt(np.mean(diff**2)), np.mean(diff), corr))


if __name__ == '__main__':
    main()
