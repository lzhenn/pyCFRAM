#!/usr/bin/env python3
"""Reproduce Fig.5: Aerosol species decomposition maps.

Two versions:
- paper_data: individual species from partial_t.nc (bc, oc, sulf, seas, dust)
- self-computed: total aerosol from cfram_result.nc (single combined term)

Usage:
    python3 scripts/plot_fig5.py --case eh13 eh22
"""
import os, sys, argparse
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case, defaults, PROJECT_ROOT

LEVELS_AER = [-3, -1.5, -0.5, 0, 0.5, 1.5, 3]
CMAP = plt.cm.RdBu_r

AER_ROWS = [
    ('bc',    'Black Carbon'),
    ('oc',    'Organic Carbon'),
    ('sulf',  'Sulfate'),
    ('seas',  'Sea Salt'),
    ('dust',  'Dust'),
    ('total', 'Total Aerosol'),
]

# paper_data case directory names (for validation plots)
PAPER_CASES = {
    'eh13': 'case_eh13_c20250102',
    'eh22': 'case_eh22_c20250118',
}


def load_paper_aerosol(case_name):
    """Load per-species aerosol dT from paper_data partial_t.nc."""
    paper_dir_name = PAPER_CASES.get(case_name)
    if not paper_dir_name:
        return None, None, None
    paper_dir = os.path.join(PROJECT_ROOT, 'paper_data', 'cfram_out', paper_dir_name)
    if not os.path.isdir(paper_dir):
        return None, None, None

    files = os.listdir(paper_dir)
    pt_file = [f for f in files if 'partial_t' in f]
    if not pt_file:
        return None, None, None

    nc = Dataset(os.path.join(paper_dir, pt_file[0]))
    lats = np.array(nc.variables['lat'][:])
    lons = np.array(nc.variables['lon'][:])
    sfc = 0  # surface level in paper_data (surface→TOA order)

    data = {}
    for species in ['bc', 'oc', 'sulf', 'seas', 'dust']:
        if species in nc.variables:
            arr = np.array(nc.variables[species][0, sfc, :, :], dtype=np.float64)
            data[species] = np.where(np.abs(arr) > 900, np.nan, arr)

    # Total = sum of all species
    data['total'] = sum(data.get(s, 0) for s in ['bc', 'oc', 'sulf', 'seas', 'dust'])
    nc.close()
    return lats, lons, data


def load_self_aerosol(case_name):
    """Load total aerosol dT from self-computed cfram_result.nc."""
    cfg = load_case(case_name)
    result_file = os.path.join(cfg['_output_dir'], 'cfram_result.nc')
    if not os.path.exists(result_file):
        return None, None, None

    nc = Dataset(result_file)
    lats = np.array(nc.variables['lat'][:])
    lons = np.array(nc.variables['lon'][:])
    sfc = -1  # surface = last index

    data = {}
    if 'dT_aerosol' in nc.variables:
        arr = np.array(nc.variables['dT_aerosol'][sfc, :, :], dtype=np.float64)
        data['total'] = np.where(np.abs(arr) > 900, np.nan, arr)
    nc.close()
    return lats, lons, data


def plot_panel(ax, lons, lats, field, title, norm, cmap, key_region=None):
    ax.set_extent([95, 125, 18, 42], crs=ccrs.PlateCarree())
    cf = ax.contourf(lons, lats, field, levels=LEVELS_AER, cmap=cmap, norm=norm,
                     transform=ccrs.PlateCarree(), extend='both')
    ax.coastlines(resolution='50m', linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, linestyle='-', edgecolor='gray')
    if key_region:
        lon0, lon1, lat0, lat1 = key_region
        ax.plot([lon0, lon1, lon1, lon0, lon0], [lat0, lat0, lat1, lat1, lat0],
                'b-', linewidth=1.5, transform=ccrs.PlateCarree())
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False; gl.right_labels = False
    gl.xlocator = plt.FixedLocator([95, 100, 105, 110, 115, 120, 125])
    gl.ylocator = plt.FixedLocator([20, 25, 30, 35, 40])
    gl.xformatter = LongitudeFormatter(); gl.yformatter = LatitudeFormatter()
    gl.xlabel_style = {'size': 7}; gl.ylabel_style = {'size': 7}
    ax.set_title(title, fontsize=9, fontweight='bold')
    return cf


def main():
    parser = argparse.ArgumentParser(description='Plot Fig.5: aerosol species decomposition')
    parser.add_argument('--case', nargs='+', default=['eh13', 'eh22'],
                        help='Case name(s)')
    args = parser.parse_args()

    case_names = args.case
    ncols = len(case_names)
    nrows = len(AER_ROWS)
    norm = mcolors.BoundaryNorm(LEVELS_AER, CMAP.N, clip=True)

    # === Paper data version ===
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.0 * nrows),
                             subplot_kw={'projection': ccrs.PlateCarree()},
                             squeeze=False)
    labels_abc = 'abcdefghijklmnopqrstuvwxyz'
    has_paper = False

    for col, cn in enumerate(case_names):
        cfg = load_case(cn)
        kr = cfg.get('plot', {}).get('key_region')
        kr_list = (kr['lon'] + kr['lat']) if kr else None
        case_label = cfg.get('case_name', cn.upper())

        lats, lons, pdata = load_paper_aerosol(cn)
        if pdata is not None:
            has_paper = True
            for row, (key, label) in enumerate(AER_ROWS):
                ax = axes[row, col]
                field = pdata.get(key, np.full((len(lats), len(lons)), np.nan))
                field = np.clip(np.nan_to_num(field, nan=0), -5, 5)
                pl = labels_abc[row * ncols + col]
                plot_panel(ax, lons, lats, field, "(%s) %s %s" % (pl, label, case_label),
                          norm, CMAP, kr_list)

    if has_paper:
        fig.subplots_adjust(bottom=0.06, top=0.96, left=0.06, right=0.94, hspace=0.25, wspace=0.15)
        cbar_ax = fig.add_axes([0.15, 0.02, 0.7, 0.015])
        fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cbar_ax,
                     orientation='horizontal', ticks=LEVELS_AER,
                     label='Partial Temperature (K)')

        first_cfg = load_case(case_names[0])
        outdir = first_cfg['_figures_dir']
        os.makedirs(outdir, exist_ok=True)
        outpath = os.path.join(outdir, 'fig5_aerosol_paper.png')
        fig.savefig(outpath, dpi=200, bbox_inches='tight')
        print("Saved: %s" % outpath)
    plt.close()

    # === Self-computed version (total aerosol only) ===
    fig2, axes2 = plt.subplots(1, ncols, figsize=(5 * ncols, 3.5),
                               subplot_kw={'projection': ccrs.PlateCarree()},
                               squeeze=False)
    has_self = False
    for col, cn in enumerate(case_names):
        cfg = load_case(cn)
        kr = cfg.get('plot', {}).get('key_region')
        kr_list = (kr['lon'] + kr['lat']) if kr else None
        case_label = cfg.get('case_name', cn.upper())

        lats, lons, sdata = load_self_aerosol(cn)
        if sdata is not None and 'total' in sdata:
            has_self = True
            ax = axes2[0, col]
            field = np.clip(np.nan_to_num(sdata['total'], nan=0), -5, 5)
            plot_panel(ax, lons, lats, field, "Total Aerosol %s [self]" % case_label,
                      norm, CMAP, kr_list)

    if has_self:
        fig2.subplots_adjust(bottom=0.15)
        cbar_ax2 = fig2.add_axes([0.15, 0.05, 0.7, 0.03])
        fig2.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cbar_ax2,
                     orientation='horizontal', ticks=LEVELS_AER,
                     label='Partial Temperature (K)')
        first_cfg = load_case(case_names[0])
        outpath2 = os.path.join(first_cfg['_figures_dir'], 'fig5_aerosol_self.png')
        fig2.savefig(outpath2, dpi=200, bbox_inches='tight')
        print("Saved: %s" % outpath2)
    plt.close()

    # === Validation: self total vs paper total ===
    if has_paper and has_self:
        print("\n=== Aerosol total validation ===")
        for cn in case_names:
            _, _, pdata = load_paper_aerosol(cn)
            _, _, sdata = load_self_aerosol(cn)
            if pdata and sdata and 'total' in pdata and 'total' in sdata:
                p = pdata['total']
                s = sdata['total']
                mask = np.isfinite(p) & np.isfinite(s) & (np.abs(p) < 900) & (np.abs(s) < 900)
                if mask.sum() > 0:
                    diff = s[mask] - p[mask]
                    corr = np.corrcoef(s[mask], p[mask])[0, 1]
                    cfg = load_case(cn)
                    print("  %s: RMSE=%.3f K, bias=%.3f K, corr=%.4f" %
                          (cfg.get('case_name', cn), np.sqrt(np.mean(diff**2)),
                           np.mean(diff), corr))

    # Domain-mean per species (paper data)
    if has_paper:
        print("\n=== Domain-mean aerosol dT (K) from paper_data ===")
        header = "%-15s" % "Species"
        for cn in case_names:
            cfg = load_case(cn)
            header += " %10s" % cfg.get('case_name', cn)
        print(header)
        for key, label in AER_ROWS:
            row_str = "%-15s" % label
            for cn in case_names:
                _, _, pdata = load_paper_aerosol(cn)
                if pdata and key in pdata:
                    row_str += " %10.3f" % np.nanmean(pdata[key])
                else:
                    row_str += " %10s" % "N/A"
            print(row_str)


if __name__ == '__main__':
    main()
