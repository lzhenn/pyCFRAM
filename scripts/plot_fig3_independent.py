#!/usr/bin/env python3
"""Plot Fig.3 comparing paper vs INDEPENDENT (ERA5+MERRA-2) CFRAM results.

Left column  : paper merra2_<case>_partial_t.nc (surface = lev[0])
Right column : cases/<case>/output/cfram_result.nc (surface = lev[-1])

Rows: Water Vapor, Cloud, Aerosol, Other Radiative (CO2+Albedo+Solar),
      Surface Process (LHFLX+SHFLX), Total.

Usage: python3 scripts/plot_fig3_independent.py eh13
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

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case, defaults

_defaults = defaults()
LEVELS_DT = _defaults['plotting']['levels_dt']
CMAP_NAME = _defaults['plotting']['cmap']
CMAP = plt.cm.get_cmap(CMAP_NAME)

PAPER_CASE_DIR = {
    'eh13': 'paper_data/cfram_out/case_eh13_c20250102/merra2_eh13_partial_t.nc',
    'eh22': 'paper_data/cfram_out/case_eh22_c20250118/merra2_eh22_partial_t.nc',
}

ROWS = [
    ('wv',    'Water Vapor'),
    ('cld',   'Cloud'),
    ('aer',   'Aerosol'),
    ('other', 'Other Radiative (CO2+Albedo+Solar)'),
    ('sfc',   'Surface Process (LHFLX+SHFLX)'),
    ('total', 'Total'),
]


def _clean(a):
    a = np.asarray(a, dtype=np.float64)
    return np.where(np.abs(a) > 900, np.nan, a)


def load_paper(case, root):
    nc = Dataset(os.path.join(root, PAPER_CASE_DIR[case]))
    lat = np.array(nc.variables['lat'][:])
    lon = np.array(nc.variables['lon'][:])
    # surface = lev[0]; variables shaped (time, lev, lat, lon)
    def g(v):
        return _clean(nc.variables[v][0, 0, :, :])

    d = {}
    d['wv']  = g('wv')
    d['cld'] = g('cld')
    aer = np.zeros_like(d['wv'])
    for v in ('bc','oc','sulf','seas','dust'):
        aer = aer + np.nan_to_num(g(v), nan=0.0)
    d['aer'] = aer
    other = np.zeros_like(d['wv'])
    for v in ('co2','albedo','solar'):
        other = other + np.nan_to_num(g(v), nan=0.0)
    d['other'] = other
    sfc = np.zeros_like(d['wv'])
    for v in ('lhflx','shflx'):
        sfc = sfc + np.nan_to_num(g(v), nan=0.0)
    d['sfc'] = sfc
    # Total = sum of all 12 named partial_t terms
    total = np.zeros_like(d['wv'])
    for v in ('cld','wv','co2','albedo','solar','bc','oc','sulf','seas','dust','lhflx','shflx'):
        total = total + np.nan_to_num(g(v), nan=0.0)
    # add ozone if present
    if 'ozone' in nc.variables:
        total = total + np.nan_to_num(g('ozone'), nan=0.0)
    d['total'] = total
    nc.close()
    return lat, lon, d


def load_independent(case, root):
    f = os.path.join(root, 'cases', case, 'output', 'cfram_result.nc')
    nc = Dataset(f)
    lat = np.array(nc.variables['lat'][:])
    lon = np.array(nc.variables['lon'][:])
    def g(v):
        return _clean(nc.variables[v][-1, :, :])  # surface = lev[-1]

    d = {}
    d['wv']  = g('dT_q')
    d['cld'] = g('dT_cloud')
    d['aer'] = g('dT_aerosol')
    other = np.zeros_like(d['wv'])
    for v in ('dT_co2','dT_albedo','dT_solar'):
        other = other + np.nan_to_num(g(v), nan=0.0)
    d['other'] = other
    sfc = np.zeros_like(d['wv'])
    for v in ('dT_lhflx','dT_shflx'):
        sfc = sfc + np.nan_to_num(g(v), nan=0.0)
    d['sfc'] = sfc
    # Total radiative+surface-process sum (to be comparable to paper total)
    total = np.zeros_like(d['wv'])
    for v in ('dT_cloud','dT_q','dT_co2','dT_albedo','dT_solar',
              'dT_aerosol','dT_lhflx','dT_shflx'):
        total = total + np.nan_to_num(g(v), nan=0.0)
    d['total'] = total
    nc.close()
    return lat, lon, d


def crop_to_paper(lat_i, lon_i, d_i, lat_p, lon_p):
    lat_mask = (lat_i >= lat_p.min() - 1e-3) & (lat_i <= lat_p.max() + 1e-3)
    lon_mask = (lon_i >= lon_p.min() - 1e-3) & (lon_i <= lon_p.max() + 1e-3)
    out = {k: v[np.ix_(lat_mask, lon_mask)] for k, v in d_i.items()}
    return lat_i[lat_mask], lon_i[lon_mask], out


def regrid_to_paper(lat_i, lon_i, field_i, lat_p, lon_p):
    """Bilinear-ish nearest regrid (indep 0.25, paper 0.25) so arrays align for corr."""
    from scipy.interpolate import RegularGridInterpolator
    interp = RegularGridInterpolator((lat_i, lon_i), np.nan_to_num(field_i, nan=0.0),
                                     bounds_error=False, fill_value=0.0)
    LatP, LonP = np.meshgrid(lat_p, lon_p, indexing='ij')
    return interp((LatP, LonP))


def spatial_corr(a, b):
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 5:
        return np.nan
    a = a[m]; b = b[m]
    if a.std() == 0 or b.std() == 0:
        return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def plot_panel(ax, lons, lats, field, title, norm, cmap):
    ax.set_extent([95, 125, 20, 40], crs=ccrs.PlateCarree())
    field = np.clip(np.nan_to_num(field, nan=0.0), -20, 20)
    ax.contourf(lons, lats, field, levels=LEVELS_DT, cmap=cmap, norm=norm,
                transform=ccrs.PlateCarree(), extend='both')
    ax.coastlines(resolution='50m', linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3, edgecolor='gray')
    gl = ax.gridlines(draw_labels=True, linewidth=0.3, color='gray', alpha=0.5, linestyle='--')
    gl.top_labels = False; gl.right_labels = False
    gl.xlocator = plt.FixedLocator([95,100,105,110,115,120,125])
    gl.ylocator = plt.FixedLocator([20,25,30,35,40])
    gl.xformatter = LongitudeFormatter(); gl.yformatter = LatitudeFormatter()
    gl.xlabel_style = {'size': 7}; gl.ylabel_style = {'size': 7}
    ax.set_title(title, fontsize=8.5, fontweight='bold')


def main():
    if len(sys.argv) < 2:
        print("Usage: plot_fig3_independent.py <case>  (eh13|eh22)")
        sys.exit(1)
    case = sys.argv[1].lower()
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    lat_p, lon_p, dp = load_paper(case, root)
    lat_i, lon_i, di = load_independent(case, root)
    # Crop to paper bounds
    lat_ic, lon_ic, dic = crop_to_paper(lat_i, lon_i, di, lat_p, lon_p)

    # Regrid each independent field to paper grid for correlation
    dic_rg = {k: regrid_to_paper(lat_ic, lon_ic, v, lat_p, lon_p) for k, v in dic.items()}

    # Consistent colorscale per row => use paper LEVELS_DT for both
    norm = mcolors.BoundaryNorm(LEVELS_DT, CMAP.N, clip=True)

    nrows = len(ROWS); ncols = 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 2.8 * nrows),
                             subplot_kw={'projection': ccrs.PlateCarree()},
                             squeeze=False)

    case_label = case.upper()
    stats = []
    for r, (key, lbl) in enumerate(ROWS):
        fp = dp[key]; fi = dic[key]; fi_rg = dic_rg[key]
        corr = spatial_corr(fp, fi_rg)
        mean_p = float(np.nanmean(fp)); mean_i = float(np.nanmean(fi))
        stats.append((lbl, mean_p, mean_i, corr))

        plot_panel(axes[r, 0], lon_p, lat_p, fp,
                   "(%s1) %s — Paper (mean=%.2fK)" % ('abcdef'[r], lbl, mean_p),
                   norm, CMAP)
        plot_panel(axes[r, 1], lon_ic, lat_ic, fi,
                   "(%s2) %s — Independent (mean=%.2fK, r=%.3f)" % ('abcdef'[r], lbl, mean_i, corr),
                   norm, CMAP)

    fig.suptitle("Fig.3 Surface DeltaT decomposition - %s  (Paper vs Independent ERA5+MERRA-2)" % case_label,
                 fontsize=12, fontweight='bold', y=0.995)
    fig.subplots_adjust(bottom=0.05, top=0.965, left=0.05, right=0.95, hspace=0.32, wspace=0.12)
    cbar_ax = fig.add_axes([0.15, 0.015, 0.7, 0.010])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=CMAP), cax=cbar_ax,
                      orientation='horizontal', ticks=LEVELS_DT)
    cb.set_label('Partial Temperature (K)', fontsize=10)

    outdir = os.path.join(root, 'cases', case, 'figures')
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, 'fig3_independent.png')
    fig.savefig(outpath, dpi=180, bbox_inches='tight')
    plt.close()
    print("Saved: %s" % outpath)
    print("\nSummary (surface dT, K):")
    print("%-36s %10s %12s %8s" % ("term", "paper_mean", "indep_mean", "corr"))
    for lbl, mp, mi, c in stats:
        print("%-36s %10.3f %12.3f %8.3f" % (lbl, mp, mi, c))


if __name__ == '__main__':
    main()
