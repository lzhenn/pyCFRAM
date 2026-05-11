#!/usr/bin/env python3
"""Compare dT_dry (Lu/Cai full-state ΔQ) vs dT_atmdyn (legacy lumping).

Two plots side-by-side: zonal-mean profile of dT_atmdyn and dT_dry, with
matching colorbar to OLD CFRAM reference (-15..+15 K with discrete levels).
"""
import os
import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case


LEVELS = [-15, -7, -2, -0.5, 0, 0.5, 2, 7, 15]
CMAP = plt.cm.RdBu_r


def panel(ax, x, y, field, title, units='K'):
    norm = mcolors.BoundaryNorm(LEVELS, ncolors=CMAP.N, extend='both')
    cf = ax.contourf(x, y, field, levels=LEVELS, cmap=CMAP, norm=norm, extend='both')
    cb = plt.colorbar(cf, ax=ax, label=units, ticks=LEVELS)
    ax.set_xlabel('Latitude')
    ax.set_ylabel('Pressure (hPa)')
    ax.set_yscale('log')
    ax.set_ylim(1000, 1)
    ax.set_yticks([1000, 800, 600, 500, 400, 300, 200, 100, 50, 20, 10, 5, 2, 1])
    ax.get_yaxis().set_major_formatter(plt.matplotlib.ticker.ScalarFormatter())
    ax.set_title(title, fontsize=11)


def main(case_name):
    cfg = load_case(case_name)
    nc_file = os.path.join(cfg['_output_dir'], 'cfram_result.nc')
    nc = Dataset(nc_file)
    lats = np.array(nc.variables['lat'][:])
    lev = np.array(nc.variables['lev'][:])

    def zm(v):
        a = np.array(nc.variables[v][:])
        a = np.where(np.abs(a) > 900, np.nan, a)
        return np.nanmean(a, axis=2)

    dT_atmdyn = zm('dT_atmdyn')
    dT_dry    = zm('dT_dry')
    dT_obs    = zm('dT_observed') if 'dT_observed' in nc.variables else None

    # Σ partials with each convention
    rad_vars = ['dT_co2', 'dT_q', 'dT_o3', 'dT_cloud_lw', 'dT_cloud_sw',
                'dT_albedo', 'dT_solar', 'dT_aerosol', 'dT_ts']
    rad_sum = np.zeros_like(dT_atmdyn)
    for v in rad_vars:
        if v in nc.variables:
            rad_sum += np.nan_to_num(zm(v))
    sum_legacy = rad_sum + np.nan_to_num(dT_atmdyn)
    sum_lucai  = rad_sum + np.nan_to_num(dT_dry)

    fig, axes = plt.subplots(2, 3, figsize=(20, 11))
    panel(axes[0, 0], lats, lev, dT_atmdyn,
          'dT_atmdyn (legacy pyCFRAM lump)')
    panel(axes[0, 1], lats, lev, dT_dry,
          'dT_dry (Lu/Cai full-state ΔQ)')
    panel(axes[0, 2], lats, lev, dT_dry - dT_atmdyn,
          'dT_dry − dT_atmdyn (= drdt·ΔT correction)')

    if dT_obs is not None:
        panel(axes[1, 0], lats, lev, dT_obs, 'dT_observed (model T2-T1)')
        panel(axes[1, 1], lats, lev, sum_legacy - dT_obs,
              'Σ(rad+atmdyn) − obs (legacy residual)')
        panel(axes[1, 2], lats, lev, sum_lucai - dT_obs,
              'Σ(rad+dT_dry) − obs (Lu/Cai residual)')
    else:
        for j in range(3):
            axes[1, j].set_visible(False)

    fig.suptitle('CESM2 4×CO2: Lu/Cai full-state ΔQ vs legacy atmdyn lumping',
                 fontsize=13)
    plt.tight_layout()
    out = os.path.join(cfg['_figures_dir'], 'diag_dry_vs_atmdyn.png')
    os.makedirs(cfg['_figures_dir'], exist_ok=True)
    plt.savefig(out, dpi=130, bbox_inches='tight')
    print('Saved:', out)
    nc.close()


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'cesm2_4xco2')
