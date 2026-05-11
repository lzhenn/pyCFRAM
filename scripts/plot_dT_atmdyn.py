#!/usr/bin/env python3
"""Plot zonal-mean dT_atmdyn (OLD CFRAM dt_atm_dyn convention) — full column.

After the OLD CFRAM alignment (which uses frc_full with T_warm), dT_atmdyn
matches OLD CFRAM's dt_atm_dyn pattern. Y-axis log scale, 1000 hPa to model top.

Color levels match OLD CFRAM reference: [-15, -7, -2, -0.5, 0, 0.5, 2, 7, 15].
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


def main(case_name):
    cfg = load_case(case_name)
    nc = Dataset(os.path.join(cfg['_output_dir'], 'cfram_result.nc'))
    lats = np.array(nc.variables['lat'][:])
    lev = np.array(nc.variables['lev'][:])

    a = np.array(nc.variables['dT_atmdyn'][:])
    a = np.where(np.abs(a) > 900, np.nan, a)
    zm = np.nanmean(a, axis=2)
    nc.close()

    fig, ax = plt.subplots(figsize=(8, 6))
    norm = mcolors.BoundaryNorm(LEVELS, ncolors=CMAP.N, extend='both')
    cf = ax.contourf(lats, lev, zm, levels=LEVELS, cmap=CMAP, norm=norm, extend='both')
    cb = plt.colorbar(cf, ax=ax, label='K', ticks=LEVELS)
    ax.set_xlabel('Latitude')
    ax.set_ylabel('Pressure (hPa)')
    # Full column: surface (1000 hPa) at bottom, model top (~1 hPa) at top.
    ax.set_yscale('log')
    ax.set_ylim(1000, lev.min())
    ax.set_yticks([1000, 500, 200, 100, 50, 20, 10, 5, 2, 1])
    ax.get_yaxis().set_major_formatter(plt.matplotlib.ticker.ScalarFormatter())
    ax.set_title('dT_atmdyn (OLD CFRAM dt_atm_dyn — full-state, sfc=0)', fontsize=12)
    plt.tight_layout()

    out = os.path.join(cfg['_figures_dir'], 'diag_zonalmean_dT_atmdyn.png')
    os.makedirs(cfg['_figures_dir'], exist_ok=True)
    plt.savefig(out, dpi=130, bbox_inches='tight')
    print('Saved:', out)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'cesm2_4xco2')
