#!/usr/bin/env python3
"""Plot zonal-mean dT_dry only, troposphere (1000-100 hPa, linear y-axis).

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

    a = np.array(nc.variables['dT_dry'][:])
    a = np.where(np.abs(a) > 900, np.nan, a)
    zm = np.nanmean(a, axis=2)
    nc.close()

    fig, ax = plt.subplots(figsize=(8, 6))
    norm = mcolors.BoundaryNorm(LEVELS, ncolors=CMAP.N, extend='both')
    cf = ax.contourf(lats, lev, zm, levels=LEVELS, cmap=CMAP, norm=norm, extend='both')
    cb = plt.colorbar(cf, ax=ax, label='K', ticks=LEVELS)
    ax.set_xlabel('Latitude')
    ax.set_ylabel('Pressure (hPa)')
    # Linear y-axis, troposphere only (1000-100 hPa, surface at bottom)
    ax.set_ylim(1000, 100)
    ax.set_title('dT_dry (Lu/Cai full-state ΔQ)', fontsize=12)
    plt.tight_layout()

    out = os.path.join(cfg['_figures_dir'], 'diag_zonalmean_dT_dry_trop.png')
    os.makedirs(cfg['_figures_dir'], exist_ok=True)
    plt.savefig(out, dpi=130, bbox_inches='tight')
    print('Saved:', out)


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'cesm2_4xco2')
