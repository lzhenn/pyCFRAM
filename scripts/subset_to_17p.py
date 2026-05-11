#!/usr/bin/env python3
"""Subset cesm2_4xco2_official input NCs from 19 plev to OLD CFRAM 17 plev
by dropping top 2 stratospheric layers (lev=1, 5 hPa).

Usage: python3 scripts/subset_to_17p.py
Reads:  cases/cesm2_4xco2_official/input/{base,perturbed}_{pres,surf}.nc
        + nonrad_forcing.nc
Writes: cases/cesm2_4xco2_official_17p_fu/input/<same names>.nc

Convention: NC stores lev in surface→TOA order (lev[0]=1000, lev[18]=1).
We keep k=0..16 (lev = 1000 down to 10 hPa) and drop k=17,18 (lev = 5, 1 hPa).
"""
import os
import shutil
from netCDF4 import Dataset
import numpy as np

ROOT = '/home/lzhenn/work/ust-jumper/pyCFRAM'
SRC = os.path.join(ROOT, 'cases/cesm2_4xco2_official/input')
DST = os.path.join(ROOT, 'cases/cesm2_4xco2_official_17p_fu/input')
KEEP_NLEV = 17  # lev[0..16] = surface 1000 down to 10 hPa


def subset_pres_nc(src_path, dst_path):
    src = Dataset(src_path)
    dst = Dataset(dst_path, 'w')
    # dimensions
    for name, dim in src.dimensions.items():
        if name == 'lev':
            dst.createDimension(name, KEEP_NLEV)
        else:
            dst.createDimension(name, len(dim))
    # variables
    for name, var in src.variables.items():
        v = dst.createVariable(name, var.datatype, var.dimensions)
        if 'lev' in var.dimensions:
            lev_axis = var.dimensions.index('lev')
            slc = [slice(None)] * var.ndim
            slc[lev_axis] = slice(0, KEEP_NLEV)
            v[:] = var[tuple(slc)]
        else:
            v[:] = var[:]
    src.close()
    dst.close()
    print(f'  Wrote {dst_path}')


def copy_surf_nc(src_path, dst_path):
    shutil.copy(src_path, dst_path)
    print(f'  Copied {dst_path}')


def main():
    os.makedirs(DST, exist_ok=True)
    print(f'SRC: {SRC}')
    print(f'DST: {DST}\n')

    print('Subsetting pres NCs (lev 19 → 17, drop top 2 strato):')
    subset_pres_nc(os.path.join(SRC, 'base_pres.nc'),
                   os.path.join(DST, 'base_pres.nc'))
    subset_pres_nc(os.path.join(SRC, 'perturbed_pres.nc'),
                   os.path.join(DST, 'perturbed_pres.nc'))
    subset_pres_nc(os.path.join(SRC, 'nonrad_forcing.nc'),
                   os.path.join(DST, 'nonrad_forcing.nc'))

    print('\nCopying surf NCs (no lev dim):')
    copy_surf_nc(os.path.join(SRC, 'base_surf.nc'),
                 os.path.join(DST, 'base_surf.nc'))
    copy_surf_nc(os.path.join(SRC, 'perturbed_surf.nc'),
                 os.path.join(DST, 'perturbed_surf.nc'))

    print('\nVerify lev ranges:')
    for f in ('base_pres.nc', 'perturbed_pres.nc', 'nonrad_forcing.nc'):
        nc = Dataset(os.path.join(DST, f))
        lev = np.array(nc.variables['lev'][:])
        print(f'  {f}: nlev={len(lev)}, lev[0]={lev[0]}, lev[-1]={lev[-1]}')
        nc.close()


if __name__ == '__main__':
    main()
