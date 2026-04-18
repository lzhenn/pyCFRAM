#!/usr/bin/env python3
"""Verify cloud LW/SW split: bulk == LW + SW to float64 precision.

Also checks regression of all pre-existing variables against a baseline.

Usage:
    python3 scripts/verify_cloud_split.py <baseline.nc> <new.nc>

The baseline should be from before the cloud-split Fortran was introduced
(i.e. a prior Phase 4 run). The new result must contain frc_cloud_lw /
frc_cloud_sw and dT_cloud_lw / dT_cloud_sw.
"""
import os, sys
import numpy as np
from netCDF4 import Dataset

COMMON_FRC = ['co2', 'q', 'ts', 'o3', 'solar', 'albedo', 'cloud', 'aerosol', 'warm',
              'bc', 'ocphi', 'ocpho', 'sulf', 'ss', 'dust']
COMMON_DT = COMMON_FRC + ['lhflx', 'shflx', 'atmdyn', 'sfcdyn', 'ocndyn']


def main():
    if len(sys.argv) != 3:
        print('Usage: verify_cloud_split.py <baseline.nc> <new.nc>')
        sys.exit(2)
    base, new = sys.argv[1], sys.argv[2]
    b = Dataset(base); n = Dataset(new)
    fail = False

    print('=== Regression vs baseline (pre-cloud-split) ===')
    for t in COMMON_FRC:
        vn = 'frc_' + t
        if vn not in b.variables or vn not in n.variables:
            continue
        fb = np.array(b.variables[vn][:]); fn = np.array(n.variables[vn][:])
        mask = (np.abs(fb) < 900) & (np.abs(fn) < 900)
        if not mask.any():
            continue
        diff = float(np.abs(fb[mask] - fn[mask]).max())
        scale = float(np.abs(fb[mask]).max())
        rel = diff / (scale + 1e-30)
        st = 'OK' if rel < 1e-12 else 'DIFF'
        if st != 'OK':
            fail = True
        print('  %-14s: rel=%.2e [%s]' % (vn, rel, st))

    for t in COMMON_DT:
        vn = 'dT_' + t
        if vn not in b.variables or vn not in n.variables:
            continue
        fb = np.array(b.variables[vn][:]); fn = np.array(n.variables[vn][:])
        mask = np.isfinite(fb) & np.isfinite(fn) & (np.abs(fb) < 900) & (np.abs(fn) < 900)
        if not mask.any():
            continue
        diff = float(np.abs(fb[mask] - fn[mask]).max())
        scale = float(np.abs(fb[mask]).max())
        rel = diff / (scale + 1e-30)
        st = 'OK' if rel < 1e-10 else 'DIFF'
        if st != 'OK':
            fail = True
        print('  %-14s: rel=%.2e [%s]' % (vn, rel, st))

    print('\n=== New cloud LW/SW vars ===')
    for vn in ('frc_cloud_lw', 'frc_cloud_sw', 'dT_cloud_lw', 'dT_cloud_sw'):
        if vn not in n.variables:
            print('  %s MISSING' % vn); fail = True; continue
        a = np.array(n.variables[vn][:])
        mask = np.isfinite(a) & (np.abs(a) < 900)
        if not mask.any():
            print('  %s all invalid' % vn); fail = True; continue
        print('  %-14s: mean=%+.4e, max|.|=%.4e' % (
            vn, float(a[mask].mean()), float(np.abs(a[mask]).max())))

    print('\n=== Additivity: bulk == LW + SW (should be EXACT) ===')
    fc = np.array(n.variables['frc_cloud'][:])
    fl = np.array(n.variables['frc_cloud_lw'][:])
    fs = np.array(n.variables['frc_cloud_sw'][:])
    mask = np.abs(fc) < 900
    diff = float(np.abs(fc[mask] - fl[mask] - fs[mask]).max())
    scale = float(np.abs(fc[mask]).max())
    rel = diff / (scale + 1e-30)
    st = 'OK' if rel < 1e-12 else 'FAIL'
    if st != 'OK':
        fail = True
    print('  frc: max|cloud - (lw+sw)|=%.4e, rel=%.2e [%s]' % (diff, rel, st))

    dc = np.array(n.variables['dT_cloud'][:])
    dl = np.array(n.variables['dT_cloud_lw'][:])
    ds = np.array(n.variables['dT_cloud_sw'][:])
    mask = np.isfinite(dc) & (np.abs(dc) < 900)
    diff = float(np.abs(dc[mask] - dl[mask] - ds[mask]).max())
    scale = float(np.abs(dc[mask]).max())
    rel = diff / (scale + 1e-30)
    st = 'OK' if rel < 1e-10 else 'FAIL'
    if st != 'OK':
        fail = True
    print('  dT:  max|cloud - (lw+sw)|=%.4e, rel=%.2e [%s]' % (diff, rel, st))

    print('\n=== Surface dT domain mean ===')
    sfc = -1
    print('  dT_cloud (bulk) : %+.4f K' % float(np.nanmean(dc[sfc])))
    print('  dT_cloud_lw     : %+.4f K' % float(np.nanmean(dl[sfc])))
    print('  dT_cloud_sw     : %+.4f K' % float(np.nanmean(ds[sfc])))
    print('  sum lw+sw       : %+.4f K' % float(np.nanmean(dl[sfc] + ds[sfc])))

    b.close(); n.close()
    print('\n=== Verification: %s ===' % ('PASS' if not fail else 'FAIL'))
    sys.exit(0 if not fail else 1)


if __name__ == '__main__':
    main()
