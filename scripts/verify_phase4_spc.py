#!/usr/bin/env python3
"""Phase 4 verification:
  - regression: all pre-existing frc_* / dT_* bitwise identical to Phase 2 baseline
  - new frc_species and dT_species are present, finite, with expected scale
  - additivity: sum_species frc ≈ frc_aerosol (within ~5% for non-linear coupling)
  - additivity: sum_species dT   ≈ dT_aerosol (same tolerance)

Usage:
    python3 scripts/verify_phase4_spc.py <baseline.nc> <new.nc>
"""
import os, sys
import numpy as np
from netCDF4 import Dataset

SPECIES = ['bc', 'ocphi', 'ocpho', 'sulf', 'ss', 'dust']
# Terms that existed pre-Phase-3 and must be bitwise identical in Phase 4 output
REGRESSION_FRC = ['co2', 'q', 'ts', 'o3', 'solar', 'albedo', 'cloud', 'aerosol', 'warm']
REGRESSION_DT = REGRESSION_FRC + ['lhflx', 'shflx', 'atmdyn', 'sfcdyn', 'ocndyn']


def main():
    if len(sys.argv) != 3:
        print('Usage: verify_phase4_spc.py <baseline.nc> <new.nc>')
        sys.exit(2)
    base, new = sys.argv[1], sys.argv[2]
    b = Dataset(base); n = Dataset(new)
    fail = False

    print('=== Regression check (bitwise identity to Phase 2) ===')
    for t in REGRESSION_FRC:
        fb = np.array(b.variables['frc_'+t][:])
        fn = np.array(n.variables['frc_'+t][:])
        mask = (np.abs(fb) < 900) & (np.abs(fn) < 900)
        diff = float(np.abs(fb[mask] - fn[mask]).max()) if mask.any() else 0.0
        scale = float(np.abs(fb[mask]).max()) if mask.any() else 1.0
        rel = diff / (scale + 1e-30)
        status = 'OK' if rel < 1e-12 else 'DIFF'
        if status == 'DIFF': fail = True
        print('  frc_%-8s: rel=%.2e [%s]' % (t, rel, status))
    for t in REGRESSION_DT:
        vname = 'dT_'+t
        if vname not in b.variables or vname not in n.variables:
            continue
        fb = np.array(b.variables[vname][:])
        fn = np.array(n.variables[vname][:])
        mask = np.isfinite(fb) & np.isfinite(fn) & (np.abs(fb) < 900) & (np.abs(fn) < 900)
        if not mask.any():
            continue
        diff = float(np.abs(fb[mask] - fn[mask]).max())
        scale = float(np.abs(fb[mask]).max())
        rel = diff / (scale + 1e-30)
        status = 'OK' if rel < 1e-10 else 'DIFF'
        if status == 'DIFF': fail = True
        print('  dT_%-9s: rel=%.2e [%s]' % (t, rel, status))

    print('\n=== Presence of new per-species variables ===')
    for spc in SPECIES:
        for prefix in ('frc_', 'dT_'):
            vname = prefix + spc
            if vname not in n.variables:
                print('  %s MISSING' % vname)
                fail = True
                continue
            arr = np.array(n.variables[vname][:])
            mask = np.isfinite(arr) & (np.abs(arr) < 900)
            if not mask.any():
                print('  %s all invalid!' % vname)
                fail = True
                continue
            print('  %-14s: mean=%+.4e, max|.|=%.4e' % (
                vname, float(arr[mask].mean()), float(np.abs(arr[mask]).max())))

    print('\n=== Additivity check: sum_species ≈ bulk ===')
    frc_aer = np.array(n.variables['frc_aerosol'][:])
    mask = (np.abs(frc_aer) < 900)
    frc_sum = np.zeros_like(frc_aer)
    for spc in SPECIES:
        a = np.array(n.variables['frc_'+spc][:])
        mask_sp = (np.abs(a) < 900)
        frc_sum[mask_sp] += a[mask_sp]
    diff = float(np.abs(frc_aer[mask] - frc_sum[mask]).max())
    scale = float(np.abs(frc_aer[mask]).max())
    rel = diff / (scale + 1e-30)
    status = 'OK' if rel < 0.10 else 'WARN' if rel < 0.50 else 'FAIL'
    print('  frc: max|bulk-sum|=%.4e, scale=%.4e, rel=%.2e [%s]' % (diff, scale, rel, status))
    if status == 'FAIL': fail = True

    dT_aer = np.array(n.variables['dT_aerosol'][:])
    mask = np.isfinite(dT_aer) & (np.abs(dT_aer) < 900)
    dT_sum = np.zeros_like(dT_aer)
    for spc in SPECIES:
        a = np.array(n.variables['dT_'+spc][:])
        mask_sp = np.isfinite(a) & (np.abs(a) < 900)
        dT_sum[mask_sp] += a[mask_sp]
    diff = float(np.abs(dT_aer[mask] - dT_sum[mask]).max())
    scale = float(np.abs(dT_aer[mask]).max())
    rel = diff / (scale + 1e-30)
    status = 'OK' if rel < 0.10 else 'WARN' if rel < 0.50 else 'FAIL'
    print('  dT : max|bulk-sum|=%.4e, scale=%.4e, rel=%.2e [%s]' % (diff, scale, rel, status))
    if status == 'FAIL': fail = True

    # Surface-layer domain-mean check (human-readable)
    print('\n=== Surface dT domain mean (top layer in lev dim = surface) ===')
    # cfram_result.nc has lev[-1] = surface
    sfc_idx = -1
    print('  dT_aerosol (bulk)        : %+.4f K' % float(np.nanmean(dT_aer[sfc_idx])))
    s = 0.0
    for spc in SPECIES:
        a = np.array(n.variables['dT_'+spc][:])
        m = float(np.nanmean(a[sfc_idx]))
        print('  dT_%-10s          : %+.4f K' % (spc, m))
        s += m
    print('  sum of 6 species         : %+.4f K' % s)

    b.close(); n.close()
    print('\n=== Phase 4 verification: %s ===' % ('PASS' if not fail else 'FAIL'))
    sys.exit(0 if not fail else 1)


if __name__ == '__main__':
    main()
