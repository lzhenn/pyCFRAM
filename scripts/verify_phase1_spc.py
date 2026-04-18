#!/usr/bin/env python3
"""Phase 1 verification: per-species AOD sum equals bulk (physics invariant).

Reads fortran/data_prep/aerosol_*_{bulk,spc}.dat and verifies:
  sum over species of AOD equals bulk AOD
  AOD-weighted ssa and g reconstructions match bulk
"""
import os, sys
import numpy as np

DP = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'fortran', 'data_prep')
NSPECIES = 6
NBND_LW = 16
NBND_SW = 14


def load_bulk(fname, nbnd, nlev, nlat, nlon):
    # ndarray.tofile() always writes C-order regardless of array flags, so
    # np.asfortranarray(...).tofile(...) in the writer is effectively a no-op.
    # Read back in C-order (numpy default).
    arr = np.fromfile(os.path.join(DP, fname), dtype=np.float64)
    return arr.reshape((nlev, nlat, nlon, nbnd))


def load_spc(fname, nbnd, nlev, nlat, nlon):
    arr = np.fromfile(os.path.join(DP, fname), dtype=np.float64)
    return arr.reshape((nlev, nlat, nlon, nbnd, NSPECIES))


def infer_grid():
    # infer (nlev, nlat, nlon) from bulk aod_lw file size
    size = os.path.getsize(os.path.join(DP, 'aerosol_aod_lw_base.dat'))
    n = size // 8 // NBND_LW  # nlev * nlat * nlon
    # Try common grids
    for nlev in (37,):
        for nlat in (81, 97):
            if n % (nlev * nlat) == 0:
                nlon = n // (nlev * nlat)
                return nlev, nlat, nlon
    raise RuntimeError('cannot infer grid, n=%d' % n)


def main():
    nlev, nlat, nlon = infer_grid()
    print('Grid inferred: nlev=%d, nlat=%d, nlon=%d' % (nlev, nlat, nlon))

    ok = True
    for tag, nbnd in [('aod_lw', NBND_LW), ('aod_sw', NBND_SW)]:
        for state in ('base', 'warm'):
            bulk = load_bulk('aerosol_%s_%s.dat' % (tag, state), nbnd, nlev, nlat, nlon)
            spc = load_spc('aerosol_%s_%s_spc.dat' % (tag, state), nbnd, nlev, nlat, nlon)
            spc_sum = spc.sum(axis=-1)
            diff = float(np.abs(bulk - spc_sum).max())
            scale = float(np.abs(bulk).max())
            rel = diff / (scale + 1e-30)
            status = 'OK' if rel < 1e-10 else 'FAIL'
            if status == 'FAIL':
                ok = False
            print('  %s %s: bulk_max=%.4e, abs_diff=%.4e, rel=%.2e [%s]' % (
                tag, state, scale, diff, rel, status))

    # SSA reconstruction check (AOD-weighted)
    for state in ('base', 'warm'):
        aod_sw = load_bulk('aerosol_aod_sw_%s.dat' % state, NBND_SW, nlev, nlat, nlon)
        aod_sw_spc = load_spc('aerosol_aod_sw_%s_spc.dat' % state, NBND_SW, nlev, nlat, nlon)
        ssa_bulk = load_bulk('aerosol_ssa_sw_%s.dat' % state, NBND_SW, nlev, nlat, nlon)
        ssa_spc = load_spc('aerosol_ssa_sw_%s_spc.dat' % state, NBND_SW, nlev, nlat, nlon)
        with np.errstate(divide='ignore', invalid='ignore'):
            recon = np.where(aod_sw > 0, (aod_sw_spc * ssa_spc).sum(axis=-1) / aod_sw, 0.0)
        diff = float(np.abs(ssa_bulk - recon).max())
        status = 'OK' if diff < 1e-10 else 'FAIL'
        if status == 'FAIL':
            ok = False
        print('  ssa %s reconstruction: max_diff=%.4e [%s]' % (state, diff, status))

        g_bulk = load_bulk('aerosol_g_sw_%s.dat' % state, NBND_SW, nlev, nlat, nlon)
        g_spc = load_spc('aerosol_g_sw_%s_spc.dat' % state, NBND_SW, nlev, nlat, nlon)
        aod_ssa = aod_sw * ssa_bulk
        with np.errstate(divide='ignore', invalid='ignore'):
            g_recon = np.where(aod_ssa > 0,
                               (aod_sw_spc * ssa_spc * g_spc).sum(axis=-1) / aod_ssa, 0.0)
        diff = float(np.abs(g_bulk - g_recon).max())
        status = 'OK' if diff < 1e-10 else 'FAIL'
        if status == 'FAIL':
            ok = False
        print('  g   %s reconstruction: max_diff=%.4e [%s]' % (state, diff, status))

    print('\n=== Phase 1 verification: %s ===' % ('PASS' if ok else 'FAIL'))
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
