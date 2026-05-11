#!/usr/bin/env python3
"""Single-column drdt diagnostic.

Picks one grid cell from cesm2_4xco2_official input (default: tropical ocean),
constructs data_prep/*.dat manually, runs cfram_fu_1col binary directly, and
inspects the resulting drdt.dat matrix.

Compares drdt diagonal against theoretical Stefan-Boltzmann surface kernel:
  drdt(nv1, nv1) ≈ -4σTs³  (W/m²/K, with σ = 5.6704e-8)

This tests candidate 4 of the dT under-estimation hypothesis: is the
pyCFRAM Fu surface drdt biased compared to expected blackbody scaling?

Usage: python3 scripts/diag_drdt_singlecol.py [ilat] [ilon]
"""
import os
import sys
import shutil
import subprocess
import tempfile
import numpy as np
from netCDF4 import Dataset

ROOT = '/home/lzhenn/work/ust-jumper/pyCFRAM'
CASE_INPUT = os.path.join(ROOT, 'cases/cesm2_4xco2_official/input')
BINARY = os.path.join(ROOT, 'fortran/cfram_fu_1col')
PLEV_HPA = np.array([1, 5, 10, 20, 30, 50, 70, 100, 150, 200, 250, 300,
                     400, 500, 600, 700, 850, 925, 1000], dtype=np.float64)
NLEV = len(PLEV_HPA)
SIGMA = 5.6704e-8


def write_bin(path, arr):
    np.asarray(arr, dtype=np.float64).tofile(path)


def read_drdt(path):
    """Read Fortran sequential record file: int8 header (nv) + double matrix
    (nv1, nv1) F-order."""
    with open(path, 'rb') as f:
        # Header record: marker(4) + int8(8) + marker(4) = 16 bytes
        rec1 = np.frombuffer(f.read(16), dtype=np.uint8)
        nv = int(np.frombuffer(rec1[4:12], dtype=np.int64)[0])
        nv1 = nv + 1
        # Body record: marker(4) + nv1*nv1 doubles + marker(4)
        body = f.read()
        # Strip 4-byte markers at start and end
        data = np.frombuffer(body[4:-4], dtype=np.float64)
        if data.size != nv1 * nv1:
            raise ValueError(f'drdt body size {data.size} != {nv1*nv1}')
        # Fortran wrote ((drdt(k, kk), kk=1, nv1), k=1, nv1):
        # outermost loop is k, innermost is kk → row-major in (k, kk).
        # Matrix shape (nv1, nv1) where drdt[k, kk] = ∂R(k)/∂T(kk).
        drdt = data.reshape((nv1, nv1))
    return drdt, nv


def extract_column_to_data_prep(ilat, ilon, dp):
    """Extract one (ilat, ilon) column from CMIP6 input NCs and write the
    Fortran data_prep directory."""
    nc_bp = Dataset(os.path.join(CASE_INPUT, 'base_pres.nc'))
    nc_ap = Dataset(os.path.join(CASE_INPUT, 'perturbed_pres.nc'))
    nc_bs = Dataset(os.path.join(CASE_INPUT, 'base_surf.nc'))
    nc_as = Dataset(os.path.join(CASE_INPUT, 'perturbed_surf.nc'))

    def get3d(nc, v):
        return np.array(nc.variables[v][0, ::-1, :, :], dtype=np.float64)

    def get2d(nc, v, default=None):
        if v in nc.variables:
            return np.array(nc.variables[v][0, :, :], dtype=np.float64)
        if default is not None:
            return np.full(nc.variables[list(nc.variables)[-1]][0].shape, default)
        raise KeyError(v)

    t_b   = get3d(nc_bp, 'ta_lay')[:, ilat, ilon]
    q_b   = get3d(nc_bp, 'q')      [:, ilat, ilon]
    o3_b  = get3d(nc_bp, 'o3')     [:, ilat, ilon]
    cc_b  = get3d(nc_bp, 'camt')   [:, ilat, ilon]
    clw_b = get3d(nc_bp, 'cliq')   [:, ilat, ilon]
    ciw_b = get3d(nc_bp, 'cice')   [:, ilat, ilon]
    co2_b = get3d(nc_bp, 'co2').mean() * 1e6  # mol/mol → ppmv (well-mixed)

    t_w   = get3d(nc_ap, 'ta_lay')[:, ilat, ilon]
    q_w   = get3d(nc_ap, 'q')      [:, ilat, ilon]
    o3_w  = get3d(nc_ap, 'o3')     [:, ilat, ilon]
    cc_w  = get3d(nc_ap, 'camt')   [:, ilat, ilon]
    clw_w = get3d(nc_ap, 'cliq')   [:, ilat, ilon]
    ciw_w = get3d(nc_ap, 'cice')   [:, ilat, ilon]
    co2_w = get3d(nc_ap, 'co2').mean() * 1e6

    ts_b    = get2d(nc_bs, 'ts')[ilat, ilon]
    ps_b    = get2d(nc_bs, 'ps')[ilat, ilon]
    solar_b = get2d(nc_bs, 'solar')[ilat, ilon]
    alb_b   = get2d(nc_bs, 'albedo')[ilat, ilon]
    huss_b  = (get2d(nc_bs, 'huss')[ilat, ilon]
               if 'huss' in nc_bs.variables else -999.0)

    ts_w    = get2d(nc_as, 'ts')[ilat, ilon]
    ps_w    = get2d(nc_as, 'ps')[ilat, ilon]
    solar_w = get2d(nc_as, 'solar')[ilat, ilon]
    alb_w   = get2d(nc_as, 'albedo')[ilat, ilon]
    huss_w  = (get2d(nc_as, 'huss')[ilat, ilon]
               if 'huss' in nc_as.variables else -999.0)

    print(f'  Cell ({ilat},{ilon}): ts_b={ts_b:.2f}K  ps_b={ps_b/100:.1f}hPa')
    print(f'    co2_b={co2_b:.1f}  co2_w={co2_w:.1f}  alb_b={alb_b:.3f}')
    print(f'    huss_b={huss_b:.4e}  q_b[bot]={q_b[NLEV-1]:.4e}')

    nc_bp.close(); nc_ap.close(); nc_bs.close(); nc_as.close()

    write_bin(os.path.join(dp, 'plev.dat'), PLEV_HPA)
    write_bin(os.path.join(dp, 't_base.dat'), t_b)
    write_bin(os.path.join(dp, 't_warm.dat'), t_w)
    write_bin(os.path.join(dp, 'hus_base.dat'), q_b)
    write_bin(os.path.join(dp, 'hus_warm.dat'), q_w)
    write_bin(os.path.join(dp, 'O3_base.dat'), o3_b)
    write_bin(os.path.join(dp, 'O3_warm.dat'), o3_w)
    write_bin(os.path.join(dp, 'cc_base.dat'), cc_b)
    write_bin(os.path.join(dp, 'cc_warm.dat'), cc_w)
    write_bin(os.path.join(dp, 'clwc_base.dat'), clw_b)
    write_bin(os.path.join(dp, 'clwc_warm.dat'), clw_w)
    write_bin(os.path.join(dp, 'ciwc_base.dat'), ciw_b)
    write_bin(os.path.join(dp, 'ciwc_warm.dat'), ciw_w)
    write_bin(os.path.join(dp, 'skt_base.dat'), [ts_b])
    write_bin(os.path.join(dp, 'skt_warm.dat'), [ts_w])
    write_bin(os.path.join(dp, 'sp_base.dat'), [ps_b])
    write_bin(os.path.join(dp, 'sp_warm.dat'), [ps_w])
    write_bin(os.path.join(dp, 'solarin_base.dat'), [solar_b])
    write_bin(os.path.join(dp, 'solarin_warm.dat'), [solar_w])
    # ssrd/ssru hardcoded same as run_parallel_python: solar_in*alb / 300.0
    write_bin(os.path.join(dp, 'ssrd_base.dat'), [300.0])
    write_bin(os.path.join(dp, 'ssrd_warm.dat'), [300.0])
    write_bin(os.path.join(dp, 'ssru_base.dat'), [300.0 * alb_b])
    write_bin(os.path.join(dp, 'ssru_warm.dat'), [300.0 * alb_w])
    write_bin(os.path.join(dp, 'co2_b.dat'), [co2_b])
    write_bin(os.path.join(dp, 'co2_w.dat'), [co2_w])
    write_bin(os.path.join(dp, 'huss_base.dat'), [huss_b])
    write_bin(os.path.join(dp, 'huss_warm.dat'), [huss_w])

    return ts_b, ps_b


def main():
    ilat = int(sys.argv[1]) if len(sys.argv) > 1 else 96   # equator-ish
    ilon = int(sys.argv[2]) if len(sys.argv) > 2 else 0    # 0 deg E
    print(f'\n=== drdt single-column diagnostic at ilat={ilat}, ilon={ilon} ===')

    tmp = tempfile.mkdtemp(prefix='drdt_diag_')
    print(f'tmp dir: {tmp}')
    dp = os.path.join(tmp, 'data_prep')
    do = os.path.join(tmp, 'data_output')
    os.makedirs(dp); os.makedirs(do)

    print('\n1. Extracting column...')
    ts, ps = extract_column_to_data_prep(ilat, ilon, dp)

    print('\n2. Running cfram_fu_1col...')
    cwd = os.getcwd()
    os.chdir(tmp)
    try:
        result = subprocess.run([BINARY], capture_output=True, text=True, timeout=60)
        if result.returncode != 0:
            print(f'BINARY FAILED rc={result.returncode}')
            print('STDOUT:', result.stdout[-500:])
            print('STDERR:', result.stderr[-500:])
            sys.exit(1)
    finally:
        os.chdir(cwd)
    print('  binary OK')

    print('\n3. Reading drdt.dat...')
    drdt, nv = read_drdt(os.path.join(do, 'drdt.dat'))
    nv1 = nv + 1
    print(f'  nv={nv}, nv1={nv1}, drdt shape={drdt.shape}')

    print('\n4. Diagonal of drdt:')
    print('  [k]  pp(k)hPa    drdt(k,k) W/m²/K   role')
    diag = np.diag(drdt)
    for k in range(nv1):
        if k < nv:
            label = f'p={PLEV_HPA[k]:.0f}hPa'
        else:
            label = f'SURFACE ps={ps/100:.1f}hPa'
        print(f'  [{k:2d}] {label:18s} {diag[k]:+10.4f}')

    sb_expected = -4.0 * SIGMA * ts**3
    print(f'\n5. Surface drdt verification:')
    print(f'  drdt(nv1, nv1)        = {diag[nv]:+.4f} W/m²/K')
    print(f'  -4σTs³ (Ts={ts:.1f}K) = {sb_expected:+.4f} W/m²/K  ← skin-only')
    print(f'  ratio drdt/(-4σTs³)   = {diag[nv]/sb_expected:.3f}')
    print('  Note: drdt is NET (skin emission ULW − atm-bottom DLW response).')
    print('        With Fix #2 (pt(nv1)=ts_use_p), atm DLW partially offsets,')
    print('        so |drdt| < |-4σTs³|. Typical ratio ~0.5-0.8 for moist BL.')

    print(f'\n6. Atm-row diagonal sample (every 4):')
    for k in range(0, nv, 4):
        print(f'  drdt(p={PLEV_HPA[k]:.0f}hPa, T+1K) = {diag[k]:+.4f} W/m²/K')

    # Save matrix for further inspection
    out_npy = os.path.join(cwd, f'/tmp/drdt_diag_{ilat}_{ilon}.npy')
    np.save(out_npy, drdt)
    print(f'\nSaved drdt matrix to: {out_npy}')

    shutil.rmtree(tmp)


if __name__ == '__main__':
    main()
