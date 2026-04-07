#!/usr/bin/env python3
"""Planck Matrix Validation (Task 1 of Phase 1 verification suite).

Runs a single-column CFRAM to obtain the Fortran-computed drdt_inv matrix,
then performs three checks:

  1. Physics check  : diagonal of drdt = inv(drdt_inv) must be negative
                      (each layer emits more LW when warmed → net flux out increases)
  2. Round-trip     : drdt_inv @ drdt ≈ I  (max off-diagonal error)
  3. Condition check: condition number of drdt_inv should be reasonable (< 1e6)

Usage:
    python3 scripts/validate_planck.py --case eh13 [--ilat 40 --ilon 60] [--plot]

The script re-uses the single-column infrastructure from run_parallel_python.py.
"""

import os, sys, argparse, struct, tempfile, shutil
import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case, defaults, get_plev, get_aerosol_map, \
    get_fortran_dir, get_lookup_dir, PROJECT_ROOT
from core.constants import NBND_LW, NBND_SW

FORTRAN_PLEV = get_plev()
NLEV = len(FORTRAN_PLEV)
_aer_cfg = get_aerosol_map()
AEROSOL_MAP = {k: (v['lw'], v['sw'], v['rd']) for k, v in _aer_cfg.items()}


# ── helpers (shared with run_parallel_python.py) ─────────────────────────────

def load_lookup_table(nc_path):
    nc = Dataset(nc_path, 'r')
    table = {v: np.array(nc.variables[v][:], dtype=np.float64)
             for v in ('bext','bsca','qext','qsca','g','rh','radius')
             if v in nc.variables}
    nc.close()
    return table


def compute_rh(q, t, p_hpa):
    w = q / (1.0 - q + 1e-30)
    e = w * p_hpa / (0.622 + w)
    es = 6.112 * np.exp(17.67 * (t - 273.15) / (t - 29.65 + 1e-10))
    return np.clip(e / (es + 1e-30), 0.0, 1.0)


def compute_aerosol_column(aer_data, rh, thick, dens, lut):
    LW_BANDS = slice(14, 30)
    SW_BANDS = slice(0, 14)
    aod_lw = np.zeros((NLEV, NBND_LW))
    aod_sw = np.zeros((NLEV, NBND_SW))
    aod_ssa = np.zeros((NLEV, NBND_SW))
    aod_g = np.zeros((NLEV, NBND_SW))
    for species, (lut_lw, lut_sw, rd) in AEROSOL_MAP.items():
        if species not in aer_data:
            continue
        mixing = aer_data[species]
        for lut_file, bands, nbnd, target in [
                (lut_lw, LW_BANDS, NBND_LW, 'lw'),
                (lut_sw, SW_BANDS, NBND_SW, 'sw')]:
            table = lut[lut_file]
            r_idx = (int(rd) - 1 if rd >= 1
                     else int(np.argmin(np.abs(table['radius'] - rd))))
            bext = table['bext'][r_idx, :, bands]
            rh_idx = np.array([np.argmin(np.abs(table['rh'] - r)) for r in rh])
            kext = np.array([bext[rh_idx[i], :] for i in range(NLEV)])
            aod_col = mixing[:, None] * dens[:, None] * 1e3 * thick[:, None] * kext * 1e-3
            aod_col = np.nan_to_num(aod_col)
            if target == 'lw':
                aod_lw += aod_col
            else:
                aod_sw += aod_col
                qext = table['qext'][r_idx, :, bands]
                qsca = table['qsca'][r_idx, :, bands]
                g_tab = table['g'][r_idx, :, bands]
                ssa_col = np.array([
                    np.where(qext[rh_idx[i], :] > 0,
                             qsca[rh_idx[i], :] / qext[rh_idx[i], :], 0.0)
                    for i in range(NLEV)])
                g_col = np.array([g_tab[rh_idx[i], :] for i in range(NLEV)])
                aod_ssa += aod_col * ssa_col
                aod_g += aod_col * ssa_col * g_col
    with np.errstate(divide='ignore', invalid='ignore'):
        eff_ssa = np.where(aod_sw > 0, aod_ssa / aod_sw, 0.0)
        eff_g = np.where(aod_ssa > 0, aod_g / aod_ssa, 0.0)
    return aod_lw, aod_sw, eff_ssa, eff_g


def write_bin(fp, data):
    np.asarray(data, dtype=np.float64).tofile(fp)


def read_fortran_seq(fp):
    with open(fp, 'rb') as f:
        raw = f.read()
    n = struct.unpack('i', raw[:4])[0]
    return np.frombuffer(raw[4:4+n], dtype=np.float64).copy()


def read_drdt_inv(path):
    with open(path, 'rb') as f:
        raw = f.read()
    pos = 0
    rec1_len = struct.unpack('i', raw[pos:pos+4])[0]; pos += 4
    if rec1_len == 8:
        nlayer_val = struct.unpack('q', raw[pos:pos+8])[0]; pos += 8
    else:
        nlayer_val = struct.unpack('i', raw[pos:pos+4])[0]; pos += 4
    pos += 4  # skip footer
    rec2_len = struct.unpack('i', raw[pos:pos+4])[0]; pos += 4
    n = nlayer_val + 1
    mat = (np.frombuffer(raw[pos:pos+rec2_len], dtype=np.float64)
             .copy().reshape(n, n, order='F'))
    return mat, nlayer_val


# ── single-column run ─────────────────────────────────────────────────────────

def run_single_column(cfg, ilat, ilon):
    """Extract column, run Fortran, return drdt_inv matrix."""
    LOOKUP_DIR = get_lookup_dir()
    lut = {}
    for species, (lut_lw, lut_sw, _) in AEROSOL_MAP.items():
        for lf in [lut_lw, lut_sw]:
            if lf not in lut:
                fp = os.path.join(LOOKUP_DIR, lf)
                if os.path.exists(fp):
                    lut[lf] = load_lookup_table(fp)

    nc_bp = Dataset(cfg['input']['base_pres'])
    nc_bs = Dataset(cfg['input']['base_surf'])
    nc_wp = Dataset(cfg['input']['perturbed_pres'])
    nc_ws = Dataset(cfg['input']['perturbed_surf'])

    def load3(nc, var):
        return np.array(nc.variables[var][0, ::-1, :, :], dtype=np.float64)
    def load2(nc, var):
        return np.array(nc.variables[var][0, :, :], dtype=np.float64)

    t_b = load3(nc_bp, 'ta_lay')[:, ilat, ilon]
    q_b = load3(nc_bp, 'q')[:, ilat, ilon]
    o3_b = load3(nc_bp, 'o3')[:, ilat, ilon]
    cc_b = load3(nc_bp, 'camt')[:, ilat, ilon]
    clwc_b = load3(nc_bp, 'cliq')[:, ilat, ilon]
    ciwc_b = load3(nc_bp, 'cice')[:, ilat, ilon]
    ts_b = load2(nc_bs, 'ts')[ilat, ilon]
    ps_b = load2(nc_bs, 'ps')[ilat, ilon]
    solar_b = load2(nc_bs, 'solar')[ilat, ilon]
    albedo_b = load2(nc_bs, 'albedo')[ilat, ilon]
    co2_base = float(np.array(nc_bp.variables['co2'][0, :, ilat, ilon]).mean() * 1e6
                     if 'co2' in nc_bp.variables
                     else defaults()['radiation'].get('co2_base', 395.0))

    t_w = load3(nc_wp, 'ta_lay')[:, ilat, ilon]
    q_w = load3(nc_wp, 'q')[:, ilat, ilon]
    o3_w = load3(nc_wp, 'o3')[:, ilat, ilon]
    cc_w = load3(nc_wp, 'camt')[:, ilat, ilon]
    clwc_w = load3(nc_wp, 'cliq')[:, ilat, ilon]
    ciwc_w = load3(nc_wp, 'cice')[:, ilat, ilon]
    ts_w = load2(nc_ws, 'ts')[ilat, ilon]
    ps_w = load2(nc_ws, 'ps')[ilat, ilon]
    solar_w = load2(nc_ws, 'solar')[ilat, ilon]
    albedo_w = load2(nc_ws, 'albedo')[ilat, ilon]
    co2_warm = float(np.array(nc_wp.variables['co2'][0, :, ilat, ilon]).mean() * 1e6
                     if 'co2' in nc_wp.variables
                     else defaults()['radiation'].get('co2_warm', 396.0))

    nc_bp.close(); nc_bs.close(); nc_wp.close(); nc_ws.close()

    # Aerosol optics
    p_hpa = FORTRAN_PLEV
    Rd = 287.05
    p_half = np.zeros(NLEV + 1)
    for i in range(NLEV - 1):
        p_half[i + 1] = 0.5 * (p_hpa[i] + p_hpa[i + 1])
    p_half[-1] = ps_b / 100.0
    delp = np.diff(p_half) * 100.0
    dens_b = (p_hpa * 100.0) / (Rd * t_b)
    thick_b = delp / (dens_b * 9.8)
    rh_b = compute_rh(q_b, t_b, p_hpa)

    # Try to load aerosol data
    nc_bp2 = Dataset(cfg['input']['base_pres'])
    nc_wp2 = Dataset(cfg['input']['perturbed_pres'])
    aer_vars = ['bc', 'oc', 'sulf', 'seas', 'dust']
    aer_b = {}
    aer_w = {}
    for s in aer_vars:
        if s in nc_bp2.variables:
            aer_b[s] = np.array(nc_bp2.variables[s][0, ::-1, ilat, ilon], dtype=np.float64)
            aer_w[s] = np.array(nc_wp2.variables[s][0, ::-1, ilat, ilon], dtype=np.float64)
    nc_bp2.close(); nc_wp2.close()

    if lut and aer_b:
        aod_lw_b, aod_sw_b, ssa_sw_b, g_sw_b = compute_aerosol_column(aer_b, rh_b, thick_b, dens_b, lut)
        p_half_w = p_half.copy(); p_half_w[-1] = ps_w / 100.0
        delp_w = np.diff(p_half_w) * 100.0
        dens_w = (p_hpa * 100.0) / (Rd * t_w)
        thick_w = delp_w / (dens_w * 9.8)
        rh_w = compute_rh(q_w, t_w, p_hpa)
        aod_lw_w, aod_sw_w, ssa_sw_w, g_sw_w = compute_aerosol_column(aer_w, rh_w, thick_w, dens_w, lut)
    else:
        aod_lw_b = np.zeros((NLEV, NBND_LW))
        aod_sw_b = np.zeros((NLEV, NBND_SW))
        ssa_sw_b = np.zeros((NLEV, NBND_SW))
        g_sw_b = np.zeros((NLEV, NBND_SW))
        aod_lw_w = aod_lw_b.copy()
        aod_sw_w = aod_sw_b.copy()
        ssa_sw_w = ssa_sw_b.copy()
        g_sw_w = g_sw_b.copy()

    tmpdir = tempfile.mkdtemp(prefix='val_planck_')
    dp = os.path.join(tmpdir, 'data_prep')
    do = os.path.join(tmpdir, 'data_output')
    os.makedirs(dp); os.makedirs(do)

    write_bin(os.path.join(dp, 't_base.dat'), t_b)
    write_bin(os.path.join(dp, 't_warm.dat'), t_w)
    write_bin(os.path.join(dp, 'hus_base.dat'), q_b)
    write_bin(os.path.join(dp, 'hus_warm.dat'), q_w)
    write_bin(os.path.join(dp, 'O3_base.dat'), o3_b)
    write_bin(os.path.join(dp, 'O3_warm.dat'), o3_w)
    write_bin(os.path.join(dp, 'cc_base.dat'), cc_b)
    write_bin(os.path.join(dp, 'cc_warm.dat'), cc_w)
    write_bin(os.path.join(dp, 'clwc_base.dat'), clwc_b)
    write_bin(os.path.join(dp, 'clwc_warm.dat'), clwc_w)
    write_bin(os.path.join(dp, 'ciwc_base.dat'), ciwc_b)
    write_bin(os.path.join(dp, 'ciwc_warm.dat'), ciwc_w)
    write_bin(os.path.join(dp, 'skt_base.dat'), np.array([ts_b]))
    write_bin(os.path.join(dp, 'skt_warm.dat'), np.array([ts_w]))
    write_bin(os.path.join(dp, 'sp_base.dat'), np.array([ps_b]))
    write_bin(os.path.join(dp, 'sp_warm.dat'), np.array([ps_w]))
    write_bin(os.path.join(dp, 'solarin_base.dat'), np.array([solar_b]))
    write_bin(os.path.join(dp, 'solarin_warm.dat'), np.array([solar_w]))
    ssrd = 300.0
    write_bin(os.path.join(dp, 'ssrd_base.dat'), np.array([ssrd]))
    write_bin(os.path.join(dp, 'ssrd_warm.dat'), np.array([ssrd]))
    write_bin(os.path.join(dp, 'ssru_base.dat'), np.array([ssrd * albedo_b]))
    write_bin(os.path.join(dp, 'ssru_warm.dat'), np.array([ssrd * albedo_w]))
    write_bin(os.path.join(dp, 'co2_b.dat'), np.array([co2_base]))
    write_bin(os.path.join(dp, 'co2_w.dat'), np.array([co2_warm]))
    write_bin(os.path.join(dp, 'aerosol_aod_lw_base.dat'), aod_lw_b)
    write_bin(os.path.join(dp, 'aerosol_aod_lw_warm.dat'), aod_lw_w)
    write_bin(os.path.join(dp, 'aerosol_aod_sw_base.dat'), aod_sw_b)
    write_bin(os.path.join(dp, 'aerosol_aod_sw_warm.dat'), aod_sw_w)
    write_bin(os.path.join(dp, 'aerosol_ssa_sw_base.dat'), ssa_sw_b)
    write_bin(os.path.join(dp, 'aerosol_ssa_sw_warm.dat'), ssa_sw_w)
    write_bin(os.path.join(dp, 'aerosol_g_sw_base.dat'), g_sw_b)
    write_bin(os.path.join(dp, 'aerosol_g_sw_warm.dat'), g_sw_w)

    exe_1col = os.path.join(get_fortran_dir(), 'cfram_rrtmg_1col')
    os.symlink(exe_1col, os.path.join(tmpdir, 'cfram_rrtmg'))
    ret = os.system('cd %s && ./cfram_rrtmg > /dev/null 2>&1' % tmpdir)

    drdt_inv = None
    nlayer = None
    if ret == 0:
        drdt_path = os.path.join(do, 'drdt_inv.dat')
        if os.path.exists(drdt_path):
            drdt_inv, nlayer = read_drdt_inv(drdt_path)

    shutil.rmtree(tmpdir)
    return drdt_inv, nlayer


# ── validation logic ──────────────────────────────────────────────────────────

def validate_matrix(drdt_inv, nlayer):
    """Three physical/mathematical checks. Returns list of (name, passed, detail)."""
    results = []
    n = nlayer + 1  # includes surface

    # 1. Physics: diagonal of drdt = inv(drdt_inv) must be negative
    drdt = np.linalg.inv(drdt_inv)
    diag = np.diag(drdt)
    n_neg = np.sum(diag < 0)
    all_neg = (n_neg == n)
    results.append((
        'Physics (diag of drdt < 0)',
        all_neg,
        f'{n_neg}/{n} diagonal elements negative  '
        f'[range: {diag.min():.4f} to {diag.max():.4f}]'
    ))

    # 2. Round-trip: drdt_inv @ drdt ≈ I
    product = drdt_inv @ drdt
    eye = np.eye(n)
    max_err = np.max(np.abs(product - eye))
    off_diag_max = np.max(np.abs(product - eye * np.diag(product)))
    passed = max_err < 1e-6
    results.append((
        'Round-trip (drdt_inv @ drdt ≈ I)',
        passed,
        f'max|product - I| = {max_err:.2e}  '
        f'off-diag max = {off_diag_max:.2e}'
    ))

    # 3. Condition number
    cond = np.linalg.cond(drdt_inv)
    passed = cond < 1e8
    results.append((
        'Condition number',
        passed,
        f'cond(drdt_inv) = {cond:.3e}  (threshold < 1e8)'
    ))

    return results, drdt


# ── optional plot ─────────────────────────────────────────────────────────────

def plot_matrix(drdt_inv, drdt, nlayer, out_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    n = nlayer + 1

    im0 = axes[0].imshow(drdt_inv, aspect='auto', cmap='RdBu_r',
                         vmin=-np.percentile(np.abs(drdt_inv), 99),
                         vmax=np.percentile(np.abs(drdt_inv), 99))
    axes[0].set_title(f'drdt_inv  ({n}×{n})')
    axes[0].set_xlabel('column j (layer j warmed)')
    axes[0].set_ylabel('row i (dT_i response)')
    plt.colorbar(im0, ax=axes[0], label='K / (W m⁻²)')

    diag = np.diag(drdt)
    axes[1].barh(np.arange(n), diag, color=np.where(diag < 0, 'steelblue', 'red'))
    axes[1].axvline(0, color='k', lw=0.8)
    axes[1].set_title('Diagonal of drdt = ∂R/∂T (should be all negative)')
    axes[1].set_xlabel('W m⁻² K⁻¹')
    axes[1].set_ylabel('level index (0=TOA, n-1=surface)')
    axes[1].invert_yaxis()

    plt.tight_layout()
    plt.savefig(out_path, dpi=120)
    print(f'Plot saved: {out_path}')


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--case', required=True, help='Case name (e.g. eh13)')
    parser.add_argument('--ilat', type=int, default=40, help='Latitude index (default=40)')
    parser.add_argument('--ilon', type=int, default=60, help='Longitude index (default=60)')
    parser.add_argument('--plot', action='store_true', help='Save diagnostic plot')
    args = parser.parse_args()

    cfg = load_case(args.case)
    OUTDIR = cfg['_output_dir']
    os.makedirs(OUTDIR, exist_ok=True)

    print(f'Case: {args.case}  column: (ilat={args.ilat}, ilon={args.ilon})')
    print('Running single-column Fortran...')
    drdt_inv, nlayer = run_single_column(cfg, args.ilat, args.ilon)

    if drdt_inv is None:
        print('FAIL: Fortran run did not produce drdt_inv.dat')
        sys.exit(1)

    print(f'drdt_inv shape: {drdt_inv.shape}  (nlayer={nlayer})\n')

    checks, drdt = validate_matrix(drdt_inv, nlayer)

    # Print results
    all_passed = True
    print(f'{"Check":<40}  {"Status":<6}  Detail')
    print('-' * 90)
    for name, passed, detail in checks:
        status = 'PASS' if passed else 'FAIL'
        if not passed:
            all_passed = False
        print(f'{name:<40}  {status:<6}  {detail}')

    print()
    if all_passed:
        print('Overall: PASS')
    else:
        print('Overall: FAIL')

    if args.plot:
        plot_path = os.path.join(OUTDIR, 'validate_planck_matrix.png')
        plot_matrix(drdt_inv, drdt, nlayer, plot_path)

    sys.exit(0 if all_passed else 1)


if __name__ == '__main__':
    main()
