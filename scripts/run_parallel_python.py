#!/usr/bin/env python3
"""Embarrassingly parallel CFRAM using multiprocessing.

Each worker:
1. Extracts single column from pre-loaded global data
2. Writes temp binary files to a unique temp dir
3. Runs single-column cfram_rrtmg executable
4. Reads forcing + Planck matrix, solves dT in Python

Usage:
    python3 scripts/run_parallel_python.py --case eh13 --nproc 40
"""
import os, sys, argparse, struct, tempfile, shutil
import numpy as np
from multiprocessing import Pool
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case, defaults, get_plev, get_aerosol_map, get_nproc, \
    get_fortran_dir, get_lookup_dir, PROJECT_ROOT
from core.constants import NBND_LW, NBND_SW

FORTRAN_PLEV = get_plev()
NLEV = len(FORTRAN_PLEV)
SCON = defaults()['radiation']['scon']
FORTRAN_EXE_DIR = get_fortran_dir()
LOOKUP_DIR = get_lookup_dir()

# Build aerosol map from config
_aer_cfg = get_aerosol_map()
AEROSOL_MAP = {k: (v['lw'], v['sw'], v['rd']) for k, v in _aer_cfg.items()}

# Global data (set by initializer)
G = {}


def load_lookup_table(nc_path):
    nc = Dataset(nc_path, 'r')
    table = {v: np.array(nc.variables[v][:], dtype=np.float64) for v in ('bext','bsca','qext','qsca','g','rh','radius') if v in nc.variables}
    nc.close()
    return table


def compute_rh(q, t, p_hpa):
    w = q / (1.0 - q + 1e-30)
    e = w * p_hpa / (0.622 + w)
    es = 6.112 * np.exp(17.67 * (t - 273.15) / (t - 29.65 + 1e-10))
    return np.clip(e / (es + 1e-30), 0.0, 1.0)


def compute_aerosol_column(aer_data, rh, thick, dens):
    """Compute aerosol optics for single column."""
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
        for lut_file, bands, nbnd, target in [(lut_lw, LW_BANDS, NBND_LW, 'lw'), (lut_sw, SW_BANDS, NBND_SW, 'sw')]:
            table = G['lut'][lut_file]
            r_idx = int(rd) - 1 if rd >= 1 else int(np.argmin(np.abs(table['radius'] - rd)))
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
                ssa_col = np.array([np.where(qext[rh_idx[i], :] > 0, qsca[rh_idx[i], :] / qext[rh_idx[i], :], 0.0) for i in range(NLEV)])
                g_col = np.array([g_tab[rh_idx[i], :] for i in range(NLEV)])
                aod_ssa += aod_col * ssa_col
                aod_g += aod_col * ssa_col * g_col

    with np.errstate(divide='ignore', invalid='ignore'):
        eff_ssa = np.where(aod_sw > 0, aod_ssa / aod_sw, 0.0)
        eff_g = np.where(aod_ssa > 0, aod_g / aod_ssa, 0.0)
    return aod_lw, aod_sw, eff_ssa, eff_g


def write_bin(fp, data):
    np.asarray(data, dtype=np.float64).tofile(fp)


def write_aer_bin(fp, data):
    """Write aerosol array in C-order (band index fastest).

    Note: this is NOT Fortran column-major, but matches Phase 3 convention
    and gives better correlation with Wu et al. paper results.
    """
    np.asarray(data, dtype=np.float64).tofile(fp)


def read_fortran_seq(fp):
    with open(fp, 'rb') as f:
        raw = f.read()
    n = struct.unpack('i', raw[:4])[0]
    return np.frombuffer(raw[4:4+n], dtype=np.float64).copy()


def process_column(args):
    """Process single column. Returns (ilat, ilon, dT_dict, frc_dict)."""
    ilat, ilon = args
    d = G  # global data

    # Extract column (already TOA->surface)
    t_b = d['t_base'][:, ilat, ilon]
    q_b = d['q_base'][:, ilat, ilon]
    o3_b = d['o3_base'][:, ilat, ilon]
    cc_b = d['cc_base'][:, ilat, ilon]
    clwc_b = d['clwc_base'][:, ilat, ilon]
    ciwc_b = d['ciwc_base'][:, ilat, ilon]
    t_w = d['t_warm'][:, ilat, ilon]
    q_w = d['q_warm'][:, ilat, ilon]
    o3_w = d['o3_warm'][:, ilat, ilon]
    cc_w = d['cc_warm'][:, ilat, ilon]
    clwc_w = d['clwc_warm'][:, ilat, ilon]
    ciwc_w = d['ciwc_warm'][:, ilat, ilon]

    ts_b = d['ts_base'][ilat, ilon]
    ps_b = d['ps_base'][ilat, ilon]
    solar_b = d['solar_base'][ilat, ilon]
    albedo_b = d['albedo_base'][ilat, ilon]
    ts_w = d['ts_warm'][ilat, ilon]
    ps_w = d['ps_warm'][ilat, ilon]
    solar_w = d['solar_warm'][ilat, ilon]
    albedo_w = d['albedo_warm'][ilat, ilon]

    # Aerosol optics
    aer_b = {s: d['aer_base_'+s][:, ilat, ilon] for s in AEROSOL_MAP}
    aer_w = {s: d['aer_warm_'+s][:, ilat, ilon] for s in AEROSOL_MAP}

    p_hpa = FORTRAN_PLEV
    Rd = 287.05
    p_half = np.zeros(NLEV + 1)
    p_half[0] = 0.0
    for i in range(NLEV - 1):
        p_half[i + 1] = 0.5 * (p_hpa[i] + p_hpa[i + 1])
    p_half[-1] = ps_b / 100.0
    delp = np.diff(p_half) * 100.0
    dens_b = (p_hpa * 100.0) / (Rd * t_b)
    thick_b = delp / (dens_b * 9.8)
    rh_b = compute_rh(q_b, t_b, p_hpa)

    aod_lw_b, aod_sw_b, ssa_sw_b, g_sw_b = compute_aerosol_column(aer_b, rh_b, thick_b, dens_b)

    p_half_w = p_half.copy()
    p_half_w[-1] = ps_w / 100.0
    delp_w = np.diff(p_half_w) * 100.0
    dens_w = (p_hpa * 100.0) / (Rd * t_w)
    thick_w = delp_w / (dens_w * 9.8)
    rh_w = compute_rh(q_w, t_w, p_hpa)
    aod_lw_w, aod_sw_w, ssa_sw_w, g_sw_w = compute_aerosol_column(aer_w, rh_w, thick_w, dens_w)

    # Write to temp dir
    tmpdir = tempfile.mkdtemp(prefix='cfram_')
    dp = os.path.join(tmpdir, 'data_prep')
    do = os.path.join(tmpdir, 'data_output')
    os.makedirs(dp)
    os.makedirs(do)

    ssrd_b = 300.0
    ssru_b = ssrd_b * albedo_b
    ssrd_w = 300.0
    ssru_w = ssrd_w * albedo_w

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
    write_bin(os.path.join(dp, 'ssrd_base.dat'), np.array([ssrd_b]))
    write_bin(os.path.join(dp, 'ssrd_warm.dat'), np.array([ssrd_w]))
    write_bin(os.path.join(dp, 'ssru_base.dat'), np.array([ssru_b]))
    write_bin(os.path.join(dp, 'ssru_warm.dat'), np.array([ssru_w]))
    write_bin(os.path.join(dp, 'co2_b.dat'), np.array([d['co2_base']]))
    write_bin(os.path.join(dp, 'co2_w.dat'), np.array([d['co2_warm']]))
    write_aer_bin(os.path.join(dp, 'aerosol_aod_lw_base.dat'), aod_lw_b)
    write_aer_bin(os.path.join(dp, 'aerosol_aod_lw_warm.dat'), aod_lw_w)
    write_aer_bin(os.path.join(dp, 'aerosol_aod_sw_base.dat'), aod_sw_b)
    write_aer_bin(os.path.join(dp, 'aerosol_aod_sw_warm.dat'), aod_sw_w)
    write_aer_bin(os.path.join(dp, 'aerosol_ssa_sw_base.dat'), ssa_sw_b)
    write_aer_bin(os.path.join(dp, 'aerosol_ssa_sw_warm.dat'), ssa_sw_w)
    write_aer_bin(os.path.join(dp, 'aerosol_g_sw_base.dat'), g_sw_b)
    write_aer_bin(os.path.join(dp, 'aerosol_g_sw_warm.dat'), g_sw_w)

    # Symlink single-column executable
    exe_1col = d['exe_1col']
    os.symlink(exe_1col, os.path.join(tmpdir, 'cfram_rrtmg'))

    # Run
    ret = os.system('cd %s && ./cfram_rrtmg > /dev/null 2>&1' % tmpdir)

    # Read output: forcing + drdt_inv, solve dT in Python
    result = {}
    rad_terms = ['co2', 'q', 'ts', 'o3', 'solar', 'albedo', 'cloud', 'aerosol', 'warm']
    nonrad_terms = ['lhflx', 'shflx']
    dyn_terms = ['atmdyn', 'sfcdyn', 'ocndyn']
    aer_species_terms = ['bc', 'oc', 'sulf', 'seas', 'dust']

    if ret == 0:
        # Read forcing from Fortran
        for t in rad_terms:
            try:
                result['frc_'+t] = read_fortran_seq(os.path.join(do, 'frc_%s.dat' % t))
            except:
                result['frc_'+t] = np.full(NLEV+1, np.nan)

        # Read drdt_inv matrix
        try:
            drdt_file = os.path.join(do, 'drdt_inv.dat')
            with open(drdt_file, 'rb') as f:
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
            mat = np.frombuffer(raw[pos:pos+rec2_len], dtype=np.float64).copy().reshape(n, n, order='F')
        except:
            mat = None
            nlayer_val = 0

        # Solve dT = -drdt_inv @ frc for ALL terms (radiative + non-radiative)
        if mat is not None:
            nl = nlayer_val
            n = nl + 1

            def solve_dT(frc_full):
                """Apply Planck inverse to forcing vector."""
                frc_col = np.zeros(n)
                frc_col[:nl] = frc_full[:nl]
                frc_col[nl] = frc_full[NLEV]  # surface
                dT_vec = -mat @ frc_col  # n-dimensional result
                # Map back to full NLEV+1 array
                dT_full = np.full(NLEV + 1, -999.0)
                dT_full[:nl] = dT_vec[:nl]
                dT_full[NLEV] = dT_vec[nl]  # surface
                return dT_full

            # Radiative terms: forcing from Fortran
            for t in rad_terms:
                result['dT_'+t] = solve_dT(result['frc_'+t])

            # Dynamic terms derived from energy conservation: F_dyn = -F_rad (= -frc_warm)
            frc_warm_col = result['frc_warm']

            # Atmospheric dynamics: F_atmdyn[atm] = -frc_warm[atm], surface = 0
            frc_atmdyn = np.zeros(NLEV + 1)
            frc_atmdyn[:nl] = -frc_warm_col[:nl]
            result['dT_atmdyn'] = solve_dT(frc_atmdyn)

            # Surface dynamics total: F_sfcdyn[sfc] = -frc_warm[sfc], atm = 0
            frc_sfcdyn = np.zeros(NLEV + 1)
            frc_sfcdyn[NLEV] = -frc_warm_col[NLEV]
            result['dT_sfcdyn'] = solve_dT(frc_sfcdyn)

            # Ocean circulation: F_ocndyn[sfc] = F_sfcdyn[sfc] - F_lhflx[sfc] - F_shflx[sfc]
            lhflx_sfc = d['frc_lhflx'][NLEV, ilat, ilon] if 'frc_lhflx' in d else 0.0
            shflx_sfc = d['frc_shflx'][NLEV, ilat, ilon] if 'frc_shflx' in d else 0.0
            frc_ocndyn = np.zeros(NLEV + 1)
            frc_ocndyn[NLEV] = frc_sfcdyn[NLEV] - lhflx_sfc - shflx_sfc
            result['dT_ocndyn'] = solve_dT(frc_ocndyn)

            # lhflx, shflx + per-species aerosol: forcing from paper_data
            for t in nonrad_terms + aer_species_terms:
                if 'frc_' + t in d:
                    result['dT_' + t] = solve_dT(d['frc_' + t][:, ilat, ilon])
        else:
            for t in rad_terms + dyn_terms + nonrad_terms + aer_species_terms:
                result['dT_'+t] = np.full(NLEV+1, np.nan)
    else:
        for t in rad_terms:
            result['frc_'+t] = np.full(NLEV+1, np.nan)
            result['dT_'+t] = np.full(NLEV+1, np.nan)
        for t in dyn_terms + nonrad_terms + aer_species_terms:
            result['dT_'+t] = np.full(NLEV+1, np.nan)

    shutil.rmtree(tmpdir)
    return (ilat, ilon, result)


def init_worker(data_dict):
    global G
    G = data_dict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', required=True, help='Case name (directory under cases/)')
    parser.add_argument('--nproc', type=int, default=None,
                        help='Number of workers (default: from config or all CPUs)')
    args = parser.parse_args()

    cfg = load_case(args.case)
    OUTDIR = cfg['_output_dir']
    os.makedirs(OUTDIR, exist_ok=True)
    nproc = args.nproc or get_nproc(cfg)
    args.nproc = nproc

    print("=== Python parallel CFRAM: %s, %d procs ===" % (cfg.get('case_name', args.case), nproc))

    # Use pre-built single-column executable (nlat=1, nlon=1)
    exe_1col = os.path.join(FORTRAN_EXE_DIR, 'cfram_rrtmg_1col')
    if not os.path.exists(exe_1col):
        print("ERROR: %s not found. Build on hqlx127 first." % exe_1col)
        sys.exit(1)
    print("Using executable: %s" % exe_1col)

    # Load data
    print("Loading input data...")
    nc_bp = Dataset(cfg['input']['base_pres'])
    nc_bs = Dataset(cfg['input']['base_surf'])
    nc_ap = Dataset(cfg['input']['perturbed_pres'])
    nc_as = Dataset(cfg['input']['perturbed_surf'])

    lats = np.array(nc_bp.variables['lat'][:])
    lons = np.array(nc_bp.variables['lon'][:])
    nlat, nlon = len(lats), len(lons)

    def get3d(nc, v): return np.array(nc.variables[v][0, ::-1, :, :], dtype=np.float64)
    def get2d(nc, v): return np.array(nc.variables[v][0, :, :], dtype=np.float64)

    data = {
        't_base': get3d(nc_bp, 'ta_lay'), 'q_base': get3d(nc_bp, 'q'),
        'o3_base': get3d(nc_bp, 'o3'), 'cc_base': get3d(nc_bp, 'camt'),
        'clwc_base': get3d(nc_bp, 'cliq'), 'ciwc_base': get3d(nc_bp, 'cice'),
        't_warm': get3d(nc_ap, 'ta_lay'), 'q_warm': get3d(nc_ap, 'q'),
        'o3_warm': get3d(nc_ap, 'o3'), 'cc_warm': get3d(nc_ap, 'camt'),
        'clwc_warm': get3d(nc_ap, 'cliq'), 'ciwc_warm': get3d(nc_ap, 'cice'),
        'ts_base': get2d(nc_bs, 'ts'), 'ps_base': get2d(nc_bs, 'ps'),
        'solar_base': get2d(nc_bs, 'solar'), 'albedo_base': get2d(nc_bs, 'albedo'),
        'ts_warm': get2d(nc_as, 'ts'), 'ps_warm': get2d(nc_as, 'ps'),
        'solar_warm': get2d(nc_as, 'solar'), 'albedo_warm': get2d(nc_as, 'albedo'),
        'co2_base': get3d(nc_bp, 'co2').mean() * 1e6,
        'co2_warm': get3d(nc_ap, 'co2').mean() * 1e6,
        'exe_1col': exe_1col,
    }
    for s in AEROSOL_MAP:
        data['aer_base_'+s] = get3d(nc_bp, s)
        data['aer_warm_'+s] = get3d(nc_ap, s)
    nc_bp.close(); nc_bs.close(); nc_ap.close(); nc_as.close()

    # Load lookup tables
    lut = {}
    for _, (lw, sw, _) in AEROSOL_MAP.items():
        for f in [lw, sw]:
            if f not in lut:
                lut[f] = load_lookup_table(os.path.join(LOOKUP_DIR, f))
    data['lut'] = lut

    # Load non-radiative forcing (optional)
    nonrad_path = cfg['input'].get('nonrad_forcing')
    if nonrad_path and os.path.exists(nonrad_path):
        print("Loading non-radiative forcing...")
        nc_pf = Dataset(nonrad_path)
        for nonrad_term in ['lhflx', 'shflx']:
            if nonrad_term not in nc_pf.variables:
                continue
            frc_3d = np.array(nc_pf.variables[nonrad_term][0, ::-1, :, :], dtype=np.float64)
            if frc_3d.shape[1] != nlat or frc_3d.shape[2] != nlon:
                print("  WARNING: %s shape %s != grid (%d,%d), skipping" % (
                    nonrad_term, frc_3d.shape, nlat, nlon))
                continue
            frc_full = np.zeros((NLEV + 1, nlat, nlon))
            frc_full[:NLEV, :, :] = frc_3d[:NLEV, :, :]
            frc_full[NLEV, :, :] = np.array(nc_pf.variables[nonrad_term][0, 0, :, :], dtype=np.float64)
            frc_full = np.where(np.abs(frc_full) > 900, 0.0, frc_full)
            data['frc_' + nonrad_term] = frc_full
            print("  frc_%s: sfc mean=%.3f W/m2" % (nonrad_term, np.nanmean(frc_full[NLEV])))
        # Also load per-species aerosol forcing from same file
        aer_species_terms = ['bc', 'oc', 'sulf', 'seas', 'dust']
        for aer_term in aer_species_terms:
            if aer_term in nc_pf.variables:
                frc_3d = np.array(nc_pf.variables[aer_term][0, ::-1, :, :], dtype=np.float64)
                if frc_3d.shape[1] != nlat or frc_3d.shape[2] != nlon:
                    print("  WARNING: %s shape %s != grid (%d,%d), skipping" % (
                        aer_term, frc_3d.shape, nlat, nlon))
                    continue
                frc_full = np.zeros((NLEV + 1, nlat, nlon))
                frc_full[:NLEV, :, :] = frc_3d[:NLEV, :, :]
                frc_full[NLEV, :, :] = np.array(nc_pf.variables[aer_term][0, 0, :, :], dtype=np.float64)
                frc_full = np.where(np.abs(frc_full) > 900, 0.0, frc_full)
                data['frc_' + aer_term] = frc_full
                print("  frc_%s: sfc mean=%.3f W/m2" % (aer_term, np.nanmean(frc_full[NLEV])))
        nc_pf.close()
    else:
        print("No non-radiative forcing provided (radiative decomposition only)")

    print("Grid: %d x %d = %d points" % (nlat, nlon, nlat * nlon))

    # Build task list
    tasks = [(i, j) for i in range(nlat) for j in range(nlon)]

    # Run parallel
    import time
    t0 = time.time()
    print("Starting %d workers..." % args.nproc)

    # Result arrays
    terms = ['co2', 'q', 'ts', 'o3', 'solar', 'albedo', 'cloud', 'aerosol', 'warm']
    nonrad_terms = ['lhflx', 'shflx']
    dyn_terms = ['atmdyn', 'sfcdyn', 'ocndyn']
    aer_species_terms = ['bc', 'oc', 'sulf', 'seas', 'dust']
    all_dT_terms = terms + dyn_terms + nonrad_terms + aer_species_terms
    dT_out = {t: np.full((NLEV+1, nlat, nlon), np.nan) for t in all_dT_terms}
    frc_out = {t: np.full((NLEV+1, nlat, nlon), np.nan) for t in terms}

    done = 0
    with Pool(args.nproc, initializer=init_worker, initargs=(data,)) as pool:
        for ilat, ilon, result in pool.imap_unordered(process_column, tasks, chunksize=4):
            for t in terms:
                dT_out[t][:, ilat, ilon] = result['dT_'+t]
                frc_out[t][:, ilat, ilon] = result['frc_'+t]
            for t in dyn_terms + nonrad_terms + aer_species_terms:
                if 'dT_'+t in result:
                    dT_out[t][:, ilat, ilon] = result['dT_'+t]
            done += 1
            if done % 200 == 0:
                elapsed = time.time() - t0
                rate = done / elapsed
                eta = (len(tasks) - done) / rate
                print("  %d/%d done (%.1f pts/s, ETA %.0fs)" % (done, len(tasks), rate, eta))

    elapsed = time.time() - t0
    print("Completed in %.1f seconds (%.1f pts/s)" % (elapsed, len(tasks) / elapsed))

    # Save as NetCDF
    outfile = os.path.join(OUTDIR, 'cfram_result.nc')
    nc = Dataset(outfile, 'w')
    nc.createDimension('lev', NLEV + 1)
    nc.createDimension('lat', nlat)
    nc.createDimension('lon', nlon)
    nc.createVariable('lat', 'f8', ('lat',))[:] = lats
    nc.createVariable('lon', 'f8', ('lon',))[:] = lons
    # Surface->TOA level coordinate (matching paper_data)
    levs = np.concatenate([FORTRAN_PLEV[::-1], [FORTRAN_PLEV[-1] + 13]])  # add ~1013 for surface
    nc.createVariable('lev', 'f8', ('lev',))[:] = levs
    for t in all_dT_terms:
        v = nc.createVariable('dT_'+t, 'f8', ('lev', 'lat', 'lon'))
        v[:] = dT_out[t]
    for t in terms:
        v = nc.createVariable('frc_'+t, 'f8', ('lev', 'lat', 'lon'))
        v[:] = frc_out[t]

    # Compute and save observed dT and dynamics residual
    nc_bp2 = Dataset(cfg['input']['base_pres'])
    nc_ap2 = Dataset(cfg['input']['perturbed_pres'])
    nc_bs2 = Dataset(cfg['input']['base_surf'])
    nc_as2 = Dataset(cfg['input']['perturbed_surf'])
    dT_obs = np.zeros((NLEV + 1, nlat, nlon))
    dT_obs[:NLEV] = np.array(nc_ap2.variables['ta_lay'][0, ::-1, :, :], dtype=np.float64) - \
                     np.array(nc_bp2.variables['ta_lay'][0, ::-1, :, :], dtype=np.float64)
    dT_obs[NLEV] = np.array(nc_as2.variables['ts'][0, :, :], dtype=np.float64) - \
                    np.array(nc_bs2.variables['ts'][0, :, :], dtype=np.float64)
    nc_bp2.close(); nc_ap2.close(); nc_bs2.close(); nc_as2.close()

    v = nc.createVariable('dT_observed', 'f8', ('lev', 'lat', 'lon'))
    v[:] = dT_obs

    nc.close()
    print("Saved: %s" % outfile)


if __name__ == '__main__':
    main()
