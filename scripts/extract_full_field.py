#!/usr/bin/env python3
"""Extract full-field data from input NetCDF for Fortran CFRAM.

Reads base/perturbed state NetCDF files (per input_spec),
computes aerosol optical properties, writes Fortran binary files.

Usage:
    python3 scripts/extract_full_field.py --case eh13
"""
import os, sys, argparse
import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case, defaults, get_plev, get_aerosol_map, get_fortran_dir, get_lookup_dir
from core.constants import NBND_LW, NBND_SW

FORTRAN_PLEV = get_plev()
SCON = defaults()['radiation']['scon']
LOOKUP_DIR = get_lookup_dir()
FORTRAN_DIR = get_fortran_dir()

# Build aerosol map from config: {species: (lw_file, sw_file, rd)}
_aer_cfg = get_aerosol_map()
AEROSOL_MAP = {k: (v['lw'], v['sw'], v['rd']) for k, v in _aer_cfg.items()}

# Canonical species order (must match run_parallel_python.py and Fortran reader).
SPECIES_ORDER = ['bc', 'ocphi', 'ocpho', 'sulf', 'ss', 'dust']
NSPECIES = len(SPECIES_ORDER)


def write_bin(filepath, data):
    np.asarray(data, dtype=np.float64).tofile(filepath)


def load_lookup_table(nc_path):
    nc = Dataset(nc_path, 'r')
    table = {}
    for v in ('bext', 'bsca', 'qext', 'qsca', 'g', 'rh', 'radius'):
        if v in nc.variables:
            table[v] = np.array(nc.variables[v][:], dtype=np.float64)
    nc.close()
    return table


def find_radius_index(table, rd):
    if rd >= 1:
        return int(rd) - 1
    return int(np.argmin(np.abs(table['radius'] - rd)))


def compute_rh(q, t, p_hpa):
    w = q / (1.0 - q + 1e-30)
    e = w * p_hpa / (0.622 + w)
    es = 6.112 * np.exp(17.67 * (t - 273.15) / (t - 29.65 + 1e-10))
    return np.clip(e / (es + 1e-30), 0.0, 1.0)


def compute_aerosol_full_field(aer_data, rh, thick, dens, lookup_dir):
    """Compute aerosol optical properties for full field.

    Args:
        aer_data: dict {species: (nlev, nlat, nlon)} mixing ratios
        rh: (nlev, nlat, nlon) relative humidity
        thick: (nlev, nlat, nlon) layer thickness in meters
        dens: (nlev, nlat, nlon) air density kg/m3
    Returns:
        aod_lw: (nlev, nlat, nlon, 16)         bulk (sum of species)
        aod_sw: (nlev, nlat, nlon, 14)         bulk
        ssa_sw: (nlev, nlat, nlon, 14)         effective (AOD-weighted)
        g_sw:   (nlev, nlat, nlon, 14)         effective
        per_species: dict {species: {'aod_lw','aod_sw','ssa_sw','g_sw'}}
                     raw per-species arrays (ssa/g unweighted).
    """
    nlev, nlat, nlon = rh.shape
    total_aod_lw = np.zeros((nlev, nlat, nlon, NBND_LW))
    total_aod_sw = np.zeros((nlev, nlat, nlon, NBND_SW))
    total_aod_ssa_sw = np.zeros((nlev, nlat, nlon, NBND_SW))
    total_aod_g_sw = np.zeros((nlev, nlat, nlon, NBND_SW))
    per_species = {}

    LW_BANDS = slice(14, 30)
    SW_BANDS = slice(0, 14)
    _cache = {}

    for species, (lut_lw_file, lut_sw_file, rd) in AEROSOL_MAP.items():
        if species not in aer_data:
            continue
        mixing = aer_data[species]  # (nlev, nlat, nlon)
        print("  aerosol %s..." % species, end=' ', flush=True)

        # LW
        lut_lw_path = os.path.join(lookup_dir, lut_lw_file)
        if lut_lw_path not in _cache:
            _cache[lut_lw_path] = load_lookup_table(lut_lw_path)
        table = _cache[lut_lw_path]
        r_idx = find_radius_index(table, rd)
        # bext: (n_radius, n_rh, n_bands=30)
        bext_lw = table['bext'][r_idx, :, LW_BANDS]  # (n_rh, 16)
        rh_table = table['rh']
        n_rh = len(rh_table)

        # Vectorized RH index lookup
        rh_flat = rh.ravel()
        rh_idx = np.argmin(np.abs(rh_flat[:, None] - rh_table[None, :]), axis=1)
        # kext: (total_points, 16)
        kext_lw = bext_lw[rh_idx, :]  # (nlev*nlat*nlon, 16)
        kext_lw = kext_lw.reshape(nlev, nlat, nlon, NBND_LW)

        aod_sp_lw = mixing[:, :, :, None] * dens[:, :, :, None] * 1e3 * \
              thick[:, :, :, None] * kext_lw * 1e-3
        aod_sp_lw = np.nan_to_num(aod_sp_lw)
        total_aod_lw += aod_sp_lw

        # SW
        lut_sw_path = os.path.join(lookup_dir, lut_sw_file)
        if lut_sw_path not in _cache:
            _cache[lut_sw_path] = load_lookup_table(lut_sw_path)
        table_sw = _cache[lut_sw_path]
        r_idx_sw = find_radius_index(table_sw, rd)
        bext_sw = table_sw['bext'][r_idx_sw, :, SW_BANDS]
        qext_sw = table_sw['qext'][r_idx_sw, :, SW_BANDS]
        qsca_sw = table_sw['qsca'][r_idx_sw, :, SW_BANDS]
        g_sw_tab = table_sw['g'][r_idx_sw, :, SW_BANDS]

        rh_table_sw = table_sw['rh']
        rh_idx_sw = np.argmin(np.abs(rh_flat[:, None] - rh_table_sw[None, :]), axis=1)

        kext = bext_sw[rh_idx_sw, :].reshape(nlev, nlat, nlon, NBND_SW)
        qe = qext_sw[rh_idx_sw, :].reshape(nlev, nlat, nlon, NBND_SW)
        qs = qsca_sw[rh_idx_sw, :].reshape(nlev, nlat, nlon, NBND_SW)
        g_out = g_sw_tab[rh_idx_sw, :].reshape(nlev, nlat, nlon, NBND_SW)

        ssa = np.where(qe > 0, qs / qe, 0.0)
        aod_sp_sw = mixing[:, :, :, None] * dens[:, :, :, None] * 1e3 * \
                 thick[:, :, :, None] * kext * 1e-3
        aod_sp_sw = np.nan_to_num(aod_sp_sw)
        total_aod_sw += aod_sp_sw
        total_aod_ssa_sw += aod_sp_sw * ssa
        total_aod_g_sw += aod_sp_sw * ssa * g_out
        per_species[species] = {
            'aod_lw': aod_sp_lw,
            'aod_sw': aod_sp_sw,
            'ssa_sw': np.nan_to_num(ssa),
            'g_sw': np.nan_to_num(g_out),
        }
        print("done")

    with np.errstate(divide='ignore', invalid='ignore'):
        eff_ssa = np.where(total_aod_sw > 0, total_aod_ssa_sw / total_aod_sw, 0.0)
        eff_g = np.where(total_aod_ssa_sw > 0, total_aod_g_sw / total_aod_ssa_sw, 0.0)

    return total_aod_lw, total_aod_sw, eff_ssa, eff_g, per_species


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', required=True, help='Case name (directory under cases/)')
    args = parser.parse_args()

    cfg = load_case(args.case)
    outdir = os.path.join(FORTRAN_DIR, "data_prep")
    os.makedirs(outdir, exist_ok=True)
    print("=== Full-field extraction: %s ===" % cfg.get('case_name', args.case))

    # ---- Load input data (standard format per input_spec) ----
    nc_bp = Dataset(cfg['input']['base_pres'])
    nc_bs = Dataset(cfg['input']['base_surf'])
    nc_ap = Dataset(cfg['input']['perturbed_pres'])
    nc_as = Dataset(cfg['input']['perturbed_surf'])

    lats = np.array(nc_bp.variables['lat'][:])
    lons = np.array(nc_bp.variables['lon'][:])
    levs = np.array(nc_bp.variables['lev'][:])
    nlat, nlon, nlev = len(lats), len(lons), len(levs)
    print("Grid: %d lat x %d lon x %d lev" % (nlat, nlon, nlev))

    # Verify levels
    assert levs[0] > levs[-1], "Expected surface->TOA"
    assert np.allclose(levs[::-1], FORTRAN_PLEV), "Level mismatch"

    def get3d(nc, varname):
        """(time, lev, lat, lon) -> (nlev, nlat, nlon), flipped TOA->sfc."""
        return np.array(nc.variables[varname][0, ::-1, :, :], dtype=np.float64)

    def get2d(nc, varname):
        """(time, lat, lon) -> (nlat, nlon)."""
        return np.array(nc.variables[varname][0, :, :], dtype=np.float64)

    # ---- 3D atmospheric profiles (TOA->surface after flip) ----
    print("Loading 3D fields...")
    t_base = get3d(nc_bp, 'ta_lay')
    q_base = get3d(nc_bp, 'q')
    o3_base = get3d(nc_bp, 'o3')
    cc_base = get3d(nc_bp, 'camt')
    clwc_base = get3d(nc_bp, 'cliq')
    ciwc_base = get3d(nc_bp, 'cice')

    t_warm = get3d(nc_ap, 'ta_lay')
    q_warm = get3d(nc_ap, 'q')
    o3_warm = get3d(nc_ap, 'o3')
    cc_warm = get3d(nc_ap, 'camt')
    clwc_warm = get3d(nc_ap, 'cliq')
    ciwc_warm = get3d(nc_ap, 'cice')

    # Aerosol mixing ratios
    aer_species = ['bc', 'ocphi', 'ocpho', 'sulf', 'ss', 'dust']
    aer_base = {s: get3d(nc_bp, s) for s in aer_species}
    aer_warm = {s: get3d(nc_ap, s) for s in aer_species}

    # ---- 2D surface ----
    print("Loading 2D fields...")
    ts_base = get2d(nc_bs, 'ts')
    ps_base = get2d(nc_bs, 'ps')
    solar_base = get2d(nc_bs, 'solar')
    albedo_base = get2d(nc_bs, 'albedo')

    ts_warm = get2d(nc_as, 'ts')
    ps_warm = get2d(nc_as, 'ps')
    solar_warm = get2d(nc_as, 'solar')
    albedo_warm = get2d(nc_as, 'albedo')

    # CO2: volume mixing ratio -> ppmv (scalar mean)
    co2_base = get3d(nc_bp, 'co2').mean() * 1e6
    co2_warm = get3d(nc_ap, 'co2').mean() * 1e6
    print("CO2: base=%.2f, warm=%.2f ppmv" % (co2_base, co2_warm))

    nc_bp.close(); nc_bs.close(); nc_ap.close(); nc_as.close()

    # ---- Derived: ssrd, ssru ----
    ssrd_base = 300.0 * np.ones((nlat, nlon))
    ssru_base = ssrd_base * albedo_base
    ssrd_warm = 300.0 * np.ones((nlat, nlon))
    ssru_warm = ssrd_warm * albedo_warm

    # ---- Aerosol optical properties ----
    print("Computing aerosol optical properties...")
    p_hpa_3d = FORTRAN_PLEV[:, None, None] * np.ones((1, nlat, nlon))
    # Interface pressures
    p_half = np.zeros((nlev + 1, nlat, nlon))
    p_half[0] = 0.0
    for i in range(nlev - 1):
        p_half[i + 1] = 0.5 * (FORTRAN_PLEV[i] + FORTRAN_PLEV[i + 1])
    p_half[-1] = ps_base / 100.0
    delp = np.diff(p_half, axis=0) * 100.0  # Pa
    Rd = 287.05
    dens_base = (p_hpa_3d * 100.0) / (Rd * t_base)
    thick_base = delp / (dens_base * 9.8)
    rh_base = compute_rh(q_base, t_base, p_hpa_3d)

    if os.path.isdir(LOOKUP_DIR):
        aod_lw_base, aod_sw_base, ssa_sw_base, g_sw_base, perspc_base = \
            compute_aerosol_full_field(aer_base, rh_base, thick_base, dens_base, LOOKUP_DIR)

        # Warm state
        p_half_w = p_half.copy()
        p_half_w[-1] = ps_warm / 100.0
        delp_w = np.diff(p_half_w, axis=0) * 100.0
        dens_warm = (p_hpa_3d * 100.0) / (Rd * t_warm)
        thick_warm = delp_w / (dens_warm * 9.8)
        rh_warm = compute_rh(q_warm, t_warm, p_hpa_3d)
        aod_lw_warm, aod_sw_warm, ssa_sw_warm, g_sw_warm, perspc_warm = \
            compute_aerosol_full_field(aer_warm, rh_warm, thick_warm, dens_warm, LOOKUP_DIR)
    else:
        print("WARNING: no lookup tables, zero aerosol")
        aod_lw_base = np.zeros((nlev, nlat, nlon, NBND_LW))
        aod_sw_base = np.zeros((nlev, nlat, nlon, NBND_SW))
        ssa_sw_base = np.zeros((nlev, nlat, nlon, NBND_SW))
        g_sw_base = np.zeros((nlev, nlat, nlon, NBND_SW))
        aod_lw_warm = aod_lw_base.copy()
        aod_sw_warm = aod_sw_base.copy()
        ssa_sw_warm = ssa_sw_base.copy()
        g_sw_warm = g_sw_base.copy()
        perspc_base = {}
        perspc_warm = {}

    # ---- Write Fortran binary ----
    # Fortran reads: (ilev, ilat, ilon) for 3D, (ilat, ilon) for 2D
    # Fortran column-major: innermost loop is ilon, then ilat, then ilev
    # numpy default is C-order (row-major), tofile writes in C order
    # Fortran direct access reads by record, and the implied DO loops are:
    #   (((var(ilev,ilat,ilon), ilon=1,nlon), ilat=1,nlat), ilev=1,nlev)
    # This is column-major order: lon fastest, then lat, then lev
    # So we need to write in Fortran order: transpose to (nlev, nlat, nlon) and use order='F'
    # OR equivalently, the data is already (nlev, nlat, nlon) in C, but Fortran reads
    # ilon innermost -> we need to swap to (nlev, nlat, nlon) Fortran-order = (nlon, nlat, nlev) C-order
    # Simplest: use .astype(np.float64).flatten(order='F').tofile()
    # since Fortran traverses arrays in column-major order

    print("\nWriting Fortran binary files...")

    def write3d(fname, data):
        """Write (nlev, nlat, nlon) array in Fortran column-major order."""
        np.asfortranarray(data, dtype=np.float64).tofile(os.path.join(outdir, fname))
        sz = os.path.getsize(os.path.join(outdir, fname))
        print("  %s (%d bytes)" % (fname, sz))

    def write2d(fname, data):
        """Write (nlat, nlon) array in Fortran column-major order."""
        np.asfortranarray(data, dtype=np.float64).tofile(os.path.join(outdir, fname))
        sz = os.path.getsize(os.path.join(outdir, fname))
        print("  %s (%d bytes)" % (fname, sz))

    def write_scalar(fname, val):
        np.array([val], dtype=np.float64).tofile(os.path.join(outdir, fname))

    # 3D atmospheric
    write3d("t_base.dat", t_base)
    write3d("t_warm.dat", t_warm)
    write3d("hus_base.dat", q_base)
    write3d("hus_warm.dat", q_warm)
    write3d("O3_base.dat", o3_base)
    write3d("O3_warm.dat", o3_warm)
    write3d("cc_base.dat", cc_base)
    write3d("cc_warm.dat", cc_warm)
    write3d("clwc_base.dat", clwc_base)
    write3d("clwc_warm.dat", clwc_warm)
    write3d("ciwc_base.dat", ciwc_base)
    write3d("ciwc_warm.dat", ciwc_warm)

    # 2D surface
    write2d("skt_base.dat", ts_base)
    write2d("skt_warm.dat", ts_warm)
    write2d("sp_base.dat", ps_base)
    write2d("sp_warm.dat", ps_warm)
    write2d("solarin_base.dat", solar_base)
    write2d("solarin_warm.dat", solar_warm)
    write2d("ssrd_base.dat", ssrd_base)
    write2d("ssrd_warm.dat", ssrd_warm)
    write2d("ssru_base.dat", ssru_base)
    write2d("ssru_warm.dat", ssru_warm)

    # CO2 scalar
    write_scalar("co2_b.dat", co2_base)
    write_scalar("co2_w.dat", co2_warm)

    # Aerosol: (nlev, nlat, nlon, nbnd) -> Fortran reads
    # (((( var(ilev,ilat,ilon,ib), ib=1,nbnd), ilon=1,nlon), ilat=1,nlat), ilev=1,nlev)
    # = column-major with ib innermost
    def write_aer(fname, data):
        """Write (nlev, nlat, nlon, nbnd) in Fortran order."""
        np.asfortranarray(data, dtype=np.float64).tofile(os.path.join(outdir, fname))
        sz = os.path.getsize(os.path.join(outdir, fname))
        print("  %s (%d bytes)" % (fname, sz))

    write_aer("aerosol_aod_lw_base.dat", aod_lw_base)
    write_aer("aerosol_aod_lw_warm.dat", aod_lw_warm)
    write_aer("aerosol_aod_sw_base.dat", aod_sw_base)
    write_aer("aerosol_aod_sw_warm.dat", aod_sw_warm)
    write_aer("aerosol_ssa_sw_base.dat", ssa_sw_base)
    write_aer("aerosol_ssa_sw_warm.dat", ssa_sw_warm)
    write_aer("aerosol_g_sw_base.dat", g_sw_base)
    write_aer("aerosol_g_sw_warm.dat", g_sw_warm)

    # Per-species optical properties. Array shape (nlev, nlat, nlon, nbnd, nspecies)
    # passed through write_aer (np.asfortranarray): Fortran reads as
    # arr(nlev, nlat, nlon, nbnd, nspecies) with nlev fastest, species slowest --
    # consistent with bulk file layout with an added trailing species dim.
    def stack_full(perspc, field, nbnd):
        arr = np.zeros((nlev, nlat, nlon, nbnd, NSPECIES), dtype=np.float64)
        for i, s in enumerate(SPECIES_ORDER):
            if s in perspc:
                arr[..., i] = perspc[s][field]
        return arr

    write_aer("aerosol_aod_lw_base_spc.dat", stack_full(perspc_base, 'aod_lw', NBND_LW))
    write_aer("aerosol_aod_lw_warm_spc.dat", stack_full(perspc_warm, 'aod_lw', NBND_LW))
    write_aer("aerosol_aod_sw_base_spc.dat", stack_full(perspc_base, 'aod_sw', NBND_SW))
    write_aer("aerosol_aod_sw_warm_spc.dat", stack_full(perspc_warm, 'aod_sw', NBND_SW))
    write_aer("aerosol_ssa_sw_base_spc.dat", stack_full(perspc_base, 'ssa_sw', NBND_SW))
    write_aer("aerosol_ssa_sw_warm_spc.dat", stack_full(perspc_warm, 'ssa_sw', NBND_SW))
    write_aer("aerosol_g_sw_base_spc.dat", stack_full(perspc_base, 'g_sw', NBND_SW))
    write_aer("aerosol_g_sw_warm_spc.dat", stack_full(perspc_warm, 'g_sw', NBND_SW))

    # Verification
    print("\n=== Verification ===")
    expected_3d = nlev * nlat * nlon * 8
    expected_2d = nlat * nlon * 8
    expected_aer_lw = nlev * nlat * nlon * NBND_LW * 8
    expected_aer_sw = nlev * nlat * nlon * NBND_SW * 8
    checks = [
        ("t_base.dat", expected_3d),
        ("skt_base.dat", expected_2d),
        ("co2_b.dat", 8),
        ("aerosol_aod_lw_base.dat", expected_aer_lw),
        ("aerosol_aod_sw_base.dat", expected_aer_sw),
    ]
    all_ok = True
    for fn, exp in checks:
        actual = os.path.getsize(os.path.join(outdir, fn))
        ok = "OK" if actual == exp else "MISMATCH"
        if ok != "OK": all_ok = False
        print("  %s: %d bytes (expected %d) %s" % (fn, actual, exp, ok))
    print("nlat=%d, nlon=%d, nlev=%d" % (nlat, nlon, nlev))
    print("Ready for Fortran: modify nlat=%d, nlon=%d in cfram_rrtmg.f90" % (nlat, nlon))


if __name__ == '__main__':
    main()
