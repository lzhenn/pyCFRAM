#!/usr/bin/env python3
"""Phase 2.1: Extract single column from paper_data and write Fortran binary.

Reads paper_data input_check.nc files (EH13), extracts one grid point
(115E, 32N), converts to Fortran cfram_rrtmg direct-access binary format.

Handles:
- Vertical axis flip (paper_data: surface->TOA, Fortran: TOA->surface)
- Variable name mapping (paper_data -> Fortran)
- Aerosol mixing ratio -> optical properties (via aerosol_optics.py)
- CO2 per-level -> scalar ppmv conversion
- albedo -> ssrd/ssru derivation

Run on hqlx204 (has netCDF4):
    python3 scripts/extract_paper_data_column.py
"""

import os
import sys
import numpy as np
from netCDF4 import Dataset

# ---- Config ----
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PAPER_DATA = os.path.join(PROJECT_ROOT, "paper_data", "cfram_out", "case_eh13_c20250102")
FORTRAN_DIR = os.path.join(PROJECT_ROOT, "fortran")
LOOKUP_DIR = os.path.join(PROJECT_ROOT, "fortran", "data_prep", "aerosol")
# Add core module path
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "core"))

TARGET_LAT = 32.0
TARGET_LON = 115.0
SCON = 1360.98  # solar constant, same as Fortran

# Fortran plev: TOA -> surface (37 levels, hPa)
FORTRAN_PLEV = np.array([1.,2.,3.,5.,7.,10.,20.,30.,50.,70.,100.,125.,150.,175.,200.,
    225.,250.,300.,350.,400.,450.,500.,550.,600.,650.,700.,750.,775.,
    800.,825.,850.,875.,900.,925.,950.,975.,1000.])

# MERRA-2 aerosol mapping: paper_data var -> list of (GOCART_bin, lookup_file, rd)
# paper_data has 6 aggregated species; we map to dominant GOCART bins
AEROSOL_MAP_LW = {
    'bc':    [('BCPHILIC', 'opticsBands_BC.v1_3.RRTMG.nc', 2)],
    'ocphi': [('OCPHILIC', 'opticsBands_OC.v1_3.RRTMG.nc', 2)],
    'ocpho': [('OCPHOBIC', 'opticsBands_OC.v1_3.RRTMG.nc', 1)],
    'sulf':  [('SO4',      'opticsBands_SU.v1_3.RRTMG.nc', 0.16e-6)],
    'ss':    [('SS003',    'opticsBands_SS.v3_5.RRTMG.nc', 3)],
    'dust':  [('DU003',    'opticsBands_DU.v15_3.RRTMG.nc', 3)],
}
AEROSOL_MAP_SW = {
    'bc':    [('BCPHILIC', 'opticsBands_BC.v1_5.RRTMG.nc', 2)],
    'ocphi': [('OCPHILIC', 'opticsBands_OC.v1_5.RRTMG.nc', 2)],
    'ocpho': [('OCPHOBIC', 'opticsBands_OC.v1_5.RRTMG.nc', 1)],
    'sulf':  [('SO4',      'opticsBands_SU.v1_5.RRTMG.nc', 0.16e-6)],
    'ss':    [('SS003',    'opticsBands_SS.v3_5.RRTMG.nc', 3)],
    'dust':  [('DU003',    'opticsBands_DU.v15_5.RRTMG.nc', 3)],
}

NBND_LW = 16
NBND_SW = 14


def find_nearest_idx(arr, val):
    return int(np.argmin(np.abs(np.array(arr) - val)))


def write_bin(filepath, data):
    """Write Fortran direct-access binary (raw float64)."""
    np.asarray(data, dtype=np.float64).tofile(filepath)
    print("  wrote %s  (%d bytes)" % (os.path.basename(filepath), os.path.getsize(filepath)))


def compute_rh(q, t, p_hpa):
    """Compute relative humidity (0-1) from specific humidity, T, and pressure.
    q: kg/kg, t: K, p_hpa: hPa
    """
    # Mixing ratio
    w = q / (1.0 - q + 1e-30)
    # Vapor pressure (hPa)
    e = w * p_hpa / (0.622 + w)
    # Saturation vapor pressure (Tetens formula, hPa)
    es = 6.112 * np.exp(17.67 * (t - 273.15) / (t - 29.65 + 1e-10))
    rh = np.clip(e / (es + 1e-30), 0.0, 1.0)
    return rh


def load_lookup_table(nc_path):
    """Load aerosol optical lookup table."""
    nc = Dataset(nc_path, 'r')
    table = {}
    for vname in ('bext', 'bsca', 'qext', 'qsca', 'g', 'rh', 'radius'):
        if vname in nc.variables:
            table[vname] = np.array(nc.variables[vname][:], dtype=np.float64)
    nc.close()
    return table


def find_radius_index(table, rd):
    if rd >= 1:
        return int(rd) - 1
    else:
        return int(np.argmin(np.abs(table['radius'] - rd)))


def find_rh_indices(table_rh, rh_profile):
    return np.array([np.argmin(np.abs(table_rh - rh_val)) for rh_val in rh_profile], dtype=int)


def compute_aerosol_optics_column(aerosol_mixing, rh, thick, dens, lookup_dir):
    """Compute per-band aerosol optical properties for single column.

    Args:
        aerosol_mixing: dict {paper_var: (nlev,) mixing ratio} in TOA->surface order
        rh: (nlev,) relative humidity
        thick: (nlev,) layer thickness in meters
        dens: (nlev,) air density kg/m3
        lookup_dir: path to opticsBands_*.nc files

    Returns:
        aod_lw: (nlev, 16) LW optical depth
        aod_sw: (nlev, 14) SW optical depth
        ssa_sw: (nlev, 14) SW single scattering albedo
        g_sw:   (nlev, 14) SW asymmetry parameter
    """
    nlev = len(rh)
    total_aod_lw = np.zeros((NBND_LW, nlev))
    total_aod_sw = np.zeros((NBND_SW, nlev))
    total_aod_ssa_sw = np.zeros((NBND_SW, nlev))
    total_aod_g_sw = np.zeros((NBND_SW, nlev))

    _cache = {}

    # Lookup table bext shape: (n_radius, n_rh, n_bands=30)
    # Bands 0-13: SW, Bands 14-29: LW
    LW_BANDS = slice(14, 30)
    SW_BANDS = slice(0, 14)

    for pvar, species_list in AEROSOL_MAP_LW.items():
        if pvar not in aerosol_mixing:
            continue
        mixing = aerosol_mixing[pvar]
        for _, lut_file, rd in species_list:
            lut_path = os.path.join(lookup_dir, lut_file)
            if lut_path not in _cache:
                _cache[lut_path] = load_lookup_table(lut_path)
            table = _cache[lut_path]
            r_idx = find_radius_index(table, rd)
            # bext(n_radius, n_rh, n_bands) -> select radius, then LW bands
            bext_lw = table['bext'][r_idx, :, LW_BANDS]  # (n_rh, 16)
            rh_idx = find_rh_indices(table['rh'], rh)
            kext = np.zeros((NBND_LW, nlev))
            for i in range(nlev):
                kext[:, i] = bext_lw[rh_idx[i], :]  # (16,)
            aod = mixing[np.newaxis, :] * dens[np.newaxis, :] * 1e3 * thick[np.newaxis, :] * kext * 1e-3
            total_aod_lw += np.nan_to_num(aod)

    for pvar, species_list in AEROSOL_MAP_SW.items():
        if pvar not in aerosol_mixing:
            continue
        mixing = aerosol_mixing[pvar]
        for _, lut_file, rd in species_list:
            lut_path = os.path.join(lookup_dir, lut_file)
            if lut_path not in _cache:
                _cache[lut_path] = load_lookup_table(lut_path)
            table = _cache[lut_path]
            r_idx = find_radius_index(table, rd)
            bext_sw = table['bext'][r_idx, :, SW_BANDS]    # (n_rh, 14)
            qext_sw = table['qext'][r_idx, :, SW_BANDS]
            qsca_sw = table['qsca'][r_idx, :, SW_BANDS]
            g_sw_tab = table['g'][r_idx, :, SW_BANDS]
            rh_idx = find_rh_indices(table['rh'], rh)
            kext = np.zeros((NBND_SW, nlev))
            ssa = np.zeros((NBND_SW, nlev))
            g_out = np.zeros((NBND_SW, nlev))
            for i in range(nlev):
                kext[:, i] = bext_sw[rh_idx[i], :]
                qe = qext_sw[rh_idx[i], :]
                qs = qsca_sw[rh_idx[i], :]
                ssa[:, i] = np.where(qe > 0, qs / qe, 0.0)
                g_out[:, i] = g_sw_tab[rh_idx[i], :]
            aod = mixing[np.newaxis, :] * dens[np.newaxis, :] * 1e3 * thick[np.newaxis, :] * kext * 1e-3
            aod = np.nan_to_num(aod)
            total_aod_sw += aod
            total_aod_ssa_sw += aod * ssa
            total_aod_g_sw += aod * ssa * g_out

    with np.errstate(divide='ignore', invalid='ignore'):
        eff_ssa = np.where(total_aod_sw > 0, total_aod_ssa_sw / total_aod_sw, 0.0)
        eff_g = np.where(total_aod_ssa_sw > 0, total_aod_g_sw / total_aod_ssa_sw, 0.0)

    return total_aod_lw.T, total_aod_sw.T, eff_ssa.T, eff_g.T


def main():
    print("=== Phase 2.1: Extract paper_data column for Fortran CFRAM ===\n")

    # ---- 1. Read paper_data ----
    files = os.listdir(PAPER_DATA)
    f_base_pres = os.path.join(PAPER_DATA, [f for f in files if 'baseline_pres' in f][0])
    f_base_surf = os.path.join(PAPER_DATA, [f for f in files if 'baseline_surf' in f][0])
    f_all_pres = os.path.join(PAPER_DATA, [f for f in files if 'all_pres' in f][0])
    f_all_surf = os.path.join(PAPER_DATA, [f for f in files if 'all_surf' in f][0])

    nc_bp = Dataset(f_base_pres)
    nc_bs = Dataset(f_base_surf)
    nc_ap = Dataset(f_all_pres)
    nc_as = Dataset(f_all_surf)

    # Get coordinates
    lats = np.array(nc_bp.variables['lat'][:])
    lons = np.array(nc_bp.variables['lon'][:])
    levs = np.array(nc_bp.variables['lev'][:])  # surface -> TOA (1000, 975, ..., 1)

    lat_idx = find_nearest_idx(lats, TARGET_LAT)
    lon_idx = find_nearest_idx(lons, TARGET_LON)
    print("Grid: lat %.2f-%.2f (%d), lon %.2f-%.2f (%d)" %
          (lats.min(), lats.max(), len(lats), lons.min(), lons.max(), len(lons)))
    print("Target: %.1fE, %.1fN -> idx (%d, %d) = (%.2fN, %.2fE)" %
          (TARGET_LON, TARGET_LAT, lat_idx, lon_idx, lats[lat_idx], lons[lon_idx]))
    print("Levels: %s (surface->TOA, %d pts)" % (levs[:5], len(levs)))

    # Verify level ordering matches Fortran (after flip)
    assert levs[0] > levs[-1], "Expected paper_data levels: surface->TOA"
    fortran_levs_check = levs[::-1]  # flip to TOA->surface
    assert np.allclose(fortran_levs_check, FORTRAN_PLEV), \
        "Level mismatch!\npaper(flipped): %s\nFortran: %s" % (fortran_levs_check, FORTRAN_PLEV)
    print("Level check PASSED: paper_data levels (flipped) match Fortran plev")

    # ---- 2. Extract column data ----
    # Helper: extract (time=0, :, lat_idx, lon_idx) and flip vertical
    def get3d(nc, varname):
        """Extract 3D variable, flip to TOA->surface."""
        data = np.array(nc.variables[varname][0, :, lat_idx, lon_idx], dtype=np.float64)
        return data[::-1]  # flip: surface->TOA => TOA->surface

    def get2d(nc, varname):
        """Extract 2D surface variable."""
        return float(nc.variables[varname][0, lat_idx, lon_idx])

    # Atmospheric profiles (37 levels, TOA->surface after flip)
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

    # Surface variables
    ts_base = get2d(nc_bs, 'ts')
    ps_base = get2d(nc_bs, 'ps')  # Pa
    solar_base = get2d(nc_bs, 'solar')  # W/m2
    albedo_base = get2d(nc_bs, 'albedo')
    emiss_base = get2d(nc_bs, 'emiss')

    ts_warm = get2d(nc_as, 'ts')
    ps_warm = get2d(nc_as, 'ps')
    solar_warm = get2d(nc_as, 'solar')
    albedo_warm = get2d(nc_as, 'albedo')
    emiss_warm = get2d(nc_as, 'emiss')

    # CO2: per-level in paper_data, scalar ppmv in Fortran
    # Check units by magnitude
    co2_profile_base = get3d(nc_bp, 'co2')
    co2_profile_warm = get3d(nc_ap, 'co2')
    co2_mean_base = co2_profile_base.mean()
    co2_mean_warm = co2_profile_warm.mean()

    if co2_mean_base < 1.0:
        # Volume mixing ratio (mol/mol) -> ppmv
        # MERRA-2 co2 is typically in volume mixing ratio
        co2_ppmv_base = co2_mean_base * 1e6
        co2_ppmv_warm = co2_mean_warm * 1e6
        # Sanity check: should be ~390-420 ppmv for 2003-2022
        if co2_ppmv_base < 300 or co2_ppmv_base > 500:
            print("\nWARNING: CO2 = %.2f ppmv seems unreasonable, check units!" % co2_ppmv_base)
        print("\nCO2 conversion: vmr %.6e -> %.2f ppmv (base)" % (co2_mean_base, co2_ppmv_base))
    else:
        # Already ppmv
        co2_ppmv_base = co2_mean_base
        co2_ppmv_warm = co2_mean_warm
        print("\nCO2 already ppmv: %.2f (base), %.2f (warm)" % (co2_ppmv_base, co2_ppmv_warm))

    # Aerosol mixing ratios
    aer_species = ['bc', 'ocphi', 'ocpho', 'sulf', 'ss', 'dust']
    aer_base = {s: get3d(nc_bp, s) for s in aer_species}
    aer_warm = {s: get3d(nc_ap, s) for s in aer_species}

    nc_bp.close()
    nc_bs.close()
    nc_ap.close()
    nc_as.close()

    # ---- 3. Diagnostics ----
    print("\n--- Extracted Column Diagnostics ---")
    print("T_base: [%.2f, %.2f] K (TOA->sfc)" % (t_base[0], t_base[-1]))
    print("T_warm: [%.2f, %.2f] K" % (t_warm[0], t_warm[-1]))
    print("q_base: [%.2e, %.2e] kg/kg" % (q_base.min(), q_base.max()))
    print("o3_base: [%.2e, %.2e] kg/kg" % (o3_base.min(), o3_base.max()))
    print("cc_base: [%.4f, %.4f]" % (cc_base.min(), cc_base.max()))
    print("clwc_base: [%.2e, %.2e]" % (clwc_base.min(), clwc_base.max()))
    print("ciwc_base: [%.2e, %.2e]" % (ciwc_base.min(), ciwc_base.max()))
    print("ts: base=%.2f, warm=%.2f K" % (ts_base, ts_warm))
    print("ps: base=%.1f, warm=%.1f Pa" % (ps_base, ps_warm))
    print("solar: base=%.2f, warm=%.2f W/m2" % (solar_base, solar_warm))
    print("albedo: base=%.4f, warm=%.4f" % (albedo_base, albedo_warm))
    print("emiss: base=%.4f, warm=%.4f" % (emiss_base, emiss_warm))
    print("NOTE: Fortran hardcodes albedo_lw=0 (emiss=1), paper_data has emiss=%.4f" % emiss_base)
    print("co2: base=%.2f, warm=%.2f ppmv" % (co2_ppmv_base, co2_ppmv_warm))
    for s in aer_species:
        print("aer_%s: base=[%.2e, %.2e]" % (s, aer_base[s].min(), aer_base[s].max()))

    # ---- 4. Derived quantities ----
    # Fortran computes: albedo_sw = ssru/ssrd, zenith = solin/scon
    # We have: solar (=solin), albedo (=albedo_sw)
    # So: solin = solar, ssrd = arbitrary nonzero, ssru = ssrd * albedo
    solarin_base = solar_base
    solarin_warm = solar_warm
    # Use physically plausible ssrd (doesn't affect calculation, only ratio matters)
    ssrd_base = 300.0  # approximate surface SW down
    ssru_base = ssrd_base * albedo_base
    ssrd_warm = 300.0
    ssru_warm = ssrd_warm * albedo_warm

    print("\nDerived: ssrd=%.1f, ssru=%.4f (albedo_sw=%.4f)" %
          (ssrd_base, ssru_base, ssru_base / ssrd_base))

    # ---- 5. Compute aerosol optical properties ----
    print("\n--- Computing Aerosol Optical Properties ---")
    # Need: RH, layer thickness, air density
    p_hpa = FORTRAN_PLEV  # pressure levels in hPa (TOA->surface)

    # Interface pressures (half levels): p_half[0]=TOA top, p_half[nlev]=surface
    p_half = np.zeros(len(p_hpa) + 1)
    p_half[0] = 0.0  # top of atmosphere
    for i in range(len(p_hpa) - 1):
        p_half[i + 1] = 0.5 * (p_hpa[i] + p_hpa[i + 1])
    p_half[-1] = ps_base / 100.0  # surface pressure in hPa

    delp = np.diff(p_half) * 100.0  # Pa
    Rd = 287.05  # J/(kg*K)
    dens = (p_hpa * 100.0) / (Rd * t_base)  # kg/m3
    thick = delp / (dens * 9.8)  # m

    rh_base = compute_rh(q_base, t_base, p_hpa)
    rh_warm = compute_rh(q_warm, t_warm, p_hpa)

    if os.path.isdir(LOOKUP_DIR):
        aod_lw_base, aod_sw_base, ssa_sw_base, g_sw_base = \
            compute_aerosol_optics_column(aer_base, rh_base, thick, dens, LOOKUP_DIR)
        # Recompute for warm state
        dens_w = (p_hpa * 100.0) / (Rd * t_warm)
        p_half_w = p_half.copy()
        p_half_w[-1] = ps_warm / 100.0
        delp_w = np.diff(p_half_w) * 100.0
        thick_w = delp_w / (dens_w * 9.8)
        aod_lw_warm, aod_sw_warm, ssa_sw_warm, g_sw_warm = \
            compute_aerosol_optics_column(aer_warm, rh_warm, thick_w, dens_w, LOOKUP_DIR)
        print("  AOD_LW total: base=%.4e, warm=%.4e" %
              (aod_lw_base.sum(), aod_lw_warm.sum()))
        print("  AOD_SW total: base=%.4e, warm=%.4e" %
              (aod_sw_base.sum(), aod_sw_warm.sum()))
    else:
        print("  WARNING: lookup dir not found, using zero aerosol")
        aod_lw_base = np.zeros((37, NBND_LW))
        aod_sw_base = np.zeros((37, NBND_SW))
        ssa_sw_base = np.zeros((37, NBND_SW))
        g_sw_base = np.zeros((37, NBND_SW))
        aod_lw_warm = np.zeros((37, NBND_LW))
        aod_sw_warm = np.zeros((37, NBND_SW))
        ssa_sw_warm = np.zeros((37, NBND_SW))
        g_sw_warm = np.zeros((37, NBND_SW))

    # ---- 6. Write Fortran binary files ----
    outdir = os.path.join(FORTRAN_DIR, "data_prep")
    os.makedirs(outdir, exist_ok=True)
    print("\n--- Writing Fortran Binary Files to %s ---" % outdir)

    # 3D atmospheric (37 values each, TOA->surface)
    write_bin(os.path.join(outdir, "t_base.dat"), t_base)
    write_bin(os.path.join(outdir, "t_warm.dat"), t_warm)
    write_bin(os.path.join(outdir, "hus_base.dat"), q_base)
    write_bin(os.path.join(outdir, "hus_warm.dat"), q_warm)
    write_bin(os.path.join(outdir, "O3_base.dat"), o3_base)
    write_bin(os.path.join(outdir, "O3_warm.dat"), o3_warm)
    write_bin(os.path.join(outdir, "cc_base.dat"), cc_base)
    write_bin(os.path.join(outdir, "cc_warm.dat"), cc_warm)
    write_bin(os.path.join(outdir, "clwc_base.dat"), clwc_base)
    write_bin(os.path.join(outdir, "clwc_warm.dat"), clwc_warm)
    write_bin(os.path.join(outdir, "ciwc_base.dat"), ciwc_base)
    write_bin(os.path.join(outdir, "ciwc_warm.dat"), ciwc_warm)

    # 2D surface scalars (1 value each)
    write_bin(os.path.join(outdir, "skt_base.dat"), np.array([ts_base]))
    write_bin(os.path.join(outdir, "skt_warm.dat"), np.array([ts_warm]))
    write_bin(os.path.join(outdir, "sp_base.dat"), np.array([ps_base]))
    write_bin(os.path.join(outdir, "sp_warm.dat"), np.array([ps_warm]))
    write_bin(os.path.join(outdir, "solarin_base.dat"), np.array([solarin_base]))
    write_bin(os.path.join(outdir, "solarin_warm.dat"), np.array([solarin_warm]))
    write_bin(os.path.join(outdir, "ssrd_base.dat"), np.array([ssrd_base]))
    write_bin(os.path.join(outdir, "ssrd_warm.dat"), np.array([ssrd_warm]))
    write_bin(os.path.join(outdir, "ssru_base.dat"), np.array([ssru_base]))
    write_bin(os.path.join(outdir, "ssru_warm.dat"), np.array([ssru_warm]))

    # CO2 scalar
    write_bin(os.path.join(outdir, "co2_b.dat"), np.array([co2_ppmv_base]))
    write_bin(os.path.join(outdir, "co2_w.dat"), np.array([co2_ppmv_warm]))

    # Aerosol optical properties
    # LW: (nlev, nbndlw) = (37, 16)
    write_bin(os.path.join(outdir, "aerosol_aod_lw_base.dat"), aod_lw_base)
    write_bin(os.path.join(outdir, "aerosol_aod_lw_warm.dat"), aod_lw_warm)
    # SW: (nlev, nbnd_sw) = (37, 14)
    write_bin(os.path.join(outdir, "aerosol_aod_sw_base.dat"), aod_sw_base)
    write_bin(os.path.join(outdir, "aerosol_aod_sw_warm.dat"), aod_sw_warm)
    write_bin(os.path.join(outdir, "aerosol_ssa_sw_base.dat"), ssa_sw_base)
    write_bin(os.path.join(outdir, "aerosol_ssa_sw_warm.dat"), ssa_sw_warm)
    write_bin(os.path.join(outdir, "aerosol_g_sw_base.dat"), g_sw_base)
    write_bin(os.path.join(outdir, "aerosol_g_sw_warm.dat"), g_sw_warm)

    # ---- 7. Verification summary ----
    print("\n=== Verification Summary ===")
    print("Files written: %d" % len([f for f in os.listdir(outdir) if f.endswith('.dat')]))
    expected_sizes = {
        't_base.dat': 296, 'hus_base.dat': 296, 'O3_base.dat': 296,
        'cc_base.dat': 296, 'clwc_base.dat': 296, 'ciwc_base.dat': 296,
        'skt_base.dat': 8, 'sp_base.dat': 8, 'solarin_base.dat': 8,
        'ssrd_base.dat': 8, 'ssru_base.dat': 8, 'co2_b.dat': 8,
        'aerosol_aod_lw_base.dat': 4736, 'aerosol_aod_sw_base.dat': 4144,
        'aerosol_ssa_sw_base.dat': 4144, 'aerosol_g_sw_base.dat': 4144,
    }
    all_ok = True
    for fname, expected in expected_sizes.items():
        fpath = os.path.join(outdir, fname)
        actual = os.path.getsize(fpath)
        if actual != expected:
            print("  SIZE MISMATCH: %s expected %d got %d" % (fname, expected, actual))
            all_ok = False
    if all_ok:
        print("  All file sizes match expected values!")

    print("\nReady to run Fortran CFRAM on hqlx74:")
    print("  cd %s && ./cfram_rrtmg" % FORTRAN_DIR)


if __name__ == '__main__':
    main()
