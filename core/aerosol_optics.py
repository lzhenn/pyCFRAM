"""Aerosol optical property computation for CFRAM-A.

Translates Matlab scripts: main_aer.m, aer_opt_lw.m, aer_opt_sw.m,
aer_writer_lw.m, aer_writer_sw.m.

Computes AOD, SSA, and asymmetry parameter per RRTMG band from
MERRA-2 aerosol mixing ratios using GOCART lookup tables.
"""

import os
import numpy as np
from netCDF4 import Dataset


# Species configuration: (MERRA2 var, lookup table, radius/index)
# When rd < 1, it's treated as effective radius (m) for nearest-neighbor lookup.
# When rd >= 1, it's treated as a direct array index.
AEROSOL_SPECIES = {
    # Dust: 5 size bins, direct index
    'DU001': ('opticsBands_DU.v15_3.RRTMG.nc', 1),
    'DU002': ('opticsBands_DU.v15_3.RRTMG.nc', 2),
    'DU003': ('opticsBands_DU.v15_3.RRTMG.nc', 3),
    'DU004': ('opticsBands_DU.v15_3.RRTMG.nc', 4),
    'DU005': ('opticsBands_DU.v15_3.RRTMG.nc', 5),
    # Sea salt: 5 size bins, direct index
    'SS001': ('opticsBands_SS.v3_5.RRTMG.nc', 1),
    'SS002': ('opticsBands_SS.v3_5.RRTMG.nc', 2),
    'SS003': ('opticsBands_SS.v3_5.RRTMG.nc', 3),
    'SS004': ('opticsBands_SS.v3_5.RRTMG.nc', 4),
    'SS005': ('opticsBands_SS.v3_5.RRTMG.nc', 5),
    # Black carbon: phobic=bin1, philic=bin2
    'BCPHOBIC': ('opticsBands_BC.v1_3.RRTMG.nc', 1),
    'BCPHILIC': ('opticsBands_BC.v1_3.RRTMG.nc', 2),
    # Organic carbon: phobic=bin1, philic=bin2
    'OCPHOBIC': ('opticsBands_OC.v1_3.RRTMG.nc', 1),
    'OCPHILIC': ('opticsBands_OC.v1_3.RRTMG.nc', 2),
    # Sulfate: effective radius 0.16 um
    'SO4': ('opticsBands_SU.v1_3.RRTMG.nc', 0.16e-6),
}

# For SW: use v1_5 versions (30 bands)
AEROSOL_SPECIES_SW = {
    'DU001': ('opticsBands_DU.v15_5.RRTMG.nc', 1),
    'DU002': ('opticsBands_DU.v15_5.RRTMG.nc', 2),
    'DU003': ('opticsBands_DU.v15_5.RRTMG.nc', 3),
    'DU004': ('opticsBands_DU.v15_5.RRTMG.nc', 4),
    'DU005': ('opticsBands_DU.v15_5.RRTMG.nc', 5),
    'SS001': ('opticsBands_SS.v3_5.RRTMG.nc', 1),
    'SS002': ('opticsBands_SS.v3_5.RRTMG.nc', 2),
    'SS003': ('opticsBands_SS.v3_5.RRTMG.nc', 3),
    'SS004': ('opticsBands_SS.v3_5.RRTMG.nc', 4),
    'SS005': ('opticsBands_SS.v3_5.RRTMG.nc', 5),
    'BCPHOBIC': ('opticsBands_BC.v1_5.RRTMG.nc', 1),
    'BCPHILIC': ('opticsBands_BC.v1_5.RRTMG.nc', 2),
    'OCPHOBIC': ('opticsBands_OC.v1_5.RRTMG.nc', 1),
    'OCPHILIC': ('opticsBands_OC.v1_5.RRTMG.nc', 2),
    'SO4': ('opticsBands_SU.v1_5.RRTMG.nc', 0.16e-6),
}

# Band indices: SW bands 1-14 (0-indexed: 0-13), LW bands 15-30 (0-indexed: 14-29)
SW_BANDS = slice(0, 14)
LW_BANDS = slice(14, 30)


def _load_lookup_table(nc_path):
    """Load aerosol optical lookup table from NetCDF.

    Returns dict with bext, bsca, qext, qsca, g, rh, radius arrays.
    """
    nc = Dataset(nc_path, 'r')
    table = {}
    for vname in ('bext', 'bsca', 'qext', 'qsca', 'g', 'rh', 'radius'):
        if vname in nc.variables:
            table[vname] = np.array(nc.variables[vname][:], dtype=np.float64)
    nc.close()
    return table


def _find_radius_index(table, rd):
    """Find radius index: direct index if rd>=1, else nearest radius."""
    if rd >= 1:
        return int(rd) - 1  # Matlab 1-based to Python 0-based
    else:
        radius = table['radius']
        return np.argmin(np.abs(radius - rd))


def _find_rh_indices(table_rh, rh_profile):
    """Find nearest RH index for each layer."""
    indices = np.zeros(len(rh_profile), dtype=int)
    for i, rh_val in enumerate(rh_profile):
        indices[i] = np.argmin(np.abs(table_rh - rh_val))
    return indices


def compute_aer_optics_lw(mixing_ratio, table, rd, rh_profile, thick, dens):
    """Compute LW aerosol optical depth for one species.

    Args:
        mixing_ratio: (nlev,) mass mixing ratio (kg/kg)
        table: lookup table dict from _load_lookup_table
        rd: radius parameter (index or effective radius in m)
        rh_profile: (nlev,) relative humidity (0-1 fraction)
        thick: (nlev,) layer thickness (m)
        dens: (nlev,) air density (kg/m3)

    Returns:
        aod: (nbnd_lw, nlev) aerosol optical depth per LW band
    """
    nlev = len(mixing_ratio)
    r_idx = _find_radius_index(table, rd)

    # bext(band, rh, radius) - LW bands 15-30 (indices 14-29)
    bext = table['bext']
    lw_bext = bext[LW_BANDS, :, r_idx]  # (16, n_rh)

    rh_idx = _find_rh_indices(table['rh'], rh_profile)

    # Lookup kext for each layer by nearest RH
    kext = np.zeros((16, nlev), dtype=np.float64)
    for i in range(nlev):
        kext[:, i] = lw_bext[:, rh_idx[i]]

    # AOD = q * rho * 1e3 * dz * kext * 1e-3
    # q: kg/kg, rho: kg/m3, thick: m, kext: m2/kg
    aod = mixing_ratio[np.newaxis, :] * dens[np.newaxis, :] * 1e3 * \
          thick[np.newaxis, :] * kext * 1e-3

    aod = np.nan_to_num(aod, nan=0.0)
    return aod


def compute_aer_optics_sw(mixing_ratio, table, rd, rh_profile, thick, dens):
    """Compute SW aerosol optical properties for one species.

    Returns:
        aod: (nbnd_sw, nlev) optical depth
        ssa: (nbnd_sw, nlev) single scattering albedo
        g:   (nbnd_sw, nlev) asymmetry parameter
    """
    nlev = len(mixing_ratio)
    r_idx = _find_radius_index(table, rd)

    bext = table['bext'][SW_BANDS, :, r_idx]   # (14, n_rh)
    qext = table['qext'][SW_BANDS, :, r_idx]
    qsca = table['qsca'][SW_BANDS, :, r_idx]
    g_tab = table['g'][SW_BANDS, :, r_idx]

    rh_idx = _find_rh_indices(table['rh'], rh_profile)

    kext = np.zeros((14, nlev), dtype=np.float64)
    ssa = np.zeros((14, nlev), dtype=np.float64)
    g_out = np.zeros((14, nlev), dtype=np.float64)

    for i in range(nlev):
        kext[:, i] = bext[:, rh_idx[i]]
        qe = qext[:, rh_idx[i]]
        qs = qsca[:, rh_idx[i]]
        with np.errstate(divide='ignore', invalid='ignore'):
            ssa[:, i] = np.where(qe > 0, qs / qe, 0.0)
        g_out[:, i] = g_tab[:, rh_idx[i]]

    aod = mixing_ratio[np.newaxis, :] * dens[np.newaxis, :] * 1e3 * \
          thick[np.newaxis, :] * kext * 1e-3
    aod = np.nan_to_num(aod, nan=0.0)
    ssa = np.nan_to_num(ssa, nan=0.0)
    g_out = np.nan_to_num(g_out, nan=0.0)

    return aod, ssa, g_out


def compute_total_aerosol_optics(aerosol_data, rh_profile, delp, dens,
                                 lookup_dir):
    """Compute total aerosol optical properties for all species.

    Args:
        aerosol_data: dict mapping MERRA2 variable name -> (nlev,) mixing ratio
        rh_profile: (nlev,) relative humidity
        delp: (nlev,) pressure thickness (Pa)
        dens: (nlev,) air density (kg/m3)
        lookup_dir: directory containing opticsBands_*.nc files

    Returns:
        aod_lw: (nlev, nbndlw) total LW aerosol optical depth
        aod_sw: (nlev, jpband) total SW aerosol optical depth
        ssa_sw: (nlev, jpband) effective single scattering albedo
        g_sw:   (nlev, jpband) effective asymmetry parameter
    """
    nlev = len(rh_profile)
    thick = delp / (dens * 9.8)  # layer thickness in meters

    # Accumulate totals
    total_aod_lw = np.zeros((16, nlev), dtype=np.float64)
    total_aod_sw = np.zeros((14, nlev), dtype=np.float64)
    total_aod_ssa_sw = np.zeros((14, nlev), dtype=np.float64)  # for weighted SSA
    total_aod_g_sw = np.zeros((14, nlev), dtype=np.float64)    # for weighted g

    # LW lookup tables (_3 versions)
    _table_cache = {}

    for var_name, (lut_file, rd) in AEROSOL_SPECIES.items():
        if var_name not in aerosol_data:
            continue
        mixing = aerosol_data[var_name]

        lut_path = os.path.join(lookup_dir, lut_file)
        if lut_path not in _table_cache:
            _table_cache[lut_path] = _load_lookup_table(lut_path)
        table = _table_cache[lut_path]

        aod = compute_aer_optics_lw(mixing, table, rd, rh_profile, thick, dens)
        total_aod_lw += aod

    # SW lookup tables (_5 versions)
    for var_name, (lut_file, rd) in AEROSOL_SPECIES_SW.items():
        if var_name not in aerosol_data:
            continue
        mixing = aerosol_data[var_name]

        lut_path = os.path.join(lookup_dir, lut_file)
        if lut_path not in _table_cache:
            _table_cache[lut_path] = _load_lookup_table(lut_path)
        table = _table_cache[lut_path]

        aod, ssa, g = compute_aer_optics_sw(
            mixing, table, rd, rh_profile, thick, dens)
        total_aod_sw += aod
        total_aod_ssa_sw += aod * ssa    # AOD-weighted SSA
        total_aod_g_sw += aod * ssa * g  # AOD*SSA-weighted g

    # Effective SSA and g (AOD-weighted)
    with np.errstate(divide='ignore', invalid='ignore'):
        eff_ssa = np.where(total_aod_sw > 0,
                           total_aod_ssa_sw / total_aod_sw, 0.0)
        eff_g = np.where(total_aod_ssa_sw > 0,
                         total_aod_g_sw / total_aod_ssa_sw, 0.0)

    # Transpose to (nlev, nbnd) for Fortran compatibility
    return (total_aod_lw.T,     # (nlev, 16)
            total_aod_sw.T,     # (nlev, 14)
            eff_ssa.T,          # (nlev, 14)
            eff_g.T)            # (nlev, 14)
