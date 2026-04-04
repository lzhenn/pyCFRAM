"""Radiation driver: prepare atmospheric profiles for RRTMG and compute fluxes.

Translates the Fortran rad_driver.f90 logic into Python.
Handles coordinate transformation (ERA top-to-bottom → RRTMG bottom-to-top),
column density calculation, cloud water path, and flux convergence.

This module prepares inputs for RRTMG but does NOT call RRTMG directly —
the actual radiation call is done by the backend (subprocess or f2py).
"""

import numpy as np
from .constants import (GRAVITY, MW_DRY, MW_H2O, AVOGADRO,
                        CLOUD_RE_ICE, CLOUD_RE_LIQ, NBND_LW, NBND_SW)

# RRTMG band ranges (from parrrsw.f90)
JPB1 = 16  # first SW band index (1-based Fortran), 0-based = 15
JPB2 = 29  # last SW band index (1-based Fortran), 0-based = 28


def compute_interface_levels(plev, t, ps, ts):
    """Compute interface pressures and temperatures.

    Args:
        plev: (nlayer,) pressure levels, hPa, top-to-bottom (ERA order)
        t: (nlayer,) temperature at each level, K
        ps: surface pressure, hPa
        ts: surface temperature, K

    Returns:
        pint: (nlayer+1,) interface pressures, hPa, top-to-bottom
              pint[0] = 0.5 hPa (TOA), pint[nlayer] = ps
        tint: (nlayer+1,) interface temperatures, K
        dp: (nlayer,) pressure thickness of each layer, hPa
    """
    nlayer = len(plev)
    pint = np.zeros(nlayer + 1)
    tint = np.zeros(nlayer + 1)

    pint[0] = 0.5       # TOA boundary
    tint[0] = t[0]
    pint[nlayer] = ps
    tint[nlayer] = ts

    for i in range(1, nlayer):
        pint[i] = (plev[i - 1] + plev[i]) / 2.0
        tint[i] = (t[i - 1] + t[i]) / 2.0

    dp = np.diff(pint)  # pint[i+1] - pint[i], length nlayer
    return pint, tint, dp


def compute_column_densities(plev, t, q, o3, co2_ppmv, ch4_ppmv, n2o_ppmv,
                             pint, dp):
    """Compute gas mixing ratios and column densities.

    Args:
        plev: (nlayer,) pressure levels, hPa
        t: (nlayer,) temperature, K
        q: (nlayer,) specific humidity, kg/kg
        o3: (nlayer,) ozone mass mixing ratio, kg/kg
        co2_ppmv: CO2 concentration, ppmv
        ch4_ppmv: CH4 concentration, ppmv
        n2o_ppmv: N2O concentration, ppmv
        pint: (nlayer+1,) interface pressures, hPa
        dp: (nlayer,) pressure thickness, hPa

    Returns:
        wkl: (7, nlayer) volume mixing ratios (before multiplication by coldry)
        coldry: (nlayer,) dry air column density, molecules/cm²
        wbrodl: (nlayer,) broadening gas column density, molecules/cm²
    """
    nlayer = len(plev)
    wkl = np.zeros((7, nlayer))

    # Gas mixing ratios (volume)
    wkl[0, :] = q * 28.966 / 18.016          # H2O: specific humidity → volume mixing ratio
    wkl[1, :] = co2_ppmv * 1e-6              # CO2
    wkl[2, :] = o3 * 28.966 / 48.0           # O3: mass mixing ratio → volume
    wkl[3, :] = n2o_ppmv * 1e-6              # N2O
    # wkl[4, :] = 0  # CO (not used)
    wkl[5, :] = ch4_ppmv * 1e-6              # CH4
    wkl[6, :] = 0.209                         # O2

    coldry = np.zeros(nlayer)
    wbrodl = np.zeros(nlayer)

    for i in range(nlayer):
        # Moist air molecular weight
        amm = (1.0 - wkl[0, i]) * 28.966 + wkl[0, i] * 18.016
        # Column density: molecules/cm²
        coldry[i] = (dp[i] * 1e3 * 6.02e23) / (1e2 * 9.8 * amm * (1.0 + wkl[0, i]))
        # Broadening gas
        summol = wkl[1, i] + wkl[2, i] + wkl[3, i] + wkl[5, i] + wkl[6, i]
        wbrodl[i] = coldry[i] * (1.0 - summol)

    return wkl, coldry, wbrodl


def compute_cloud_water_path(cldfrac, cldiwc, cldlwc, dp):
    """Compute cloud ice/liquid water path from mixing ratios.

    Args:
        cldfrac: (nlayer,) cloud fraction
        cldiwc: (nlayer,) cloud ice water content, kg/kg
        cldlwc: (nlayer,) cloud liquid water content, kg/kg
        dp: (nlayer,) pressure thickness, hPa

    Returns:
        ciwp: (nlayer,) ice water path, g/m²
        clwp: (nlayer,) liquid water path, g/m²
    """
    coef = (1.0 / 9.8) * 1e2 * 1e3  # = 10204.08
    ciwp = np.where(cldfrac > 1e-5, coef * cldiwc * dp, 0.0)
    clwp = np.where(cldfrac > 1e-5, coef * cldlwc * dp, 0.0)
    return ciwp, clwp


def compute_pwvcm(coldry, wkl_col, pz_surface):
    """Compute precipitable water vapor in cm.

    Args:
        coldry: (nlayer,) dry air column density
        wkl_col: (7, nlayer) column amounts (coldry × mixing ratio)
        pz_surface: surface pressure in RRTMG order (= pz[0])

    Returns:
        pwvcm: precipitable water vapor, cm
    """
    amttl = np.sum(coldry + wkl_col[0, :])
    wvttl = np.sum(wkl_col[0, :])
    wvsh = (18.016 * wvttl) / (28.966 * amttl)
    pwvcm = wvsh * (1e3 * pz_surface) / (1e2 * 9.8066)
    return pwvcm


def flip_to_rrtmg(arr_era, nlayer):
    """Flip array from ERA order (top-to-bottom) to RRTMG order (bottom-to-top).

    For 1D layer arrays: rrtmg[i] = era[nlayer - i - 1]  (0-based)
    """
    return arr_era[:nlayer][::-1].copy()


def flip_interfaces_to_rrtmg(pint_era, tint_era, nlayer):
    """Flip interface arrays from ERA to RRTMG order.

    ERA:   pint[0]=TOA, pint[nlayer]=surface
    RRTMG: pz[0]=surface, pz[nlayer]=TOA
    """
    pz = pint_era[:nlayer + 1][::-1].copy()
    tz = tint_era[:nlayer + 1][::-1].copy()
    return pz, tz


def compute_flux_convergence(fd, fu, nlayer):
    """Compute net radiative flux convergence per layer.

    Convergence = energy absorbed by the layer = (flux_in_top + flux_in_bottom) - (flux_out_top + flux_out_bottom)

    In ERA order (index 1=TOA, nlayer+1=surface):
        conv[i] = fd[i] + fu[i+1] - fd[i+1] - fu[i]

    Args:
        fd: (nlayer+1,) downward flux, ERA order
        fu: (nlayer+1,) upward flux, ERA order
        nlayer: number of layers

    Returns:
        conv: (nlayer,) flux convergence per layer, W/m²
    """
    conv = np.zeros(nlayer)
    for i in range(nlayer):
        conv[i] = fd[i] + fu[i + 1] - fd[i + 1] - fu[i]
    return conv


def prepare_rrtmg_inputs(nlayer, plev, t, q, o3, co2_ppmv, ch4_ppmv, n2o_ppmv,
                          ps, ts, zenith, albedo_sw, albedo_lw,
                          cldfrac, cldlwc, cldiwc,
                          tauaer_sw, ssaaer_sw, asmaer_sw, tauaer_lw,
                          icld=2, iaer=10):
    """Prepare all RRTMG inputs from ERA-order atmospheric profiles.

    This replicates the full rad_driver.f90 preprocessing logic.

    All input arrays are in ERA order (top-to-bottom, index 0 = TOA).
    Returns a dict of arrays in RRTMG order (bottom-to-top, index 0 = surface).
    """
    # Interface levels (ERA order)
    pint, tint, dp = compute_interface_levels(plev[:nlayer], t[:nlayer], ps, ts)

    # Column densities (ERA order)
    wkl_mix, coldry_era, wbrodl_era = compute_column_densities(
        plev[:nlayer], t[:nlayer], q[:nlayer], o3[:nlayer],
        co2_ppmv, ch4_ppmv, n2o_ppmv, pint, dp)

    # Cloud water path (ERA order)
    ciwp_era, clwp_era = compute_cloud_water_path(
        cldfrac[:nlayer], cldiwc[:nlayer], cldlwc[:nlayer], dp)

    # Flip to RRTMG order (bottom-to-top)
    pave = flip_to_rrtmg(plev, nlayer)
    tave = flip_to_rrtmg(t, nlayer)
    pdp = flip_to_rrtmg(dp, nlayer)
    coldry = flip_to_rrtmg(coldry_era, nlayer)
    wbrodl = flip_to_rrtmg(wbrodl_era, nlayer)

    # wkl: flip each gas species
    wkl = np.zeros_like(wkl_mix)
    for imol in range(7):
        wkl[imol, :nlayer] = flip_to_rrtmg(wkl_mix[imol, :nlayer], nlayer)

    # Convert mixing ratios to column amounts
    for imol in range(7):
        wkl[imol, :nlayer] *= coldry[:nlayer]

    # Interface levels (RRTMG order)
    pz, tz = flip_interfaces_to_rrtmg(pint, tint, nlayer)

    # Precipitable water vapor
    pwvcm = compute_pwvcm(coldry[:nlayer], wkl[:, :nlayer], pz[0])

    # Surface emissivity
    semiss_lw = np.full(NBND_LW, 1.0 - albedo_lw)

    # Cloud properties (RRTMG order)
    cloud_frac = flip_to_rrtmg(cldfrac, nlayer) if icld != 0 else np.zeros(nlayer)
    ciwp = flip_to_rrtmg(ciwp_era, nlayer) if icld != 0 else np.zeros(nlayer)
    clwp = flip_to_rrtmg(clwp_era, nlayer) if icld != 0 else np.zeros(nlayer)
    rei = np.full(nlayer, CLOUD_RE_ICE)
    rel = np.full(nlayer, CLOUD_RE_LIQ)

    # Aerosol properties (RRTMG order)
    if iaer != 0 and tauaer_lw is not None:
        tau_lw = tauaer_lw[:nlayer][::-1].copy()
        tau_sw = tauaer_sw[:nlayer][::-1].copy()
        ssa_sw = ssaaer_sw[:nlayer][::-1].copy()
        asm_sw = asmaer_sw[:nlayer][::-1].copy()
    else:
        tau_lw = np.zeros((nlayer, NBND_LW))
        tau_sw = np.zeros((nlayer, NBND_SW))
        ssa_sw = np.zeros((nlayer, NBND_SW))
        asm_sw = np.zeros((nlayer, NBND_SW))

    return {
        'nlayer': nlayer,
        'iaer': iaer,
        'icld': icld,
        'pave': pave,
        'tave': tave,
        'pz': pz,
        'tz': tz,
        'pdp': pdp,
        'tbound': ts,
        'coldry': coldry,
        'wbrodl': wbrodl,
        'wkl': wkl,
        'pwvcm': pwvcm,
        'semiss_lw': semiss_lw,
        'zenith': zenith,
        'albedo_sw': albedo_sw,
        'cldfrac': cloud_frac,
        'ciwp': ciwp,
        'clwp': clwp,
        'rei': rei,
        'rel': rel,
        'tauaer_lw': tau_lw,
        'tauaer_sw': tau_sw,
        'ssaaer_sw': ssa_sw,
        'asmaer_sw': asm_sw,
    }
