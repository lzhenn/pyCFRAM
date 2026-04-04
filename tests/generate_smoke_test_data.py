"""Generate synthetic input data for smoke testing cfram_rrtmg.

Creates a single-column atmosphere based on US Standard Atmosphere 1976,
writes all required binary input files to a data_prep/ directory,
then runs the Fortran executable.

Usage:
    python tests/generate_smoke_test_data.py [--exe PATH] [--outdir PATH]
"""

import os
import sys
import numpy as np
import argparse

# Add parent dir to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from core.constants import ERA5_PLEVELS, NLEV, NBND_LW, NBND_SW


def us_standard_temperature(p_hpa):
    """Approximate US Standard Atmosphere 1976 temperature profile.

    Args:
        p_hpa: pressure in hPa (array)

    Returns:
        temperature in K (array)
    """
    # Simple piecewise linear approximation
    t = np.zeros_like(p_hpa, dtype=np.float64)
    for i, p in enumerate(p_hpa):
        if p >= 226:   # troposphere (surface to ~11km)
            t[i] = 288.15 - 6.5 * (1 - (p / 1013.25) ** 0.19) * 44.33
            # Simplified: T decreases from 288K at surface
            t[i] = max(t[i], 216.65)
        elif p >= 54:  # lower stratosphere (11-20km), isothermal
            t[i] = 216.65
        elif p >= 8.7:  # upper stratosphere (20-32km)
            z_approx = 44.33 * (1 - (p / 1013.25) ** 0.19)
            t[i] = 216.65 + 1.0 * (z_approx - 20.0)
            t[i] = min(t[i], 270.65)
        else:  # mesosphere
            t[i] = 270.65 - 2.8 * (44.33 * (1 - (p / 1013.25) ** 0.19) - 32.0)
            t[i] = max(t[i], 186.87)
    return t


def us_standard_humidity(p_hpa, t_k):
    """Approximate specific humidity profile.

    Roughly exponential decrease with altitude, capped by saturation.
    """
    # Simple exponential decrease from surface
    q = np.zeros_like(p_hpa, dtype=np.float64)
    q_surface = 0.008  # ~8 g/kg at surface
    for i, p in enumerate(p_hpa):
        q[i] = q_surface * (p / 1013.25) ** 3
        q[i] = max(q[i], 3e-6)  # minimum stratospheric value
    return q


def generate_test_data(outdir):
    """Generate all input binary files for cfram_rrtmg single-column test."""
    os.makedirs(outdir, exist_ok=True)

    nlev = NLEV
    nlat = 1
    nlon = 1

    # Pressure levels (top to bottom as in ERA5)
    plev = ERA5_PLEVELS  # 1, 2, 3, ... 1000 hPa

    # Temperature profile
    t_base = us_standard_temperature(plev)
    # Warm case: +2K uniform warming
    t_warm = t_base + 2.0

    # Humidity profile
    q_base = us_standard_humidity(plev, t_base)
    # Warm: +10% humidity (Clausius-Clapeyron-like)
    q_warm = q_base * 1.10

    # Ozone (simple stratospheric peak, kg/kg mass mixing ratio)
    o3_base = np.zeros(nlev, dtype=np.float64)
    for i, p in enumerate(plev):
        if 10 < p < 100:
            o3_base[i] = 8e-6 * np.exp(-((np.log(p) - np.log(30)) ** 2) / 2.0)
        elif p <= 10:
            o3_base[i] = 2e-6
        else:
            o3_base[i] = 3e-8
    o3_warm = o3_base * 0.98  # slight decrease

    # Surface quantities
    ts_base = np.array([288.15], dtype=np.float64)  # surface temperature
    ts_warm = np.array([290.15], dtype=np.float64)

    ps_base = np.array([101325.0], dtype=np.float64)  # Pa (Fortran divides by 100)
    ps_warm = np.array([101325.0], dtype=np.float64)

    # TOA solar irradiance (W/m2, daily mean at ~30N summer)
    solin_base = np.array([450.0], dtype=np.float64)
    solin_warm = np.array([450.0], dtype=np.float64)

    # Surface SW down/up (W/m2)
    ssrd_base = np.array([250.0], dtype=np.float64)
    ssrd_warm = np.array([248.0], dtype=np.float64)
    ssru_base = np.array([50.0], dtype=np.float64)  # albedo ~0.2
    ssru_warm = np.array([49.6], dtype=np.float64)

    # CO2 (ppmv)
    co2_base = np.array([380.0], dtype=np.float64)
    co2_warm = np.array([400.0], dtype=np.float64)

    # Cloud properties (clear sky for smoke test)
    cc_base = np.zeros(nlev, dtype=np.float64)
    cc_warm = np.zeros(nlev, dtype=np.float64)
    clwc_base = np.zeros(nlev, dtype=np.float64)
    clwc_warm = np.zeros(nlev, dtype=np.float64)
    ciwc_base = np.zeros(nlev, dtype=np.float64)
    ciwc_warm = np.zeros(nlev, dtype=np.float64)

    # Aerosol optical properties (zero for smoke test)
    nbndlw = NBND_LW  # 16
    jpband = NBND_SW   # 14
    aod_lw = np.zeros((nlev, nlat, nlon, nbndlw), dtype=np.float64)
    aod_sw = np.zeros((nlev, nlat, nlon, jpband), dtype=np.float64)
    ssa_sw = np.zeros((nlev, nlat, nlon, jpband), dtype=np.float64)
    g_sw = np.zeros((nlev, nlat, nlon, jpband), dtype=np.float64)

    # Write all files
    def write(fname, data):
        data.astype(np.float64).tofile(os.path.join(outdir, fname))

    write('t_base.dat', t_base)
    write('t_warm.dat', t_warm)
    write('hus_base.dat', q_base)
    write('hus_warm.dat', q_warm)
    write('O3_base.dat', o3_base)
    write('O3_warm.dat', o3_warm)
    write('skt_base.dat', ts_base)
    write('skt_warm.dat', ts_warm)
    write('sp_base.dat', ps_base)
    write('sp_warm.dat', ps_warm)
    write('solarin_base.dat', solin_base)
    write('solarin_warm.dat', solin_warm)
    write('ssrd_base.dat', ssrd_base)
    write('ssrd_warm.dat', ssrd_warm)
    write('ssru_base.dat', ssru_base)
    write('ssru_warm.dat', ssru_warm)
    write('co2_b.dat', co2_base)
    write('co2_w.dat', co2_warm)
    write('cc_base.dat', cc_base)
    write('cc_warm.dat', cc_warm)
    write('clwc_base.dat', clwc_base)
    write('clwc_warm.dat', clwc_warm)
    write('ciwc_base.dat', ciwc_base)
    write('ciwc_warm.dat', ciwc_warm)
    write('aerosol_aod_lw_base.dat', aod_lw)
    write('aerosol_aod_lw_warm.dat', aod_lw)
    write('aerosol_aod_sw_base.dat', aod_sw)
    write('aerosol_aod_sw_warm.dat', aod_sw)
    write('aerosol_ssa_sw_base.dat', ssa_sw)
    write('aerosol_ssa_sw_warm.dat', ssa_sw)
    write('aerosol_g_sw_base.dat', g_sw)
    write('aerosol_g_sw_warm.dat', g_sw)

    print(f"Generated {len(os.listdir(outdir))} files in {outdir}")
    print(f"  T range: {t_base.min():.1f} - {t_base.max():.1f} K")
    print(f"  q range: {q_base.min():.2e} - {q_base.max():.2e} kg/kg")
    print(f"  O3 peak: {o3_base.max():.2e} kg/kg at {plev[np.argmax(o3_base)]:.0f} hPa")
    print(f"  Ts: {ts_base[0]:.1f} K, Ps: {ps_base[0]:.0f} Pa")
    print(f"  CO2: {co2_base[0]:.0f} / {co2_warm[0]:.0f} ppmv")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--outdir', default='data_prep',
                        help='Output directory for binary files')
    args = parser.parse_args()
    generate_test_data(args.outdir)
