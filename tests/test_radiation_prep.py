#!/usr/bin/env python3
"""Test radiation.py preprocessing against known values.

Validates coordinate transformations, column densities, cloud water paths
using the synthetic standard atmosphere from generate_smoke_test_data.py.
"""

import sys, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from core.constants import ERA5_PLEVELS, NLEV
from core.radiation import (compute_interface_levels, compute_column_densities,
                            compute_cloud_water_path, flip_to_rrtmg,
                            flip_interfaces_to_rrtmg, compute_flux_convergence,
                            prepare_rrtmg_inputs)
from tests.generate_smoke_test_data import (us_standard_temperature,
                                             us_standard_humidity)


def test_interface_levels():
    """Test interface pressure/temperature computation."""
    plev = ERA5_PLEVELS  # 1, 2, 3, ..., 1000 hPa
    t = us_standard_temperature(plev)
    ps = 1013.25  # hPa
    ts = 288.15   # K

    pint, tint, dp = compute_interface_levels(plev, t, ps, ts)

    # Check dimensions
    assert len(pint) == NLEV + 1, f"pint length {len(pint)} != {NLEV+1}"
    assert len(dp) == NLEV, f"dp length {len(dp)} != {NLEV}"

    # Check boundaries
    assert pint[0] == 0.5, f"TOA pressure {pint[0]} != 0.5"
    assert pint[NLEV] == ps, f"Surface pressure {pint[NLEV]} != {ps}"
    assert tint[0] == t[0], f"TOA temperature {tint[0]} != {t[0]}"
    assert tint[NLEV] == ts, f"Surface temperature {tint[NLEV]} != {ts}"

    # dp should be positive (pressure increases downward)
    assert np.all(dp > 0), "Some dp values are negative"

    # Sum of dp should equal ps - 0.5
    assert abs(np.sum(dp) - (ps - 0.5)) < 1e-10, \
        f"Sum(dp)={np.sum(dp)} != {ps - 0.5}"

    print("  interface_levels: PASS")


def test_column_densities():
    """Test column density computation."""
    plev = ERA5_PLEVELS
    t = us_standard_temperature(plev)
    q = us_standard_humidity(plev, t)
    o3 = np.full(NLEV, 1e-7)  # simple ozone
    ps = 1013.25
    ts = 288.15

    pint, tint, dp = compute_interface_levels(plev, t, ps, ts)
    wkl, coldry, wbrodl = compute_column_densities(
        plev, t, q, o3, 380.0, 1.6, 0.28, pint, dp)

    # Column densities should be positive
    assert np.all(coldry > 0), "Some coldry values are negative"
    assert np.all(wbrodl > 0), "Some wbrodl values are negative"

    # wkl[0] = H2O mixing ratio, should be < 0.05 (5% volume)
    assert np.all(wkl[0] < 0.05), f"H2O mixing ratio too high: {wkl[0].max()}"

    # wkl[6] = O2, should be 0.209
    assert np.allclose(wkl[6], 0.209), f"O2 mixing ratio wrong: {wkl[6, 0]}"

    # Total column: order of magnitude ~10^25 molecules/cm² for whole atmosphere
    total_col = np.sum(coldry)
    assert 1e24 < total_col < 1e26, f"Total coldry {total_col} out of range"

    print("  column_densities: PASS")


def test_cloud_water_path():
    """Test cloud water path computation."""
    nlayer = 5
    cldfrac = np.array([0.0, 0.5, 0.8, 0.3, 0.0])
    cldiwc = np.array([0.0, 1e-5, 2e-5, 1e-5, 0.0])  # kg/kg
    cldlwc = np.array([0.0, 5e-5, 3e-5, 2e-5, 0.0])
    dp = np.array([10.0, 50.0, 100.0, 50.0, 25.0])  # hPa

    ciwp, clwp = compute_cloud_water_path(cldfrac, cldiwc, cldlwc, dp)

    # Clear layers should have zero
    assert ciwp[0] == 0.0 and ciwp[4] == 0.0
    assert clwp[0] == 0.0 and clwp[4] == 0.0

    # Cloudy layers should have positive values
    assert ciwp[1] > 0 and clwp[1] > 0

    # Check formula: coef = (1/9.8) * 100 * 1000 = 10204.08
    coef = (1.0 / 9.8) * 1e2 * 1e3
    expected_ciwp_1 = coef * 1e-5 * 50.0
    assert abs(ciwp[1] - expected_ciwp_1) < 1e-6, \
        f"ciwp[1]={ciwp[1]} != {expected_ciwp_1}"

    print("  cloud_water_path: PASS")


def test_coordinate_flip():
    """Test ERA-to-RRTMG coordinate transformation."""
    nlayer = 5
    era = np.array([1.0, 2.0, 3.0, 4.0, 5.0])  # top to bottom
    rrtmg = flip_to_rrtmg(era, nlayer)

    # RRTMG should be bottom to top
    expected = np.array([5.0, 4.0, 3.0, 2.0, 1.0])
    assert np.allclose(rrtmg, expected), f"Flip failed: {rrtmg} != {expected}"

    # Interface flip
    pint_era = np.array([0.5, 1.5, 2.5, 3.5, 4.5, 1013.0])  # 0..nlayer
    tint_era = np.array([200, 210, 220, 250, 270, 288])
    pz, tz = flip_interfaces_to_rrtmg(pint_era, tint_era, nlayer)

    assert pz[0] == 1013.0, f"pz[0]={pz[0]} should be surface"
    assert pz[nlayer] == 0.5, f"pz[{nlayer}]={pz[nlayer]} should be TOA"

    print("  coordinate_flip: PASS")


def test_flux_convergence():
    """Test flux convergence computation."""
    nlayer = 3
    # fd[0]=TOA, fd[3]=surface
    fd = np.array([0.0, 50.0, 100.0, 200.0])
    fu = np.array([250.0, 200.0, 150.0, 50.0])

    conv = compute_flux_convergence(fd, fu, nlayer)

    # conv[0] = fd[0] + fu[1] - fd[1] - fu[0] = 0+200-50-250 = -100
    assert conv[0] == -100.0, f"conv[0]={conv[0]}"
    # conv[1] = fd[1] + fu[2] - fd[2] - fu[1] = 50+150-100-200 = -100
    assert conv[1] == -100.0, f"conv[1]={conv[1]}"

    print("  flux_convergence: PASS")


def test_prepare_rrtmg_inputs():
    """Test full input preparation pipeline."""
    plev = ERA5_PLEVELS
    t = us_standard_temperature(plev)
    q = us_standard_humidity(plev, t)
    o3 = np.full(NLEV, 1e-7)
    nlayer = 37  # all levels active for standard atmosphere

    inputs = prepare_rrtmg_inputs(
        nlayer=nlayer, plev=plev, t=t, q=q, o3=o3,
        co2_ppmv=380.0, ch4_ppmv=1.6, n2o_ppmv=0.28,
        ps=1013.25, ts=288.15, zenith=0.5, albedo_sw=0.2, albedo_lw=0.0,
        cldfrac=np.zeros(NLEV), cldlwc=np.zeros(NLEV), cldiwc=np.zeros(NLEV),
        tauaer_sw=None, ssaaer_sw=None, asmaer_sw=None, tauaer_lw=None,
        icld=0, iaer=0,
    )

    # Check RRTMG order: pz[0] should be surface (highest pressure)
    assert inputs['pz'][0] > inputs['pz'][nlayer], \
        f"pz not in RRTMG order: pz[0]={inputs['pz'][0]}, pz[N]={inputs['pz'][nlayer]}"

    # tave[0] should be near surface temperature
    assert inputs['tave'][0] > 270, \
        f"tave[0]={inputs['tave'][0]} too cold for surface"

    # coldry should be positive everywhere
    assert np.all(inputs['coldry'][:nlayer] > 0), "Negative coldry"

    # wkl should be column amounts (large numbers)
    assert inputs['wkl'][0, 0] > 1e15, \
        f"wkl[H2O,0]={inputs['wkl'][0,0]} too small for column amount"

    print("  prepare_rrtmg_inputs: PASS")


if __name__ == '__main__':
    print("=== Testing radiation preprocessing ===")
    test_interface_levels()
    test_column_densities()
    test_cloud_water_path()
    test_coordinate_flip()
    test_flux_convergence()
    test_prepare_rrtmg_inputs()
    print("\n=== ALL TESTS PASSED ===")
