"""Self-tests for cfram_core: sign conventions, shape, energy balance.

Run: python test_cfram_core.py
All tests should print PASS.
"""
import numpy as np
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cfram_core import (R_per_layer, S_per_layer, compute_drdt, solve_dT,
                        lump_atm_residual)


def test_R_layer_signs():
    """Test R_per_layer sign: with N+1 half-levels, output has N atm + 1 sfc."""
    # 3 half-levels → 2 atm layers + 1 surface row (3 entries in R)
    # Uniform F_net_up across the column ⇒ no atm cooling.
    lw_up = np.array([200., 250., 300.])
    lw_dn = np.array([0.,    50., 100.])
    R = R_per_layer(lw_up, lw_dn)
    assert R.shape == (3,), f'shape {R.shape}'
    assert np.isclose(R[0], 0.) and np.isclose(R[1], 0.), \
        f'atm cooling should be 0 with uniform F_net_up, got {R[:2]}'
    assert np.isclose(R[2], 200.), f'sfc emission - absorption = 300-100 = 200, got {R[2]}'
    print('PASS: R_per_layer shape (N+1) and uniform-flux → 0 atm cooling')

    # Top layer emits more than it receives from below (F_net_up decreases ↑):
    lw_up = np.array([300., 200., 200.])
    lw_dn = np.array([0.,    50., 100.])
    R = R_per_layer(lw_up, lw_dn)
    # F_net_up = [300, 150, 100]
    # R_atm[0] = 300 - 150 = 150  (top atm layer cooling)
    # R_atm[1] = 150 - 100 = 50   (lower atm layer cooling)
    # R_sfc = 200 - 100 = 100
    assert np.isclose(R[0], 150.) and np.isclose(R[1], 50.) and np.isclose(R[2], 100.)
    print('PASS: R_per_layer cooling sign (R > 0 when F_net_up decreases ↓)')


def test_S_layer_signs():
    sw_up = np.array([100., 50., 0.])
    sw_dn = np.array([400., 350., 300.])
    S = S_per_layer(sw_up, sw_dn)
    # F_net_dn = sw_dn - sw_up = [300, 300, 300]
    # S_atm = [0, 0]; S_sfc = 300 - 0 = 300
    assert np.isclose(S[0], 0.) and np.isclose(S[1], 0.)
    assert np.isclose(S[2], 300.)
    print('PASS: S_per_layer (uniform net dn → 0 atm; sfc absorbs 300)')


def test_drdt_diagonal_sign():
    """Heating a layer ⇒ more LW emitted ⇒ R increases ⇒ drdt[k,k] > 0."""
    # Mock 3-layer atm using simple Stefan-Boltzmann per layer (not realistic
    # but enough to test sign of finite-difference matrix).
    sigma = 5.67e-8
    nlev = 3

    def R_func(T_atm, Ts):
        # Each atm layer cools by 4*sigma*T^3 per K (Planck linear).
        # Surface emits sigma*Ts^4 - 0 (treat downward LW as 0).
        R = np.empty(nlev + 1)
        for k in range(nlev):
            R[k] = 0.5 * sigma * T_atm[k]**4   # half goes up (toy)
        R[nlev] = sigma * Ts**4
        return R

    T = np.array([220., 250., 280.])
    Ts = 290.
    drdt, R0 = compute_drdt(R_func, T, Ts, dT=0.5)
    # Diagonal must be positive
    diag = np.diag(drdt)
    assert (diag > 0).all(), f'drdt diagonal must be positive: {diag}'
    # Off-diagonal: in this toy, layers are independent, so off-diag ≈ 0
    off_diag = drdt - np.diag(diag)
    assert np.max(np.abs(off_diag)) < 1e-3
    print('PASS: compute_drdt diagonal > 0 (Planck radiates more when warmer)')


def test_solve_dT_signs():
    """For positive forcing (extra energy in), dT > 0 (warming)."""
    # Diagonal drdt_inv (uncoupled layers)
    drdt_inv = np.diag([0.3, 0.3, 0.3, 0.3])  # K per W/m²
    # Cooling perturbation (frc > 0 means X causes more layer cooling)
    frc = np.array([0., 1., 0., 0.])  # 1 W/m² extra cooling at layer 1
    dT = solve_dT(drdt_inv, frc)
    assert dT[1] < 0, 'extra cooling forcing → layer cools'
    # Warming perturbation via ΔS
    delta_S = np.array([0., 1., 0., 0.])
    dT2 = solve_dT(drdt_inv, np.zeros(4), delta_S)
    assert dT2[1] > 0, 'extra solar absorption → layer warms'
    print('PASS: solve_dT signs (cooling/warming forcings → cooling/warming dT)')


def test_solve_dT_co2_doubling():
    """Realistic-ish CO2 doubling sanity: greenhouse warms surface."""
    # Mock: doubling CO2 reduces R at surface (less LW escapes) by 4 W/m²
    # 30-layer column, drdt is approximately Stefan-Boltzmann at base.
    nlev = 30
    n = nlev + 1
    diag = np.full(n, 1.0)  # ≈ 1 W/m²/K Planck
    diag[-1] = 4.0          # surface stronger Planck (σ·4·T^3 ≈ 5 at 300K)
    drdt = np.diag(diag)
    drdt_inv = np.linalg.inv(drdt)

    # CO2 reduces R_sfc by 4 W/m² (less LW escapes from surface)
    # ΔR_sfc < 0 ⇒ frc[sfc] < 0
    frc = np.zeros(n); frc[-1] = -4.0
    dT = solve_dT(drdt_inv, frc)
    # Expect surface warming +4/4 = +1 K
    assert dT[-1] > 0, f'surface should warm, got {dT[-1]}'
    assert np.isclose(dT[-1], 1.0), f'expected ~1K, got {dT[-1]}'
    print('PASS: solve_dT CO2 surface warming (-4 W/m² ΔR ⇒ +1 K dT)')


def test_lump_atm_residual():
    frc = np.array([1., 2., 3., 5.])  # nlev=3 + sfc
    frc_atm, _ = lump_atm_residual(frc)
    assert np.allclose(frc_atm[:3], -frc[:3])
    assert frc_atm[3] == 0
    print('PASS: lump_atm_residual (atm rows negated, surface zeroed)')


if __name__ == '__main__':
    test_R_layer_signs()
    test_S_layer_signs()
    test_drdt_diagonal_sign()
    test_solve_dT_signs()
    test_solve_dT_co2_doubling()
    test_lump_atm_residual()
    print('\nAll cfram_core tests PASSED.')
