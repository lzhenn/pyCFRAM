"""Single-column CFRAM decomposition — pure math, no climlab dependency.

Conventions (verified against Lu & Cai 2009 Part II)
====================================================

Indexing
--------
Layer index k runs top-to-bottom: k=0 is TOA, k=nlev-1 is the lowest atm
layer, k=nlev is the surface "layer". Vectors of length nlev+1 cover atm + sfc.
Flux profiles (LW_up, LW_dn, SW_up, SW_dn) are at HALF levels, length nlev+1
(0 = TOA, nlev = surface), so each atm layer k sits between half-levels k and
k+1.

Per-layer R and S
-----------------
R_k = "net LW radiation flux LEAVING layer k" (W/m²), positive ⇒ layer cools.
  Atm:     R_k = LW_net_up(half k) - LW_net_up(half k+1)
            where LW_net_up = LW_up - LW_dn.
  Surface: R_sfc = LW_up(sfc) - LW_dn(sfc)  (= upward emission − absorbed dn).

S_k = "SW absorbed by layer k" (W/m²), positive ⇒ layer warms.
  Atm:     S_k = SW_net_dn(half k) - SW_net_dn(half k+1)
            where SW_net_dn = SW_dn - SW_up.
  Surface: S_sfc = SW_dn(sfc) - SW_up(sfc).

Energy balance (RCE steady state)
---------------------------------
R = S + Q_dyn   (Q_dyn = non-radiative heating: dynamics + latent + sensible)

Linearisation (Lu/Cai 2009 Eq.4 rewritten)
------------------------------------------
For a perturbation X (e.g. doubling CO2, replacing q1→q2), with T held at base:
  ΔR_X[k] := R(perturbed)[k] - R(base)[k]    (positive ⇒ X cools layer k)
  ΔS_X[k] := S(perturbed)[k] - S(base)[k]    (positive ⇒ X warms layer k)

In response, the column re-equilibrates with ΔT_X(p):
   drdt · ΔT_X = ΔS_X - ΔR_X     (= "extra energy received" by each layer)
   ⇒ ΔT_X = drdt⁻¹ · (ΔS_X - ΔR_X)
          = -drdt⁻¹ · (ΔR_X - ΔS_X)

`solve_dT` accepts (frc, delta_S=None) where frc=ΔR, delta_S=ΔS. For
LW-only perturbations (e.g. CO2), pass delta_S=None.

pyCFRAM Fortran's `solve_dT(frc) = -drdt_inv @ frc` matches this when
delta_S is implicitly zero (LW-only) and frc=ΔR.
"""

import numpy as np


def R_per_layer(lw_up, lw_dn):
    """Compute R per layer + surface from half-level LW flux profiles.

    Parameters
    ----------
    lw_up, lw_dn : (nlev+1,) arrays
        Upward and downward LW flux at half levels (TOA at index 0).

    Returns
    -------
    R : (nlev+1,) array
        Atm LW divergence (cooling) followed by surface row.
    """
    F_net_up = lw_up - lw_dn
    R_atm = F_net_up[:-1] - F_net_up[1:]
    R_sfc = lw_up[-1] - lw_dn[-1]
    return np.concatenate([R_atm, [R_sfc]])


def S_per_layer(sw_up, sw_dn):
    """Compute S per layer + surface from half-level SW flux profiles."""
    F_net_dn = sw_dn - sw_up
    S_atm = F_net_dn[:-1] - F_net_dn[1:]
    S_sfc = sw_dn[-1] - sw_up[-1]
    return np.concatenate([S_atm, [S_sfc]])


def compute_drdt(R_func, T_atm, Ts, dT=1.0):
    """Compute Planck matrix drdt[i, k] = ∂R_i/∂T_k via one-sided finite diff.

    Parameters
    ----------
    R_func : callable
        f(T_atm, Ts) -> R per layer (size nlev+1). Other state baked in.
    T_atm : (nlev,) array
        Base atmospheric T profile.
    Ts : float
        Base surface T.
    dT : float
        Perturbation magnitude (K). Default 1.0.

    Returns
    -------
    drdt : (nlev+1, nlev+1) array
        drdt[:, k] for k < nlev: response to atm layer k perturbation.
        drdt[:, nlev]: response to surface T perturbation.
    R_base : (nlev+1,) array
    """
    nlev = len(T_atm)
    R_base = R_func(T_atm, Ts)
    drdt = np.empty((nlev + 1, nlev + 1))
    for k in range(nlev):
        Tp = T_atm.copy(); Tp[k] += dT
        drdt[:, k] = (R_func(Tp, Ts) - R_base) / dT
    drdt[:, nlev] = (R_func(T_atm, Ts + dT) - R_base) / dT
    return drdt, R_base


def solve_dT(drdt_inv, frc, delta_S=None):
    """Linear CFRAM solve: dT_X = drdt⁻¹ · (ΔS_X - frc_X).

    Parameters
    ----------
    drdt_inv : (n+1, n+1) array
        Inverse of Planck matrix.
    frc : (n+1,) array
        ΔR per layer due to perturbation X (positive = X cools layer).
    delta_S : (n+1,) array or None
        ΔS per layer (positive = X warms layer via SW). None ≡ zero
        (LW-only perturbations like CO2).

    Returns
    -------
    dT : (n+1,) array
    """
    if delta_S is None:
        return -drdt_inv @ frc
    return drdt_inv @ (delta_S - frc)


def lump_atm_residual(frc, delta_S=None):
    """Build pyCFRAM-style "atmdyn" forcing: -frc[atm], surface row = 0.

    Used to compute dT_atmdyn = drdt⁻¹ · (ΔS_atm - (-ΔR_atm))
                              = drdt⁻¹ · (ΔS_atm + ΔR_atm)  [in atm rows]
    with the surface row contributing zero.
    """
    nlev_p1 = len(frc)
    frc_atm = np.zeros_like(frc)
    frc_atm[:-1] = -frc[:-1]
    if delta_S is None:
        return frc_atm, None
    dS_atm = np.zeros_like(delta_S)
    dS_atm[:-1] = -delta_S[:-1]
    return frc_atm, dS_atm
