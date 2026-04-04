"""Planck feedback matrix construction and CFRAM solver.

Translates drdt.f90 and the matrix solve logic from cfram_rrtmg.f90.

The Planck matrix ∂R/∂T is constructed by perturbing each atmospheric layer
and the surface temperature by +1K, running LW radiation, and recording
the change in LW flux convergence.

The CFRAM equation:
    ΔT_i = -(∂R/∂T)⁻¹ × ΔR_i
"""

import numpy as np
from scipy import linalg


def build_planck_matrix(rad_lw_func, base_inputs, lw_base, fdl_base, ful_base,
                        nlayer):
    """Construct the Planck feedback matrix ∂R/∂T.

    Perturbs each atmospheric layer by +1K and the surface by +1K,
    runs LW-only radiation, and computes the Jacobian.

    Args:
        rad_lw_func: callable(inputs_dict) -> (lw, fdl, ful)
            Function that runs LW radiation and returns:
            - lw: (nlayer,) LW flux convergence per layer (ERA order)
            - fdl: (nlayer+1,) LW downward flux (ERA order)
            - ful: (nlayer+1,) LW upward flux (ERA order)
        base_inputs: dict of RRTMG input arrays for base state
        lw_base: (nlayer,) baseline LW flux convergence
        fdl_base: (nlayer+1,) baseline LW downward flux
        ful_base: (nlayer+1,) baseline LW upward flux
        nlayer: number of active atmospheric layers

    Returns:
        drdt: (nlayer+1, nlayer+1) Planck matrix
              Rows: which layer's radiation changes
              Columns: which layer's temperature was perturbed
              Last row/col = surface
    """
    drdt = np.zeros((nlayer + 1, nlayer + 1))

    # Surface net LW flux at baseline
    sfc_net_base = fdl_base[nlayer] - ful_base[nlayer]

    # Perturb each atmospheric layer by +1K
    for k in range(nlayer):
        perturbed = _perturb_layer_temperature(base_inputs, k, nlayer, dT=1.0)
        lw_1k, fdl_1k, ful_1k = rad_lw_func(perturbed)

        # Row k: how radiation changes when layer k is warmed by 1K
        drdt[k, :nlayer] = lw_1k[:nlayer] - lw_base[:nlayer]
        drdt[k, nlayer] = (fdl_1k[nlayer] - ful_1k[nlayer]) - sfc_net_base

    # Perturb surface temperature by +1K
    perturbed = _perturb_surface_temperature(base_inputs, dT=1.0)
    lw_1k, fdl_1k, ful_1k = rad_lw_func(perturbed)

    drdt[nlayer, :nlayer] = lw_1k[:nlayer] - lw_base[:nlayer]
    drdt[nlayer, nlayer] = (fdl_1k[nlayer] - ful_1k[nlayer]) - sfc_net_base

    return drdt


def _perturb_layer_temperature(inputs, layer_idx, nlayer, dT=1.0):
    """Create a copy of inputs with one atmospheric layer perturbed by dT.

    Note: layer_idx is in ERA order (0 = TOA). The RRTMG-order arrays
    in inputs need the corresponding index flipped.

    In RRTMG order: rrtmg_idx = nlayer - 1 - layer_idx
    """
    perturbed = {}
    for key, val in inputs.items():
        if isinstance(val, np.ndarray):
            perturbed[key] = val.copy()
        else:
            perturbed[key] = val

    # Perturb tave in RRTMG order
    rrtmg_idx = nlayer - 1 - layer_idx
    perturbed['tave'][rrtmg_idx] += dT

    # Also perturb interface temperatures
    # In RRTMG order: tz[0]=surface, tz[nlayer]=TOA
    # Interface between rrtmg layers rrtmg_idx and rrtmg_idx+1
    if rrtmg_idx < nlayer:
        perturbed['tz'][rrtmg_idx] += dT * 0.5
    if rrtmg_idx + 1 <= nlayer:
        perturbed['tz'][rrtmg_idx + 1] += dT * 0.5

    return perturbed


def _perturb_surface_temperature(inputs, dT=1.0):
    """Create a copy of inputs with surface temperature perturbed by dT."""
    perturbed = {}
    for key, val in inputs.items():
        if isinstance(val, np.ndarray):
            perturbed[key] = val.copy()
        else:
            perturbed[key] = val

    perturbed['tbound'] = inputs['tbound'] + dT
    # tz[0] is surface in RRTMG order
    perturbed['tz'][0] += dT

    return perturbed


def solve_cfram(drdt, forcing, nlayer):
    """Solve CFRAM equation for partial temperature changes.

    ΔT = -(∂R/∂T)⁻¹_atm × ΔR

    Uses only the atmospheric part of the Planck matrix (excluding surface),
    as done in the Fortran code.

    Args:
        drdt: (nlayer+1, nlayer+1) full Planck matrix
        forcing: (nlayer+1,) or (nlayer,) radiative forcing vector
              If nlayer+1: includes surface forcing (only atm part used)
              If nlayer: atmospheric forcing only

    Returns:
        dT: (nlayer,) partial temperature change per layer
    """
    # Extract atmospheric part
    drdt_atm = drdt[:nlayer, :nlayer]

    # Use only atmospheric forcing
    frc_atm = forcing[:nlayer]

    # Solve: drdt_atm × dT = -frc_atm
    # Equivalent to: dT = drdt_atm_inv × (-frc_atm)
    try:
        dT = linalg.solve(drdt_atm, -frc_atm)
    except linalg.LinAlgError:
        # Fallback to pseudo-inverse for singular matrices
        drdt_atm_inv = linalg.pinv(drdt_atm)
        dT = drdt_atm_inv @ (-frc_atm)

    return dT


def solve_cfram_multi(drdt, forcings_dict, nlayer):
    """Solve CFRAM for multiple forcing terms efficiently.

    Inverts the Planck matrix once and applies to all forcing vectors.

    Args:
        drdt: (nlayer+1, nlayer+1) full Planck matrix
        forcings_dict: dict of {term_name: (nlayer+1,) forcing vector}
        nlayer: number of layers

    Returns:
        dT_dict: dict of {term_name: (nlayer,) temperature change}
    """
    drdt_atm = drdt[:nlayer, :nlayer]

    try:
        drdt_atm_inv = linalg.inv(drdt_atm)
    except linalg.LinAlgError:
        drdt_atm_inv = linalg.pinv(drdt_atm)

    dT_dict = {}
    for term, frc in forcings_dict.items():
        frc_atm = frc[:nlayer]
        dT_dict[term] = drdt_atm_inv @ (-frc_atm)

    return dT_dict
