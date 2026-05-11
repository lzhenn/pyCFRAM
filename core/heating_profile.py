"""Vertical distribution of column-integrated non-radiative fluxes.

Use case: in CFRAM, the non-radiative forcing (latent / sensible heat) needs to
appear as a vertical *profile* of layer heating rates (W/m^2 per layer), not as
a single surface flux. When a tendency-resolving GCM is used as input, the
column profile can come directly from model diagnostics (DTCOND, ZMDT, DTV).
When only a surface bulk flux is available (reanalysis, idealized RCM), we fall
back to a parameterized vertical shape.

Two parameterizations are provided:

* ``lu_cai_profile`` - condensation heating profile from Cai & Lu (2009)
  Part II Eq.(1): quadratic shape, peak at p_mid=560 hPa, zero outside
  [p_min=292, p_max=831] hPa. Vertically integrates to surface evaporation
  rate. This is what the original RCM uses for Q^lh.
* ``pbl_exp_profile`` - sensible heat / vertical diffusion profile,
  exponential decay above surface with PBL scale height (default
  H = 100 hPa). Vertically integrates to surface SH flux.

Both profiles are layer-integrated, normalized so summing across levels
at a grid point equals 1 (so multiplying by surface flux conserves the
column-integrated flux).

Public API
----------
``distribute_to_column(surface_flux_2d, plev, ps, profile='lu_cai', **kwargs)``
    Take a 2D surface flux field and return a 3D column-distributed forcing
    that integrates back to the surface flux at every grid point.

Tendency-direct path (preferred when available) lives in
``data/cesm_history.py`` (or any GCM-specific reader). It produces the same
shape 3D field directly from heating tendencies (DTCOND/ZMDT/DTV), bypassing
this fallback.
"""

import numpy as np


# Cai & Lu (2009) Part II Eq.(1) bounds (hPa)
LU_CAI_PMIN = 292.0
LU_CAI_PMAX = 831.0


def _layer_thickness_hpa(plev):
    """Layer thicknesses (hPa) from level midpoints; preserves input order."""
    plev = np.asarray(plev, dtype=np.float64)
    order = np.argsort(plev)
    p_sorted = plev[order]
    p_half = np.zeros(len(p_sorted) + 1)
    p_half[1:-1] = 0.5 * (p_sorted[:-1] + p_sorted[1:])
    p_half[0] = max(0.0, p_sorted[0] - 0.5 * (p_sorted[1] - p_sorted[0]))
    p_half[-1] = p_sorted[-1] + 0.5 * (p_sorted[-1] - p_sorted[-2])
    dp_sorted = np.diff(p_half)
    dp = np.empty_like(plev)
    dp[order] = dp_sorted
    return dp


def lu_cai_profile(plev, ps_hpa=None):
    """Normalized condensation heating profile, Cai & Lu 2009 Eq.(1).

    Quadratic shape, zero outside [pmin, pmax], peaking at
    pmid = (pmin + pmax) / 2 ≈ 560 hPa::

        f(p) ∝ 1 - ((p - pmid) / dp)^2,    pmin <= p <= pmax

    Returns layer-integrated weights summing to 1.

    Parameters
    ----------
    plev : (nlev,) array
        Pressure levels in hPa, any order. Output preserves input order.
    ps_hpa : ignored
        API-symmetric placeholder (Lu/Cai profile is terrain-independent).

    Returns
    -------
    weights : (nlev,) array, sum == 1.
    """
    plev = np.asarray(plev, dtype=np.float64)
    pmid = 0.5 * (LU_CAI_PMIN + LU_CAI_PMAX)
    dp_band = 0.5 * (LU_CAI_PMAX - LU_CAI_PMIN)

    inside = (plev >= LU_CAI_PMIN) & (plev <= LU_CAI_PMAX)
    shape = np.zeros_like(plev)
    shape[inside] = 1.0 - ((plev[inside] - pmid) / dp_band) ** 2

    dp = _layer_thickness_hpa(plev)
    w = shape * dp
    s = w.sum()
    return w / s if s > 0 else w


def pbl_exp_profile(plev, ps_hpa, scale_hpa=100.0):
    """Exponential PBL profile: peaks at surface, decays upward.

    weight(p) ∝ exp(-(ps - p) / scale)  for p <= ps,  else 0.

    Layer-integrated weights summing to 1 (per column when ps is 2D).

    Parameters
    ----------
    plev : (nlev,) array
        Pressure levels in hPa.
    ps_hpa : scalar or (nlat, nlon) array
        Surface pressure in hPa.
    scale_hpa : float
        Decay scale (hPa). Default 100 hPa ~ 1 km PBL height.

    Returns
    -------
    weights : (nlev,) if ps_hpa scalar; (nlev, nlat, nlon) if 2D.
    """
    plev = np.asarray(plev, dtype=np.float64)
    ps_hpa = np.asarray(ps_hpa, dtype=np.float64)
    dp = _layer_thickness_hpa(plev)

    if ps_hpa.ndim == 0:
        shape = np.where(plev <= ps_hpa,
                         np.exp(-(ps_hpa - plev) / scale_hpa),
                         0.0)
        w = shape * dp
        s = w.sum()
        return w / s if s > 0 else w

    p3 = plev[:, None, None]
    ps3 = ps_hpa[None, :, :]
    shape = np.where(p3 <= ps3, np.exp(-(ps3 - p3) / scale_hpa), 0.0)
    w = shape * dp[:, None, None]
    s = w.sum(axis=0, keepdims=True)
    return np.where(s > 0, w / s, 0.0)


PROFILES = {
    'lu_cai': lu_cai_profile,
    'pbl_exp': pbl_exp_profile,
}


def expand_surface_to_column_energy_conserved(
        surface_flux_2d, plev_hpa, ps_hpa=None,
        profile='lu_cai', surface_level_index=None, **kwargs):
    """Build energy-conserving column forcing from a surface flux.

    Physical model (Lu & Cai 2009 Eq.(1)+(2)): a surface latent flux removes
    energy from the surface and deposits it back into the atmosphere column
    via condensation, so the COLUMN-INTEGRAL of the resulting Q^lh forcing
    vector is **zero**. Surface row carries the original sign; atmosphere
    rows carry the opposite sign with vertical distribution given by
    ``profile``.

    For sensible heat (pbl_exp), same logic: surface flux into PBL becomes
    column-integrated heating in the lower atmosphere, with surface row
    holding the negative.

    Parameters
    ----------
    surface_flux_2d : (nlat, nlon) array
        Surface flux in W/m^2 (sign convention as provided by user; output
        preserves it).
    plev_hpa : (nlev,) array
        Pressure levels.
    ps_hpa : (nlat, nlon) array or None
        Surface pressure for terrain-aware profiles.
    profile : str or callable
    surface_level_index : int or None
        Index in ``plev_hpa`` where the surface row sits. If None, defaults
        to the highest-pressure level (argmax of plev).

    Returns
    -------
    column_forcing : (nlev, nlat, nlon) array
        Energy-conserving column distribution. Sums to ~zero over each
        grid point's column (within roundoff).
    """
    surface_flux_2d = np.asarray(surface_flux_2d, dtype=np.float64)
    plev_hpa = np.asarray(plev_hpa, dtype=np.float64)

    if surface_level_index is None:
        surface_level_index = int(np.argmax(plev_hpa))

    atm_part = -distribute_to_column(surface_flux_2d, plev_hpa,
                                     ps_hpa=ps_hpa, profile=profile,
                                     **kwargs)
    out = atm_part.copy()
    out[surface_level_index, :, :] += surface_flux_2d
    return out


def distribute_to_column(surface_flux_2d, plev_hpa, ps_hpa=None,
                         profile='lu_cai', **kwargs):
    """Distribute a 2D surface flux field over the vertical column.

    Output sign matches input: an upward latent flux of -80 W/m^2 distributes
    to a column whose layer-summed flux equals -80 W/m^2 at every column.

    Parameters
    ----------
    surface_flux_2d : (nlat, nlon) array
        Surface flux in W/m^2.
    plev_hpa : (nlev,) array
        Pressure levels in hPa.
    ps_hpa : (nlat, nlon) array or None
        Surface pressure. Required for terrain-aware profiles ('pbl_exp');
        ignored by 'lu_cai'.
    profile : str or callable
        One of {'lu_cai', 'pbl_exp'} or a callable
        ``f(plev, ps, **kwargs) -> weights`` returning normalized weights.
    **kwargs
        Forwarded to the profile callable (e.g. ``scale_hpa`` for pbl_exp).

    Returns
    -------
    column_flux : (nlev, nlat, nlon) array
        Layer-integrated forcing (W/m^2 per layer). At each grid point,
        ``sum_k column_flux[k, i, j] == surface_flux_2d[i, j]`` within
        roundoff.
    """
    surface_flux_2d = np.asarray(surface_flux_2d, dtype=np.float64)
    plev_hpa = np.asarray(plev_hpa, dtype=np.float64)

    if isinstance(profile, str):
        if profile not in PROFILES:
            raise ValueError("unknown profile %r, choose from %s" %
                             (profile, list(PROFILES)))
        profile_name = profile
        profile_fn = PROFILES[profile]
    else:
        profile_name = getattr(profile, '__name__', 'custom')
        profile_fn = profile

    if profile_name == 'lu_cai_profile' or profile_name == 'lu_cai':
        weights = profile_fn(plev_hpa, **kwargs)         # (nlev,)
        column = weights[:, None, None] * surface_flux_2d[None, :, :]
    elif profile_name == 'pbl_exp_profile' or profile_name == 'pbl_exp':
        if ps_hpa is None:
            raise ValueError("pbl_exp requires ps_hpa")
        weights = profile_fn(plev_hpa, ps_hpa, **kwargs)  # (nlev, nlat, nlon)
        column = weights * surface_flux_2d[None, :, :]
    else:
        # Custom callable: try with both signatures
        try:
            weights = profile_fn(plev_hpa, ps_hpa, **kwargs)
        except TypeError:
            weights = profile_fn(plev_hpa, **kwargs)
        if weights.ndim == 1:
            column = weights[:, None, None] * surface_flux_2d[None, :, :]
        else:
            column = weights * surface_flux_2d[None, :, :]

    return column
