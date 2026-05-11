"""Unified configuration loader for pyCFRAM.

Usage:
    from core.config import load_case, defaults, get_plev, get_aerosol_map, PROJECT_ROOT
"""
import os
import yaml
import numpy as np

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

_DEFAULTS = None


def defaults():
    """Load and cache global defaults from configs/defaults.yaml."""
    global _DEFAULTS
    if _DEFAULTS is None:
        path = os.path.join(PROJECT_ROOT, 'configs', 'defaults.yaml')
        with open(path) as f:
            _DEFAULTS = yaml.safe_load(f)
    return _DEFAULTS


def load_case(case_name):
    """Load case configuration from cases/<name>/case.yaml.

    Returns dict with resolved input paths and directory paths.
    """
    case_dir = os.path.join(PROJECT_ROOT, 'cases', case_name)
    cfg_path = os.path.join(case_dir, 'case.yaml')
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)

    # Inject directory paths
    cfg['_case_dir'] = case_dir
    cfg['_output_dir'] = os.path.join(case_dir, 'output')
    cfg['_figures_dir'] = os.path.join(case_dir, 'figures')

    # Resolve input file paths relative to case_dir
    input_keys = ['base_pres', 'base_surf', 'perturbed_pres', 'perturbed_surf',
                  'nonrad_forcing']
    for key in input_keys:
        if key in cfg.get('input', {}):
            val = cfg['input'][key]
            if not os.path.isabs(val):
                cfg['input'][key] = os.path.join(case_dir, val)

    return cfg


def get_nproc(case_cfg=None):
    """Get number of parallel workers. Default = all CPUs."""
    if case_cfg:
        n = case_cfg.get('run', {}).get('nproc', 'auto')
    else:
        n = defaults().get('run', {}).get('nproc', 'auto')
    if n == 'auto' or n is None:
        return os.cpu_count()
    return int(n)


def get_plev(case_cfg=None):
    """Get pressure levels (hPa, TOA→surface).

    Resolution order:
      1. case.yaml `grid.pressure_levels` — explicit override (e.g. 17/19-plev
         for CMIP6 sub-grids).
      2. input NetCDF `lev` variable — auto-derive from the actual data so a
         climlab 30-level RCE column or any other custom grid Just Works
         without yaml duplication.
      3. configs/defaults.yaml `grid.pressure_levels` — last-resort default.
    """
    if case_cfg and 'grid' in case_cfg and 'pressure_levels' in case_cfg['grid']:
        return np.array(case_cfg['grid']['pressure_levels'], dtype=np.float64)
    # Auto-derive from input NetCDF lev variable when available.
    if case_cfg:
        bp = case_cfg.get('input', {}).get('base_pres')
        if bp and os.path.exists(bp):
            from netCDF4 import Dataset    # local import: avoid hard dep at import time
            nc = Dataset(bp)
            try:
                plev = np.array(nc.variables['lev'][:], dtype=np.float64)
            finally:
                nc.close()
            # Convention: pyCFRAM internal arrays are TOA→surface. NetCDF stores
            # surface→TOA (matching CMIP6 plev). Detect and reverse if needed.
            if plev[0] > plev[-1]:
                plev = plev[::-1]
            return plev
    return np.array(defaults()['grid']['pressure_levels'], dtype=np.float64)


# Mapping from logical radiation scheme name → built Fortran binary in fortran/.
# Both engines run as single-column workers; nlev is inferred at runtime from the
# size of data_prep/plev.dat, so no per-grid recompilation is needed.
RADIATION_SCHEMES = {
    'fu':    'cfram_fu_1col',
    'rrtmg': 'cfram_rrtmg_1col',
    # extend here when additional schemes are added (e.g. 'cam', 'lw_only')
}


def get_executable(case_cfg=None):
    """Resolve Fortran binary name for a case.

    Resolution order:
      1. case.yaml `radiation.scheme` (preferred, schema-driven)
      2. case.yaml `run.executable`   (legacy escape hatch — direct binary name)
      3. default: cfram_rrtmg_1col

    Examples:
        radiation:
          scheme: fu                 # → cfram_fu_1col
        radiation:
          scheme: rrtmg              # → cfram_rrtmg_1col
        run:
          executable: my_custom_bin  # legacy: pass-through
    """
    if case_cfg:
        scheme = case_cfg.get('radiation', {}).get('scheme')
        if scheme:
            if scheme not in RADIATION_SCHEMES:
                raise ValueError(
                    "Unknown radiation.scheme '%s'. Known: %s"
                    % (scheme, sorted(RADIATION_SCHEMES))
                )
            return RADIATION_SCHEMES[scheme]
        exe = case_cfg.get('run', {}).get('executable')
        if exe:
            return exe
    return 'cfram_rrtmg_1col'


def get_drdt_eval(case_cfg=None):
    """Where to evaluate the Planck matrix `∂R/∂T`.

    Returns
    -------
    str : 'base' (default, 1st-order CFRAM) or 'midstate'
        'midstate' = (T_base + T_warm)/2 — equivalent to a 2-term Taylor
        expansion centered at the midpoint, reducing the linearisation
        residual from O(ΔT²) to O((ΔT/2)²) (~4× improvement). Currently
        supported by the RRTMG engine only; the Fu engine ignores the flag
        and always uses base.

    Per-case via case.yaml `radiation.drdt_eval`.
    """
    if case_cfg:
        v = case_cfg.get('radiation', {}).get('drdt_eval')
        if v:
            if v not in ('base', 'midstate'):
                raise ValueError(
                    "radiation.drdt_eval must be 'base' or 'midstate', got %r" % v)
            return v
    return 'base'


def get_drdt_probe(case_cfg=None):
    """Finite-difference scheme used to build the Planck Jacobian.

    Returns
    -------
    str : 'onesided' (default) or 'centered'
        'onesided' = standard: J_kj ≈ [R(T_j + 1K) - R(T_j)] / 1K
            Leading error: ½ R_TT|_T · 1K = O(R_TT) per column.
        'centered' = ±0.5K probe: J_kj ≈ [R(T_j + 0.5K) - R(T_j - 0.5K)] / 1K
            Cancels R_TT exactly, leading error (1/24)·R_TTT. Costs 2×
            rad_driver_lw calls in calc_drdt vs 1×, so the Planck-Jacobian
            stage is ~2× slower (small fraction of total cost). RRTMG-only.

    Per-case via case.yaml `radiation.drdt_probe`.
    """
    if case_cfg:
        v = case_cfg.get('radiation', {}).get('drdt_probe')
        if v:
            if v not in ('onesided', 'centered'):
                raise ValueError(
                    "radiation.drdt_probe must be 'onesided' or 'centered', got %r" % v)
            return v
    return 'onesided'


def get_q_handling(case_cfg=None):
    """How to compute the water-vapour partial perturbation `frc_q`.

    Returns
    -------
    str : 'independent' (default), 'feedback', or 'midstate'
        'independent' = standard CFRAM:
            frc_q = R(q_warm + T_base + others_base) - R(base)
            Treats q as an independent radiative variable. Correct for
            real atmospheres (ERA5, CESM2) where q has its own observed
            variability decoupled from local T.

        'feedback' = Manabe RH-fixed (one-sided warm path):
            frc_q = R(q_warm + T_warm + others_base)
                  - R(q_base + T_warm + others_base)
            Computes q radiative impact in the warm-T atmosphere.
            Avoids the supersaturation artifact at T_base. Costs +1
            rad_driver call/cell. Fixes surface closure but does NOT
            cancel the ∂²R/∂T∂q cross-term (cross-term sign is flipped).

        'midstate' = 2nd-order CFRAM q path:
            frc_q = R(q_warm + T_mid + q_mid_other + others_mid)
                  - R(q_base + T_mid + q_mid_other + others_mid)
            q perturbation evaluated in the midstate atmosphere
            (T_mid = (T_base+T_warm)/2, ...). Cancels the ∂²R/∂T∂q
            cross-term, analogous to co2_handling=midstate. Costs +2
            rad_driver calls/cell (q_warm and q_base, both in midstate
            atmosphere). RRTMG-only.

    RRTMG-only — Fu engine ignores the flag.
    """
    if case_cfg:
        v = case_cfg.get('radiation', {}).get('q_handling')
        if v:
            if v not in ('independent', 'feedback', 'midstate'):
                raise ValueError(
                    "radiation.q_handling must be 'independent', 'feedback', or 'midstate', got %r" % v)
            return v
    return 'independent'


def get_co2_handling(case_cfg=None):
    """How to compute the CO2 partial perturbation `frc_co2`.

    Returns
    -------
    str : 'base' (default) or 'midstate'
        'base' = standard CFRAM:
            frc_co2 = R(co2_warm + T_base + others_base) - R(base)
            CO2 perturbation evaluated in cold base atmosphere.

        'midstate' = 2nd-order CFRAM CO2 path:
            frc_co2 = R(co2_warm + T_mid + q_mid + others_mid)
                    - R(co2_base + T_mid + q_mid + others_mid)
            CO2 perturbation evaluated in the *midstate* atmosphere
            (T_mid = (T_base+T_warm)/2, q_mid = ..., etc.). This cancels
            the ∂²R/∂T∂C cross-term between CO2 and temperature, which is
            the dominant residual source in the 300-500 hPa CO2 15-μm
            band saturation transition. Costs +2 rad_driver calls per
            cell (warm & base CO2 in midstate atmosphere).

    RRTMG-only — Fu engine ignores the flag.
    """
    if case_cfg:
        v = case_cfg.get('radiation', {}).get('co2_handling')
        if v:
            if v not in ('base', 'midstate'):
                raise ValueError(
                    "radiation.co2_handling must be 'base' or 'midstate', got %r" % v)
            return v
    return 'base'


def get_output_terms(case_cfg=None):
    """List of dT_*/frc_* terms to write to cfram_result.nc, or None for all.

    Per-case via case.yaml `radiation.output_terms`. Useful for idealized
    experiments (e.g. clear-sky climlab RCE) where only a subset of partial
    perturbations are physically meaningful — listing only `[co2, q]` strips
    the always-zero `cloud / aerosol / ts / albedo / o3 / solar` rows from
    the output NetCDF.

    Note: derived dynamics terms (`dry/atmdyn/sfcdyn/ocndyn/observed/full/warm`
    plus `lhflx/shflx`) are *not* filtered — they are physically required for
    energy-balance closure and always written.

    Returns
    -------
    list[str] or None
        None means "write everything" (default behavior).
    """
    if case_cfg:
        terms = case_cfg.get('radiation', {}).get('output_terms')
        if terms:
            return list(terms)
    return None


def get_aerosol_map():
    """Get aerosol species → lookup table mapping from config."""
    return defaults()['aerosol']['species']


def get_fortran_dir():
    """Get path to Fortran executable directory."""
    return os.path.join(PROJECT_ROOT, 'fortran')


def get_lookup_dir():
    """Get path to aerosol optical property lookup tables."""
    return os.path.join(PROJECT_ROOT, 'fortran', 'data_prep', 'aerosol')
