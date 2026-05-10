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
