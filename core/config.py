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


def get_plev():
    """Get pressure levels array (hPa, TOA→surface)."""
    return np.array(defaults()['grid']['pressure_levels'], dtype=np.float64)


def get_aerosol_map():
    """Get aerosol species → lookup table mapping from config."""
    return defaults()['aerosol']['species']


def get_fortran_dir():
    """Get path to Fortran executable directory."""
    return os.path.join(PROJECT_ROOT, 'fortran')


def get_lookup_dir():
    """Get path to aerosol optical property lookup tables."""
    return os.path.join(PROJECT_ROOT, 'fortran', 'data_prep', 'aerosol')
