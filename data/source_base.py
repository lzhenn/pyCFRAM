"""Abstract data source and factory registry for pyCFRAM input generation.

Each data source (ERA5, MERRA-2, CMIP6, ...) implements DataSource.build_states()
to produce base and perturbed atmospheric states in the standard pyCFRAM format
defined by docs/input_spec.md.

Usage:
    from data.source_base import get_source
    source = get_source(case_cfg)
    base_state, perturbed_state = source.build_states()
"""

import numpy as np


class DataSource:
    """Abstract base class for reanalysis/model data sources.

    Subclasses must implement build_states() which returns two dicts
    (base_state, perturbed_state). Each dict contains:

        3D variables: np.array(nlev, nlat, nlon)  -- surface→TOA level order
        2D variables: np.array(nlat, nlon)
        Coordinates:  lat (1D), lon (1D), lev (1D, hPa, surface→TOA)

    Variable names must match docs/input_spec.md exactly:
        3D: ta_lay, q, o3, camt, cliq, cice, co2, bc, ocphi, ocpho, sulf, ss, dust
        2D: ts, ps, solar, albedo
    """

    def __init__(self, case_cfg):
        self.cfg = case_cfg
        self.source_cfg = case_cfg.get('source', {})

    def build_states(self):
        """Compute base and perturbed atmospheric states.

        Returns:
            (base_state, perturbed_state, nonrad_forcing) where:
            - base_state, perturbed_state: dicts of arrays + coordinates
            - nonrad_forcing: dict with 'lhflx', 'shflx' (nlat, nlon) W/m²
              or empty dict if not available
        """
        raise NotImplementedError

    def get_co2(self, period='base'):
        """Get CO2 volume mixing ratio (mol/mol) for a period.

        Args:
            period: 'base' or 'perturbed'
        """
        co2_cfg = self.source_cfg.get('co2', {})
        src = co2_cfg.get('source', 'constant')
        if src == 'constant':
            if period == 'base':
                ppmv = co2_cfg.get('base_ppmv', 415.0)
            else:
                ppmv = co2_cfg.get('perturbed_ppmv', 418.0)
            return ppmv * 1e-6  # ppmv → mol/mol
        raise ValueError(f"Unknown CO2 source: {src}")

    def get_aerosol(self, shape, target_plev_hpa=None,
                    target_lat=None, target_lon=None,
                    data_dir=None, years=None, month=None,
                    warm_days=None, period='clim'):
        """Get aerosol mixing ratios for all 6 GOCART species.

        Args:
            shape: (nlev, nlat, nlon) for 3D arrays (used by zero-fill)
            target_plev_hpa: pressure levels in hPa (for merra2 interpolation)
            target_lat, target_lon: target grid coordinates
            data_dir: MERRA-2 data directory
            years, month, warm_days: temporal parameters for MERRA-2
            period: 'clim' or year number

        Returns:
            dict of {species_name: (nlev, nlat, nlon)} arrays
        """
        aer_cfg = self.source_cfg.get('aerosol', {})
        src = aer_cfg.get('source', 'zero')
        if src == 'zero':
            species = ['bc', 'ocphi', 'ocpho', 'sulf', 'ss', 'dust']
            return {s: np.zeros(shape, dtype=np.float64) for s in species}
        if src == 'merra2':
            from .merra2_aerosol import load_merra2_aerosol
            import os
            m2_dir = aer_cfg.get('data_dir', 'era5_data/merra2')
            if not os.path.isabs(m2_dir):
                from core.config import PROJECT_ROOT
                m2_dir = os.path.join(PROJECT_ROOT, m2_dir)
            return load_merra2_aerosol(
                data_dir=m2_dir, years=years, month=month,
                warm_days=warm_days,
                target_plev_hpa=target_plev_hpa,
                target_lat=target_lat, target_lon=target_lon,
                period=period)
        raise ValueError(f"Unknown aerosol source: {src}")


# ── Registry ─────────────────────────────────────────────────────────────────

_SOURCES = {}


def register_source(name):
    """Decorator to register a DataSource subclass under a name."""
    def decorator(cls):
        _SOURCES[name] = cls
        return cls
    return decorator


def get_source(case_cfg):
    """Factory: instantiate the right DataSource from case config.

    The source type is read from case_cfg['source']['type'].
    """
    src_type = case_cfg.get('source', {}).get('type')
    if src_type is None:
        raise ValueError("No 'source.type' in case config. "
                         "Add a 'source:' section to case.yaml.")
    if src_type not in _SOURCES:
        available = list(_SOURCES.keys())
        raise ValueError(f"Unknown source type '{src_type}'. "
                         f"Available: {available}")
    return _SOURCES[src_type](case_cfg)
