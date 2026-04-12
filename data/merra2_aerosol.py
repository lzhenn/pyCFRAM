"""MERRA-2 aerosol mixing ratio loader for pyCFRAM.

Loads 3D aerosol mixing ratios from MERRA-2 M2I3NVAER (inst3_3d_aer_Nv),
performs vertical interpolation from 72 hybrid sigma levels to 37 ERA5
pressure levels, and horizontal interpolation from 0.5°×0.625° to the
target ERA5 grid (0.25°×0.25°).

Species mapping:
    MERRA-2 variable    → pyCFRAM variable
    BCPHILIC + BCPHOBIC → bc
    OCPHILIC            → ocphi
    OCPHOBIC            → ocpho
    SO4                 → sulf
    SS001+...+SS005     → ss
    DU001+...+DU005     → dust

All output units are kg/kg (mass mixing ratio), matching ERA5 convention.
"""

import os
import glob
import numpy as np
from netCDF4 import Dataset
from scipy.interpolate import RegularGridInterpolator


# MERRA-2 fixed top-of-atmosphere pressure (Pa)
PTOP = 1.0  # 0.01 hPa = 1 Pa


def _load_merra2_day(filepath):
    """Load one day of MERRA-2 aerosol data.

    Returns:
        aer_vars: dict of {varname: (ntime, nlev, nlat, nlon)} arrays
        delp: (ntime, nlev, nlat, nlon) layer pressure thickness (Pa)
        lat, lon, lev: coordinate arrays
    """
    nc = Dataset(filepath, 'r')

    aer_names = ['BCPHILIC', 'BCPHOBIC', 'OCPHILIC', 'OCPHOBIC', 'SO4',
                 'SS001', 'SS002', 'SS003', 'SS004', 'SS005',
                 'DU001', 'DU002', 'DU003', 'DU004', 'DU005']

    aer_vars = {}
    for v in aer_names:
        if v in nc.variables:
            aer_vars[v] = np.array(nc.variables[v][:], dtype=np.float64)

    delp = np.array(nc.variables['DELP'][:], dtype=np.float64)
    lat = np.array(nc.variables['lat'][:], dtype=np.float64)
    lon = np.array(nc.variables['lon'][:], dtype=np.float64)
    lev = np.array(nc.variables['lev'][:], dtype=np.float64) if 'lev' in nc.variables else np.arange(72)

    nc.close()
    return aer_vars, delp, lat, lon, lev


def _compute_pressure_midpoints(delp):
    """Compute mid-layer pressures from DELP.

    MERRA-2 layers: L=0 is TOA, L=71 is surface.
    Edge pressure: p_edge[0] = PTOP, p_edge[k+1] = p_edge[k] + DELP[k]
    Mid-layer pressure: p_mid[k] = (p_edge[k] + p_edge[k+1]) / 2

    Args:
        delp: (nlev, nlat, nlon) or (ntime, nlev, nlat, nlon)

    Returns:
        p_mid: same shape as delp, in Pa
    """
    # Cumulative sum along level axis to get edge pressures
    if delp.ndim == 4:
        lev_axis = 1
    else:
        lev_axis = 0

    p_edge_bottom = PTOP + np.cumsum(delp, axis=lev_axis)
    # p_edge_top[k] = p_edge_bottom[k] - delp[k]
    p_edge_top = p_edge_bottom - delp
    p_mid = (p_edge_top + p_edge_bottom) / 2.0

    return p_mid  # Pa


def _interp_vertical(data, p_mid, target_plev_pa):
    """Interpolate from MERRA-2 hybrid levels to target pressure levels.

    Vectorized log-pressure interpolation using searchsorted for speed.

    Args:
        data: (nlev_m2, nlat_m2, nlon_m2) mixing ratio
        p_mid: (nlev_m2, nlat_m2, nlon_m2) mid-layer pressure in Pa
        target_plev_pa: (nlev_target,) target pressure levels in Pa

    Returns:
        result: (nlev_target, nlat_m2, nlon_m2)
    """
    nlev_m2, nlat, nlon = data.shape
    nlev_tgt = len(target_plev_pa)

    # Work in log-pressure space
    # MERRA-2: lev[0]=TOA (low pressure), lev[-1]=surface (high pressure)
    # p_mid increases from TOA to surface → log_p also increases
    log_p = np.log(np.clip(p_mid, 1e-6, None))    # (nlev_m2, nlat, nlon)
    log_tgt = np.log(target_plev_pa)                # (nlev_tgt,)

    # Reshape for broadcasting: (nlev_tgt, 1, nlat, nlon) vs (nlev_m2, nlat, nlon)
    log_tgt_4d = log_tgt[:, np.newaxis, np.newaxis, np.newaxis]  # (nlev_tgt, 1, nlat, nlon)
    log_p_exp = log_p[np.newaxis, :, :, :]                        # (1, nlev_m2, nlat, nlon)

    # For each target level, find bracket indices in source (per column)
    # diff[t, k, j, i] = log_tgt[t] - log_p[k, j, i]
    # We want the last k where log_p[k] <= log_tgt[t]
    above = log_p_exp <= log_tgt_4d  # (nlev_tgt, nlev_m2, nlat, nlon) bool

    # idx_lo: index of last source level <= target (i.e. highest pressure ≤ target)
    # Use sum of True values = number of source levels below target
    idx_lo = np.sum(above, axis=1) - 1  # (nlev_tgt, nlat, nlon)
    idx_lo = np.clip(idx_lo, 0, nlev_m2 - 2)
    idx_hi = idx_lo + 1

    # Gather log_p and data at bracket indices using advanced indexing
    # jj/ii broadcast: (1,nlat,1) and (1,1,nlon) to (nlev_tgt, nlat, nlon)
    jj = np.arange(nlat)[np.newaxis, :, np.newaxis]  # (1, nlat, 1)
    ii = np.arange(nlon)[np.newaxis, np.newaxis, :]  # (1, 1, nlon)

    lp_lo = log_p[idx_lo, jj, ii]   # (nlev_tgt, nlat, nlon)
    lp_hi = log_p[idx_hi, jj, ii]
    d_lo = data[idx_lo, jj, ii]
    d_hi = data[idx_hi, jj, ii]

    # Linear interpolation weight in log-pressure space
    denom = lp_hi - lp_lo
    safe_denom = np.where(np.abs(denom) < 1e-10, 1.0, denom)
    w = (log_tgt[:, np.newaxis, np.newaxis] - lp_lo) / safe_denom
    w = np.clip(w, 0.0, 1.0)

    result = d_lo + w * (d_hi - d_lo)
    return result  # (nlev_tgt, nlat, nlon)


def _interp_horizontal(data_3d, lat_src, lon_src, lat_tgt, lon_tgt):
    """Bilinear interpolation from MERRA-2 grid to target grid.

    Args:
        data_3d: (nlev, nlat_src, nlon_src)
        lat_src, lon_src: 1D source coordinates
        lat_tgt, lon_tgt: 1D target coordinates

    Returns:
        result: (nlev, nlat_tgt, nlon_tgt)
    """
    nlev = data_3d.shape[0]
    result = np.zeros((nlev, len(lat_tgt), len(lon_tgt)), dtype=np.float64)

    for k in range(nlev):
        interp = RegularGridInterpolator(
            (lat_src, lon_src), data_3d[k],
            method='linear', bounds_error=False, fill_value=None)
        # Build target mesh
        lat_grid, lon_grid = np.meshgrid(lat_tgt, lon_tgt, indexing='ij')
        pts = np.column_stack([lat_grid.ravel(), lon_grid.ravel()])
        result[k] = interp(pts).reshape(len(lat_tgt), len(lon_tgt))

    return result


def load_merra2_aerosol(data_dir, years, month, warm_days,
                        target_plev_hpa, target_lat, target_lon,
                        period='clim'):
    """Load and process MERRA-2 aerosol for a climatological or event period.

    Args:
        data_dir: directory containing M2I3NVAER_{YYYYMMDD}.nc4 files
        years: list of years to average over
        month: month number
        warm_days: list of 0-based day indices
        target_plev_hpa: (nlev,) target pressure levels in hPa (surface→TOA)
        target_lat: (nlat,) target latitude array
        target_lon: (nlon,) target longitude array
        period: 'clim' (average all years) or year number (single year)

    Returns:
        dict of {species: (nlev, nlat, nlon)} on target grid
    """
    target_plev_pa = target_plev_hpa * 100.0  # hPa → Pa
    nlev_tgt = len(target_plev_hpa)
    nlat_tgt = len(target_lat)
    nlon_tgt = len(target_lon)
    shape_tgt = (nlev_tgt, nlat_tgt, nlon_tgt)

    # Species aggregation mapping
    species_map = {
        'bc': ['BCPHILIC', 'BCPHOBIC'],
        'ocphi': ['OCPHILIC'],
        'ocpho': ['OCPHOBIC'],
        'sulf': ['SO4'],
        'ss': ['SS001', 'SS002', 'SS003', 'SS004', 'SS005'],
        'dust': ['DU001', 'DU002', 'DU003', 'DU004', 'DU005'],
    }

    # Accumulate over years
    accum = {sp: np.zeros(shape_tgt, dtype=np.float64) for sp in species_map}
    count = 0

    if period == 'clim':
        proc_years = years
    else:
        proc_years = [int(period)]

    for yr_idx, yr in enumerate(proc_years):
        print(f"    MERRA-2: year {yr} ({yr_idx+1}/{len(proc_years)})", flush=True)
        # Collect files for warm days
        day_data = {sp: [] for sp in species_map}
        n_days_loaded = 0

        for d in warm_days:
            day_num = d + 1  # 0-based → 1-based
            fname = os.path.join(data_dir,
                                 f'M2I3NVAER_{yr}{month:02d}{day_num:02d}.nc4')
            if not os.path.exists(fname):
                # Try alternate pattern
                alt = glob.glob(os.path.join(data_dir,
                                             f'*{yr}{month:02d}{day_num:02d}*.nc4'))
                if alt:
                    fname = alt[0]
                else:
                    continue

            aer_vars, delp, lat_m2, lon_m2, _ = _load_merra2_day(fname)

            # Daily mean over 3-hourly steps (axis=0)
            delp_daily = delp.mean(axis=0)
            p_mid = _compute_pressure_midpoints(delp_daily)

            for sp, m2_vars in species_map.items():
                # Sum contributing MERRA-2 variables, then daily mean
                sp_data = np.zeros_like(delp)
                for mv in m2_vars:
                    if mv in aer_vars:
                        sp_data += aer_vars[mv]
                sp_daily = sp_data.mean(axis=0)  # (nlev, nlat, nlon)

                # Vertical interpolation
                sp_vinterp = _interp_vertical(sp_daily, p_mid,
                                               target_plev_pa)

                # Horizontal interpolation
                # Ensure lat is ascending for RegularGridInterpolator
                if lat_m2[0] > lat_m2[-1]:
                    sp_vinterp = sp_vinterp[:, ::-1, :]
                    lat_sorted = lat_m2[::-1]
                else:
                    lat_sorted = lat_m2

                sp_hinterp = _interp_horizontal(sp_vinterp,
                                                 lat_sorted, lon_m2,
                                                 target_lat, target_lon)
                day_data[sp].append(sp_hinterp)

            n_days_loaded += 1

        if n_days_loaded == 0:
            print(f"    Warning: no MERRA-2 data for {yr}")
            continue

        # Average over warm days for this year
        for sp in species_map:
            if day_data[sp]:
                yr_mean = np.mean(day_data[sp], axis=0)
                accum[sp] += yr_mean
        count += 1

    if count == 0:
        print("  WARNING: No MERRA-2 data loaded. Using zero aerosol.")
        return {sp: np.zeros(shape_tgt, dtype=np.float64) for sp in species_map}

    # Average over years
    result = {}
    for sp in species_map:
        result[sp] = accum[sp] / count
        # Ensure non-negative mixing ratios
        result[sp] = np.maximum(result[sp], 0.0)

    return result
