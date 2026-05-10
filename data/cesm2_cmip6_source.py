"""CMIP6 CESM2 raw output reader + hybrid→plev re-projection.

Loads piControl / abrupt-4xCO2 monthly NetCDF, computes annual climatology
across selected year range, and re-projects cl/clw/cli from 32-layer hybrid
sigma-pressure to user-supplied plev (e.g. CMIP6 plev19).

Convention (input files, CMIP6 standard)
----------------------------------------
- plev (ta/hus): surface→TOA in NetCDF, units Pa
- hybrid (cl/clw/cli): TOA→sfc in NetCDF; pressure of layer k =
                       a[k]*p0 + b[k]*ps(lat,lon)
- lat: -90→90 (S→N), lon: 0→358.75 (W→E)
- time: monthly mean, calendar=noleap, time_origin=0001-01-01

Convention (output, pyCFRAM input)
----------------------------------
NetCDF stores arrays in **surface→TOA** lev order (same as CMIP6 plev). pyCFRAM
applies `[::-1]` internally to flip to TOA→sfc for Fortran processing.
"""
import os
import glob
import numpy as np
from netCDF4 import Dataset


# CMIP6 plev (surface→TOA in NetCDF) for ta/hus — fixed
PLEV19_PA = np.array([100000, 92500, 85000, 70000, 60000, 50000, 40000,
                      30000, 25000, 20000, 15000, 10000, 7000, 5000,
                      3000, 2000, 1000, 500, 100], dtype=np.float64)

# pyCFRAM TOA→sfc convention (reverse of CMIP6 NetCDF order, in hPa)
PLEV19_HPA_TOP_DOWN = (PLEV19_PA[::-1] / 100.0).copy()  # [1, 5, ..., 925, 1000]


def list_files(raw_dir, exp_subdir):
    """Discover all CMIP6 monthly NetCDFs in raw_dir/exp_subdir."""
    files = sorted(glob.glob(os.path.join(raw_dir, exp_subdir, '*.nc')))
    var_files = {}
    for f in files:
        # filename: <var>_Amon_CESM2_<exp>_r1i1p1f1_gn_<period>.nc
        var = os.path.basename(f).split('_')[0]
        var_files[var] = f
    return var_files


def years_to_month_indices(time_var, year_start, year_end):
    """Return slice over months whose YEAR (per noleap calendar) is in
    [year_start, year_end] inclusive.

    CESM2 noleap time stored in days since 0001-01-01. 365 days/year exactly.
    Month index: floor(time / 365 * 12) gives integer year offset directly.
    """
    days = np.asarray(time_var[:], dtype=np.float64)
    # noleap: each year is exactly 365 days. time stored at month midpoint.
    # We can recover (year, month) approx via days / 365.0 → year fraction.
    year_frac = days / 365.0
    year_int = year_frac.astype(int) + 1   # +1 because CMIP6 year origin is yr 0001
    # NB: year_int = 1 corresponds to days 0..364
    # Standard CESM convention: time is mid-month → year = floor(days/365)+1
    mask = (year_int >= year_start) & (year_int <= year_end)
    return np.where(mask)[0]


def annual_climo_from_monthly(field, time_indices):
    """field shape (time, ...) → annual climo, day-weighted.

    NOTE: CMIP6 plev19 data has fillvalue (~1e20) for cells at pressure
    levels below local surface (e.g., 1000 hPa over high orography). Some
    grid points have valid data in some months and fillvalue in others
    (seasonal ps variation). A simple weighted mean would produce polluted
    values like 8e16 (fillvalue × month_fraction).

    Fix: convert fillvalues to NaN first, then per-cell day-weighted mean
    over **only valid months**. Cells with no valid data in any month
    remain NaN (handled later by mask_subsurface_layers.py).
    """
    NOLEAP_DAYS = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31],
                           dtype=np.float64)
    n = len(time_indices)
    if n == 0:
        raise ValueError('No time records selected')

    sub = np.asarray(field[time_indices], dtype=np.float64)
    # Mask fillvalues (CMIP6 typical 1e20; covers anything > 1e15 safely)
    sub = np.where(np.abs(sub) > 1e15, np.nan, sub)

    months = (time_indices % 12)
    w = NOLEAP_DAYS[months]
    w_b = w.reshape((n,) + (1,) * (sub.ndim - 1))   # broadcast shape

    # Per-cell weighted mean ignoring NaN: sum(w*v where valid) / sum(w where valid)
    valid = ~np.isnan(sub)
    num = np.nansum(np.where(valid, w_b * sub, 0.0), axis=0)
    den = np.sum(np.where(valid, w_b, 0.0), axis=0)
    with np.errstate(divide='ignore', invalid='ignore'):
        out = np.where(den > 0, num / den, np.nan)
    return out


def load_climo_pres(raw_dir, exp_subdir, year_start, year_end):
    """Load all pres-level + 2D variables, return (data_dict, lat, lon, plev).

    Returns dict with keys: ta, hus, cl, clw, cli, ts, ps, rsdt, rsds, rsus,
                            hfls, hfss, plus hybrid: a, b, p0
    Each value is 3D (lev, lat, lon) for upper-air, 2D (lat, lon) for surface.
    """
    files = list_files(raw_dir, exp_subdir)
    print('  Loading climo year %d-%d from %s/' % (year_start, year_end, exp_subdir))

    # Time indices come from any variable (all aligned)
    f = Dataset(files['ta'])
    idx = years_to_month_indices(f.variables['time'], year_start, year_end)
    print('  selected %d months' % len(idx))
    lat = np.array(f.variables['lat'][:])
    lon = np.array(f.variables['lon'][:])
    plev = np.array(f.variables['plev'][:])
    f.close()

    out = {'lat': lat, 'lon': lon, 'plev_pres': plev}

    # Variables on plev grid (ta, hus): shape (nlev_plev=19, nlat, nlon)
    for var in ('ta', 'hus'):
        f = Dataset(files[var])
        out[var] = annual_climo_from_monthly(f.variables[var], idx)
        f.close()

    # Variables on hybrid grid (cl, clw, cli): shape (nlev_hyb=32, nlat, nlon)
    f = Dataset(files['cl'])
    out['hybrid_a'] = np.array(f.variables['a'][:], dtype=np.float64)
    out['hybrid_b'] = np.array(f.variables['b'][:], dtype=np.float64)
    out['hybrid_p0'] = float(f.variables['p0'][...])
    out['cl'] = annual_climo_from_monthly(f.variables['cl'], idx) / 100.0  # %→fraction
    f.close()

    f = Dataset(files['clw'])
    out['clw'] = annual_climo_from_monthly(f.variables['clw'], idx)
    f.close()

    f = Dataset(files['cli'])
    out['cli'] = annual_climo_from_monthly(f.variables['cli'], idx)
    f.close()

    # 2D surface fields. huss = 2m specific humidity (CMIP6 standard, kg/kg);
    # used by Fu RT for ph(nv1) — apple-to-apple OLD CFRAM raw/CFRAM.zip
    # GW-base.f L322-330 reads huss_base.dat for the surface row of /atmosp/.
    for var in ('ts', 'ps', 'rsdt', 'rsds', 'rsus', 'hfls', 'hfss', 'huss'):
        f = Dataset(files[var])
        out[var] = annual_climo_from_monthly(f.variables[var], idx)
        f.close()

    return out


def hybrid_to_plev(field_hyb, a, b, p0, ps_2d, plev_target_pa):
    """Re-project field on hybrid sigma-pressure to fixed pressure levels.

    Args:
        field_hyb: shape (nlev_hyb, nlat, nlon), TOA→sfc order
        a, b: shape (nlev_hyb,), hybrid coefficients
        p0: scalar reference pressure (Pa)
        ps_2d: shape (nlat, nlon), surface pressure (Pa)
        plev_target_pa: shape (nlev_target,), target pressure levels (Pa)
                        — order doesn't matter, output matches input order

    Returns:
        field_plev: shape (nlev_target, nlat, nlon)

    Method: log-p linear interpolation per column. Above hybrid TOA: 0.
    Below hybrid bottom: extend bottom value to surface.
    """
    nlev_hyb, nlat, nlon = field_hyb.shape
    nlev_target = len(plev_target_pa)
    out = np.zeros((nlev_target, nlat, nlon), dtype=np.float64)

    # log of target levels (stays constant per column)
    log_pt = np.log(plev_target_pa)

    for j in range(nlat):
        for i in range(nlon):
            # Compute hybrid layer pressures at this column
            p_hyb = a * p0 + b * ps_2d[j, i]   # shape (nlev_hyb,), TOA→sfc
            log_phyb = np.log(p_hyb)
            field_col = field_hyb[:, j, i]      # shape (nlev_hyb,)

            # np.interp requires increasing x (log_phyb is already increasing
            # since p_hyb goes TOA→sfc = small→large, log strictly increasing).
            # For target levels above hybrid top: extrapolate to 0.
            # For target levels below bottom: extrapolate to bottom value.
            interp_vals = np.interp(log_pt, log_phyb, field_col,
                                    left=0.0, right=field_col[-1])
            out[:, j, i] = interp_vals

    return out


def hybrid_to_plev_mass_conserving(field_hyb, a, b, p0, ps_2d, plev_target_pa):
    """Mass-conserving re-projection of mixing ratio from hybrid to plev.

    Builds cumulative mass M(p) = ∫_0^p f * dp on the hybrid grid (trapezoidal),
    interpolates M at target plev boundaries, and differences to extract
    layer-mean mixing ratios. By construction, ∑(target_field·Δp_target)
    over a column equals ∑(hybrid_field·Δp_hybrid) — total cloud mass conserved.

    Replaces the lossy single-point log-linear sampling in `hybrid_to_plev`,
    which was found (via scripts/diag_cloud_column.py) to lose ~5.5% of
    column-integrated liquid water on average, with up to 24% loss in
    boundary-layer cloud regions. Use this for cl, clw, cli; hus/ta are
    already on plev directly from CMIP6 CMOR and don't need re-projection.

    Args:
        field_hyb: (nlev_hyb, nlat, nlon) mixing ratio on hybrid (TOA→sfc)
        a, b: hybrid coefficients (TOA→sfc)
        p0: reference pressure (Pa)
        ps_2d: (nlat, nlon) surface pressure (Pa)
        plev_target_pa: target plev (Pa), TOA→sfc ordering (ascending)

    Returns:
        field_plev: (nlev_target, nlat, nlon) layer-mean mixing ratio,
                    TOA→sfc ordering. The value at index k represents the
                    mass-mean mixing ratio for the layer between
                    plev_target[k] (top, smaller p) and plev_target[k+1]
                    (bottom, larger p). For the LAST index, the layer is
                    between plev_target[-2] and ps (capped).
    """
    nlev_hyb = len(a)
    nlev_t = len(plev_target_pa)
    nlat, nlon = ps_2d.shape

    out = np.zeros((nlev_t, nlat, nlon), dtype=np.float64)

    for j in range(nlat):
        for i in range(nlon):
            ps = ps_2d[j, i]
            p_h = a * p0 + b * ps              # (nlev_hyb,) TOA→sfc, ascending
            f_h = field_hyb[:, j, i]
            f_h = np.where(np.isnan(f_h), 0.0, f_h)

            # Cumulative mass from TOA via trapezoidal. Anchors: M(p=0)=0,
            # M(p=ps) = M(p_h[-1]) + f_h[-1] * (ps - p_h[-1]).
            M_hyb = np.zeros(nlev_hyb)
            for kh in range(1, nlev_hyb):
                dp = p_h[kh] - p_h[kh-1]
                f_avg = 0.5 * (f_h[kh] + f_h[kh-1])
                M_hyb[kh] = M_hyb[kh-1] + f_avg * dp
            M_at_ps = M_hyb[-1] + f_h[-1] * (ps - p_h[-1])

            p_anchor = np.concatenate([[0.0], p_h, [ps]])
            M_anchor = np.concatenate([[0.0], M_hyb, [M_at_ps]])

            # Target layer k is between plev_target[k] (top) and plev_target[k+1]
            # (bottom). Last index: layer between plev_target[-1] and ps.
            for kt in range(nlev_t):
                p_top = plev_target_pa[kt]
                if kt < nlev_t - 1:
                    p_bot_nominal = plev_target_pa[kt + 1]
                else:
                    p_bot_nominal = ps    # last index: integrate to surface
                p_bot = min(p_bot_nominal, ps)
                if p_bot <= p_top:
                    out[kt, j, i] = 0.0   # subsurface or zero-thickness
                    continue
                M_top = np.interp(p_top, p_anchor, M_anchor)
                M_bot = np.interp(p_bot, p_anchor, M_anchor)
                out[kt, j, i] = (M_bot - M_top) / (p_bot - p_top)

    return out


def compute_albedo(rsus, rsds):
    """Surface albedo = rsus / rsds, clipped [0, 1]. Where rsds≈0 (polar
    night), set albedo=0 (no SW input → albedo undefined, doesn't matter)."""
    with np.errstate(divide='ignore', invalid='ignore'):
        alb = np.where(rsds > 1.0, rsus / rsds, 0.0)
    return np.clip(alb, 0.0, 1.0)


def reorder_for_pycfram_input(arr_top_down):
    """If input is TOA→sfc (numpy convention used internally in this module),
    flip to sfc→TOA for pyCFRAM input NetCDF (which expects sfc→TOA, with
    `[::-1]` applied inside run_parallel_python.py to recover TOA→sfc).
    """
    return arr_top_down[::-1]
