#!/usr/bin/env python3
"""Diagnostic: compare cloud water/ice column-integrated mass between
CMIP6 hybrid (32 levels) and pyCFRAM plev (19 levels) representations.

Tests candidate 2 of the WV/CO2/CLD dT under-estimation hypothesis: does the
hybrid → plev re-projection in data/cesm2_cmip6_source.py:hybrid_to_plev()
lose cloud column mass?

Reads raw CMIP6:
  - cl, clw, cli, ps for piControl 500-599 climo
  - The hybrid coefficients a, b, p0 from cl file
Reads pyCFRAM input:
  - cases/cesm2_4xco2_official/input/base_pres.nc — cliq, cice on plev19
  - cases/cesm2_4xco2_official/input/base_surf.nc — ps

Computes for each cell (lat, lon):
  W_hybrid_clw = sum_k clw_hyb(k) * dp_hyb(k) / g    (kg/m^2)
  W_plev_clw   = sum_k cliq_plev(k) * dp_plev(k) / g  (kg/m^2)

Compares ratio W_plev / W_hybrid. Should be ~1.0 if mass-conserving.

Usage:  python3 scripts/diag_cloud_column.py
"""
import os
import sys
import numpy as np
from netCDF4 import Dataset

ROOT = '/home/lzhenn/work/ust-jumper/pyCFRAM'
G = 9.80616  # m/s^2

# Match build_cesm2_official.py climo period
CMIP6_RAW = os.path.join(ROOT, 'raw_data/cesm2_cmip6/piControl')
PYCFRAM_INPUT = os.path.join(ROOT, 'cases/cesm2_4xco2_official/input')
YEAR_START, YEAR_END = 500, 599


def years_to_month_indices(time_var, y0, y1):
    days = np.asarray(time_var[:], dtype=np.float64)
    year = (days / 365.0).astype(int) + 1
    return np.where((year >= y0) & (year <= y1))[0]


def annual_climo(field, idx):
    """NaN-aware day-weighted annual climatology over selected months."""
    NOLEAP = np.array([31,28,31,30,31,30,31,31,30,31,30,31], dtype=np.float64)
    sub = np.asarray(field[idx], dtype=np.float64)
    sub = np.where(np.abs(sub) > 1e15, np.nan, sub)
    months = (idx % 12)
    w = NOLEAP[months]
    w_b = w.reshape((len(idx),) + (1,) * (sub.ndim - 1))
    weighted = sub * w_b
    valid = ~np.isnan(sub)
    den = (w_b * valid).sum(axis=0)
    num = np.where(valid, weighted, 0.0).sum(axis=0)
    return np.where(den > 0, num / den, np.nan)


def hybrid_column_mass(clw_hyb, a, b, p0, ps_2d):
    """Compute column-integrated cloud water mass (kg/m^2) on hybrid grid.

    Layer k spans hybrid level k to k+1. Uses TOA→sfc ordering as stored.
    Pressure at hybrid level k: p(k) = a(k)*p0 + b(k)*ps.
    Layer thickness: dp(k) = |p(k+1) - p(k)|.
    Column mass: W = sum_k clw(k) * dp(k) / g.
    """
    # CMIP6 stores hybrid TOA→sfc (a decreases, b increases). Verify.
    nlev_hyb = len(a)
    nlat, nlon = ps_2d.shape
    # Build pressure per cell, per hybrid level (nlev, nlat, nlon)
    a_b = a.reshape(nlev_hyb, 1, 1)
    b_b = b.reshape(nlev_hyb, 1, 1)
    p_hyb = a_b * p0 + b_b * ps_2d[None, :, :]   # (nlev, nlat, nlon), Pa
    # Layer thickness: |p(k+1) - p(k)|. Need TOA→sfc to be increasing or
    # sfc→TOA to be decreasing — either way, |diff| works.
    dp = np.abs(np.diff(p_hyb, axis=0))   # (nlev-1, nlat, nlon)
    # Treat clw_hyb as defined per LEVEL (not layer). Use trapezoidal:
    # layer k mass ≈ 0.5*(clw(k)+clw(k+1)) * dp(k). For simplicity use
    # mid-point clw (still mass-conserving in trapezoidal sense).
    clw_layer = 0.5 * (clw_hyb[:-1] + clw_hyb[1:])   # (nlev-1, nlat, nlon)
    W = np.nansum(clw_layer * dp / G, axis=0)        # (nlat, nlon), kg/m^2
    return W


def plev_column_mass(clw_plev, plev_pa, ps_2d):
    """Compute column-integrated cloud water on pressure-level grid.

    plev_pa: (nlev,) shape, in surface→TOA NetCDF order (1000, 925, ..., 1, Pa).
    clw_plev: (nlev, nlat, nlon), same lev order.

    For each cell, only include layers where p(k) <= ps. Cap the lowest
    layer to ps. Layer k mass ≈ clw(k) * dp(k) / g where dp(k) is the
    layer thickness, summed over atmospheric layers only.
    """
    nlev = len(plev_pa)
    nlat, nlon = ps_2d.shape
    # Treat plev_pa as level boundaries, not layer centers (CMIP6 standard).
    # For cell-by-cell: build effective level list with surface inserted.
    # Simpler: use trapezoidal between adjacent plev levels for atm layers.
    # In sfc→TOA order, layer between plev[k] and plev[k+1] has thickness
    # |plev[k+1] - plev[k]| and clw_layer = 0.5*(clw[k]+clw[k+1]).
    # Mask layers where plev[k] > ps (subsurface) or plev[k+1] > ps.
    W = np.zeros((nlat, nlon))
    for k in range(nlev - 1):
        p_lower = max(plev_pa[k], plev_pa[k+1])     # the larger pressure (closer to surface)
        # Subsurface mask: layer is subsurface if both edges above ps
        # (i.e., both have p > ps; in pressure terms p>ps means below sfc)
        # We assume sfc→TOA order so plev[k] > plev[k+1].
        # Layer is "above surface" iff plev[k] (lower edge, higher p) <= ps.
        in_atm = (plev_pa[k] <= ps_2d).astype(np.float64)   # 1 if layer is above sfc
        dp = abs(plev_pa[k+1] - plev_pa[k])
        clw_layer = 0.5 * (clw_plev[k] + clw_plev[k+1])
        # NaN-safe
        contrib = np.where(np.isnan(clw_layer), 0.0, clw_layer) * dp / G * in_atm
        W += contrib
    return W


def main():
    print('=== Loading CMIP6 hybrid raw (piControl) ===')
    f_cl = Dataset(os.path.join(CMIP6_RAW, 'cl_Amon_CESM2_piControl_r1i1p1f1_gn_050001-059912.nc'))
    a = np.array(f_cl.variables['a'][:], dtype=np.float64)
    b = np.array(f_cl.variables['b'][:], dtype=np.float64)
    p0 = float(f_cl.variables['p0'][...])
    print(f'  hybrid nlev={len(a)}, p0={p0} Pa')
    print(f'  a[0]={a[0]:.4e} (TOA), a[-1]={a[-1]:.4e} (sfc)')
    print(f'  b[0]={b[0]:.4e} (TOA), b[-1]={b[-1]:.4e} (sfc)')

    idx = years_to_month_indices(f_cl.variables['time'], YEAR_START, YEAR_END)
    print(f'  selected {len(idx)} months ({YEAR_START}-{YEAR_END})')

    # Read clw, cli on hybrid (lazy load via climo)
    f_clw = Dataset(os.path.join(CMIP6_RAW, 'clw_Amon_CESM2_piControl_r1i1p1f1_gn_050001-059912.nc'))
    f_cli = Dataset(os.path.join(CMIP6_RAW, 'cli_Amon_CESM2_piControl_r1i1p1f1_gn_050001-059912.nc'))
    f_ps  = Dataset(os.path.join(CMIP6_RAW, 'ps_Amon_CESM2_piControl_r1i1p1f1_gn_050001-059912.nc'))

    print('\nComputing hybrid annual climo (clw, cli, ps)...')
    clw_hyb = annual_climo(f_clw.variables['clw'], idx)   # kg/kg, (32, 192, 288)
    cli_hyb = annual_climo(f_cli.variables['cli'], idx)
    ps_hyb  = annual_climo(f_ps.variables['ps'],  idx)    # Pa
    print(f'  clw_hyb shape={clw_hyb.shape}, max={np.nanmax(clw_hyb):.2e}')
    print(f'  cli_hyb shape={cli_hyb.shape}, max={np.nanmax(cli_hyb):.2e}')
    print(f'  ps_hyb min/max: {np.nanmin(ps_hyb)/100:.1f} / {np.nanmax(ps_hyb)/100:.1f} hPa')

    f_cl.close(); f_clw.close(); f_cli.close(); f_ps.close()

    # Convert clw_hyb to mass mixing ratio if needed (CMIP6 spec is kg/kg)
    print('\n=== Loading pyCFRAM plev input (cesm2_4xco2_official base_pres.nc) ===')
    nc_b = Dataset(os.path.join(PYCFRAM_INPUT, 'base_pres.nc'))
    nc_s = Dataset(os.path.join(PYCFRAM_INPUT, 'base_surf.nc'))
    lev_hpa = np.array(nc_b.variables['lev'][:])    # (19,) sfc→TOA, hPa
    cliq_plev = np.array(nc_b.variables['cliq'][0]) # (19, 192, 288), kg/kg, sfc→TOA
    cice_plev = np.array(nc_b.variables['cice'][0])
    ps_plev   = np.array(nc_s.variables['ps'][0])   # Pa
    nc_b.close(); nc_s.close()
    print(f'  plev nlev={len(lev_hpa)}, lev[0]={lev_hpa[0]}, lev[-1]={lev_hpa[-1]} hPa')
    print(f'  cliq max={np.nanmax(cliq_plev):.2e}')
    print(f'  cice max={np.nanmax(cice_plev):.2e}')

    plev_pa = lev_hpa * 100.0   # hPa → Pa

    print('\n=== Computing column-integrated cloud water mass (kg/m^2) ===')
    W_clw_hyb = hybrid_column_mass(clw_hyb, a, b, p0, ps_hyb)
    W_cli_hyb = hybrid_column_mass(cli_hyb, a, b, p0, ps_hyb)
    W_clw_plev = plev_column_mass(cliq_plev, plev_pa, ps_plev)
    W_cli_plev = plev_column_mass(cice_plev, plev_pa, ps_plev)

    print('\n--- Cloud LIQUID water column (W_clw, kg/m^2) ---')
    print(f'  hybrid  global mean: {np.nanmean(W_clw_hyb):.4f}  median: {np.nanmedian(W_clw_hyb):.4f}')
    print(f'  plev    global mean: {np.nanmean(W_clw_plev):.4f}  median: {np.nanmedian(W_clw_plev):.4f}')
    rel_diff = (W_clw_plev - W_clw_hyb) / np.where(W_clw_hyb > 1e-6, W_clw_hyb, np.nan)
    print(f'  ratio plev/hybrid: mean={np.nanmean(W_clw_plev/np.where(W_clw_hyb>1e-6,W_clw_hyb,np.nan)):.3f}, '
          f'median={np.nanmedian(W_clw_plev/np.where(W_clw_hyb>1e-6,W_clw_hyb,np.nan)):.3f}')
    print(f'  rel diff (plev-hyb)/hyb: mean={np.nanmean(rel_diff)*100:+.1f}%, '
          f'p10/p50/p90={np.nanpercentile(rel_diff,10)*100:+.1f}%/{np.nanpercentile(rel_diff,50)*100:+.1f}%/{np.nanpercentile(rel_diff,90)*100:+.1f}%')

    print('\n--- Cloud ICE water column (W_cli, kg/m^2) ---')
    print(f'  hybrid  global mean: {np.nanmean(W_cli_hyb):.4f}  median: {np.nanmedian(W_cli_hyb):.4f}')
    print(f'  plev    global mean: {np.nanmean(W_cli_plev):.4f}  median: {np.nanmedian(W_cli_plev):.4f}')
    rel_diff = (W_cli_plev - W_cli_hyb) / np.where(W_cli_hyb > 1e-6, W_cli_hyb, np.nan)
    print(f'  ratio plev/hybrid: mean={np.nanmean(W_cli_plev/np.where(W_cli_hyb>1e-6,W_cli_hyb,np.nan)):.3f}, '
          f'median={np.nanmedian(W_cli_plev/np.where(W_cli_hyb>1e-6,W_cli_hyb,np.nan)):.3f}')
    print(f'  rel diff (plev-hyb)/hyb: mean={np.nanmean(rel_diff)*100:+.1f}%, '
          f'p10/p50/p90={np.nanpercentile(rel_diff,10)*100:+.1f}%/{np.nanpercentile(rel_diff,50)*100:+.1f}%/{np.nanpercentile(rel_diff,90)*100:+.1f}%')

    print('\nDONE.')


if __name__ == '__main__':
    main()
