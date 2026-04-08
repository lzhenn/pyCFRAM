#!/usr/bin/env python3
"""Additivity / Linearization Error Validation (Task 2 of Phase 1 verification suite).

Reads cfram_result.nc and checks that the sum of all decomposed temperature
contributions equals the observed temperature change:

    dT_co2 + dT_q + dT_ts + dT_o3 + dT_solar + dT_albedo
    + dT_cloud + dT_aerosol + dT_atmdyn + dT_sfcdyn  ≈  dT_observed

Note on sub-decompositions:
  sfcdyn = ocndyn + lhflx + shflx   (sub-decompositions, NOT added separately)
  aerosol = bc + oc + sulf + seas + dust  (sub-decompositions, NOT added separately)

Two thresholds are checked:
  - STRICT : max|residual| < 0.1 K  (should hold for forcing-based design)
  - LIBERAL: max|residual| < 0.5 K  (fallback; still scientifically acceptable)

Usage:
    python3 scripts/validate_additivity.py --case eh13 [--plot]

Output:
    - Console: per-level and global residual statistics, PASS/FAIL
    - (optional) PNG: residual spatial map + per-level profile
"""

import os, sys, argparse
import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case, get_plev, PROJECT_ROOT

FORTRAN_PLEV = get_plev()
NLEV = len(FORTRAN_PLEV)

# Primary additive terms (radiative + top-level dynamic)
# NOTE: 'ts' (surface skin temperature LW re-emission) is EXCLUDED from the sum.
# dT_ts represents the radiative forcing of ΔTs on the atmosphere, which is already
# embedded in the Planck matrix (∂R/∂T includes the surface layer). Including it
# would double-count the surface LW feedback. Wu et al. (2025) correctly excludes
# this term from their 10-term decomposition.
RAD_TERMS = ['co2', 'q', 'o3', 'solar', 'albedo', 'cloud', 'aerosol']
DYN_TERMS = ['atmdyn', 'sfcdyn']
ALL_TERMS = RAD_TERMS + DYN_TERMS

# Sub-decomposition terms (informational only, not summed)
SUB_TERMS = ['ocndyn', 'lhflx', 'shflx', 'bc', 'oc', 'sulf', 'seas', 'dust']


def load_nc(path, varname):
    nc = Dataset(path, 'r')
    if varname not in nc.variables:
        nc.close()
        return None
    data = np.array(nc.variables[varname][:], dtype=np.float64)
    nc.close()
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--case', required=True, help='Case name (e.g. eh13)')
    parser.add_argument('--plot', action='store_true', help='Save diagnostic plots')
    parser.add_argument('--roi', action='store_true',
                        help='Restrict statistics to key_region defined in case.yaml')
    args = parser.parse_args()

    cfg = load_case(args.case)
    OUTDIR = cfg['_output_dir']
    outfile = os.path.join(OUTDIR, 'cfram_result.nc')

    if not os.path.exists(outfile):
        print(f'ERROR: output file not found: {outfile}')
        print('Run run_parallel_python.py first to generate it.')
        sys.exit(1)

    print(f'Loading: {outfile}')

    # Auto-detect all dT_* terms in file, then decide what to sum.
    # Primary additive set: rad terms + atmdyn + sfcdyn + lhflx + shflx
    # (sfcdyn = ocndyn+lhflx+shflx for new code; lhflx/shflx separate for old code)
    CANDIDATE_TERMS = RAD_TERMS + ['atmdyn', 'sfcdyn', 'lhflx', 'shflx']

    dT = {}
    for t in CANDIDATE_TERMS + ['observed']:
        data = load_nc(outfile, 'dT_' + t)
        if data is not None:
            dT[t] = data

    # Determine correct sum: avoid double-counting.
    # New code (forcing-based, commit a959143+): outputs dT_ocndyn, and
    #   sfcdyn = ocndyn + lhflx + shflx  → sum only sfcdyn (not lhflx/shflx)
    # Old code (residual-based): no dT_ocndyn, sfcdyn/lhflx/shflx independent
    #   → must sum all three.
    # Detection: presence of dT_ocndyn is the definitive indicator.
    has_ocndyn = load_nc(outfile, 'dT_ocndyn') is not None
    if 'sfcdyn' in dT and has_ocndyn:
        SUM_TERMS = RAD_TERMS + ['atmdyn', 'sfcdyn']
        print('Detected: new forcing-based output (sfcdyn includes lhflx+shflx)')
    else:
        SUM_TERMS = [t for t in CANDIDATE_TERMS if t in dT]
        print(f'Using available terms: {SUM_TERMS}')

    available = [t for t in SUM_TERMS if t in dT]
    missing = [t for t in SUM_TERMS if t not in dT]
    if missing:
        print(f'WARNING: missing variables in sum: {missing}')

    if not available:
        print('ERROR: no terms available to sum.')
        sys.exit(1)

    print(f'Summing terms: {available}')

    if 'observed' not in dT:
        print('ERROR: dT_observed not in output file.')
        sys.exit(1)

    # Compute sum of available terms
    dT_sum = np.zeros_like(dT['observed'])
    for t in available:
        dT_sum += dT[t]

    residual = dT_sum - dT['observed']

    # Mask fill values (-999)
    valid = (np.abs(dT['observed']) < 900) & (np.abs(dT_sum) < 900)
    residual_valid = np.where(valid, residual, np.nan)

    # ── ROI mask (optional) ───────────────────────────────────────────────────
    roi_mask = None
    roi_label = 'full domain'
    if args.roi:
        kr = cfg.get('plot', {}).get('key_region', None)
        if kr is None:
            print('WARNING: --roi specified but no key_region in case.yaml; using full domain.')
        else:
            # Read lat/lon from input file
            input_file = cfg['input'].get('base_pres')
            nc_in = Dataset(input_file, 'r')
            lat = np.array(nc_in.variables['lat'][:])
            lon = np.array(nc_in.variables['lon'][:])
            nc_in.close()
            lat_min, lat_max = kr['lat']
            lon_min, lon_max = kr['lon']
            lat_idx = np.where((lat >= lat_min) & (lat <= lat_max))[0]
            lon_idx = np.where((lon >= lon_min) & (lon <= lon_max))[0]
            # Build 2D spatial mask broadcastable to (nlev, nlat, nlon)
            roi_2d = np.zeros((len(lat), len(lon)), dtype=bool)
            roi_2d[np.ix_(lat_idx, lon_idx)] = True
            roi_mask = roi_2d[np.newaxis, :, :]  # (1, nlat, nlon)
            roi_label = (f'ROI {lat_min}–{lat_max}°N, {lon_min}–{lon_max}°E '
                         f'({len(lat_idx)}×{len(lon_idx)} pts)')
            residual_valid = np.where(roi_mask, residual_valid, np.nan)
            print(f'ROI mask applied: {roi_label}')

    # ── Statistics ────────────────────────────────────────────────────────────
    max_abs  = np.nanmax(np.abs(residual_valid))
    mean_abs = np.nanmean(np.abs(residual_valid))
    p95      = np.nanpercentile(np.abs(residual_valid), 95)
    p99      = np.nanpercentile(np.abs(residual_valid), 99)

    # Wu et al. style: signed domain-mean of sum and obs (surface level only)
    sfc_idx = residual_valid.shape[0] - 1   # last level = surface in file
    sum_sfc_mean = np.nanmean(np.where(
        valid[sfc_idx] & (roi_mask[0] if roi_mask is not None else True),
        dT_sum[sfc_idx], np.nan))
    obs_sfc_mean = np.nanmean(np.where(
        valid[sfc_idx] & (roi_mask[0] if roi_mask is not None else True),
        dT['observed'][sfc_idx], np.nan))
    signed_res   = sum_sfc_mean - obs_sfc_mean

    print(f'\nAdditivity residual  (sum({", ".join(available)}) - dT_observed)')
    print(f'  Region          : {roi_label}')
    print(f'  --- Wu et al. style (surface, signed domain-mean) ---')
    print(f'  Sum  (sfc mean) : {sum_sfc_mean:+.4f} K')
    print(f'  Obs  (sfc mean) : {obs_sfc_mean:+.4f} K')
    print(f'  Signed residual : {signed_res:+.4f} K')
    print(f'  --- Full-column absolute stats ---')
    print(f'  max |residual|  : {max_abs:.4f} K')
    print(f'  mean|residual|  : {mean_abs:.4f} K')
    print(f'  P95 |residual|  : {p95:.4f} K')
    print(f'  P99 |residual|  : {p99:.4f} K')

    # Per-level breakdown
    nlev_total = residual_valid.shape[0]
    print(f'\nPer-level max|residual| (K):')
    print(f'  {"Level":<8}  {"P(hPa)":<9}  {"max|res|":<10}  {"mean|res|":<10}')
    print(f'  {"-"*45}')
    # Pressure coordinate: first NLEV levels are atmosphere, last is surface
    plev_labels = list(FORTRAN_PLEV[::-1]) + ['sfc']  # surface->TOA order in file
    for k in range(nlev_total):
        lev_res = np.abs(residual_valid[k])
        if np.all(np.isnan(lev_res)):
            continue
        lv_max = np.nanmax(lev_res)
        lv_mean = np.nanmean(lev_res)
        if lv_max > 0.01:  # only print levels with notable residual
            try:
                plabel = f'{plev_labels[k]:.1f}' if k < len(plev_labels)-1 else 'sfc'
            except Exception:
                plabel = str(k)
            print(f'  {k:<8}  {plabel:<9}  {lv_max:<10.4f}  {lv_mean:<10.4f}')

    # ── Sub-decomposition consistency ─────────────────────────────────────────
    print('\nSub-decomposition checks:')
    # sfcdyn = ocndyn + lhflx + shflx
    sub_sfc = ['ocndyn', 'lhflx', 'shflx']
    avail_sfc = [t for t in sub_sfc if t in
                 {k: None for k in [load_nc(outfile, 'dT_'+s) is not None and s
                                    for s in sub_sfc] if k}]
    # simpler: just try to load and sum
    sub_sum = None
    for t in sub_sfc:
        data = load_nc(outfile, 'dT_' + t)
        if data is not None:
            sub_sum = data if sub_sum is None else sub_sum + data
    if sub_sum is not None and 'sfcdyn' in dT:
        sub_res = np.where(valid, sub_sum - dT['sfcdyn'], np.nan)
        sub_max = np.nanmax(np.abs(sub_res))
        sub_mean = np.nanmean(np.abs(sub_res))
        print(f'  ocndyn+lhflx+shflx vs sfcdyn: max={sub_max:.4f}K  mean={sub_mean:.4f}K')
    else:
        print('  sfcdyn sub-check: skipped (missing variables)')

    # ── PASS/FAIL ─────────────────────────────────────────────────────────────
    print('\n' + '='*55)
    if missing:
        print(f'NOTE: {len(missing)} terms missing ({missing}) - sum is partial')
    strict_pass = max_abs < 0.1
    liberal_pass = max_abs < 0.5
    if strict_pass:
        print(f'Additivity: PASS (strict, max|res|={max_abs:.4f}K < 0.1K)')
    elif liberal_pass:
        print(f'Additivity: PASS-liberal (max|res|={max_abs:.4f}K < 0.5K)')
        print('  Note: residual > 0.1K; may reflect non-linearity or missing terms')
    else:
        print(f'Additivity: FAIL (max|res|={max_abs:.4f}K ≥ 0.5K)')
    print('='*55)

    # ── Optional plot ─────────────────────────────────────────────────────────
    if args.plot:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt

            fig, axes = plt.subplots(1, 2, figsize=(14, 5))

            # Left: spatial map of max|residual| across levels
            max_map = np.nanmax(np.abs(residual_valid), axis=0)
            im = axes[0].imshow(max_map, aspect='auto', origin='lower',
                                cmap='Reds', vmin=0)
            axes[0].set_title(f'max_lev |residual| (K)  [{args.case}]')
            axes[0].set_xlabel('lon index'); axes[0].set_ylabel('lat index')
            plt.colorbar(im, ax=axes[0], label='K')

            # Right: per-level profile of max and mean |residual|
            lev_max_profile = np.nanmax(np.abs(residual_valid),
                                        axis=(1, 2))
            lev_mean_profile = np.nanmean(np.abs(residual_valid),
                                          axis=(1, 2))
            lev_idx = np.arange(nlev_total)
            axes[1].plot(lev_max_profile, lev_idx, label='max', color='red')
            axes[1].plot(lev_mean_profile, lev_idx, label='mean', color='blue',
                         linestyle='--')
            axes[1].axvline(0.1, color='gray', lw=0.8, linestyle=':', label='0.1K threshold')
            axes[1].set_xlabel('|residual| (K)')
            axes[1].set_ylabel('level index (0=TOA)')
            axes[1].set_title('Per-level residual profile')
            axes[1].invert_yaxis()
            axes[1].legend()

            plt.tight_layout()
            plot_path = os.path.join(OUTDIR, 'validate_additivity.png')
            plt.savefig(plot_path, dpi=120)
            print(f'Plot saved: {plot_path}')
        except ImportError:
            print('matplotlib not available, skipping plot.')

    sys.exit(0 if liberal_pass else 1)


if __name__ == '__main__':
    main()
