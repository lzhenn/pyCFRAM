#!/usr/bin/env python3
"""Layer-1 Spatial Correlation Analysis of Additivity Residual.

For each grid column, computes residual = sum(dT_k) - dT_observed.
Then asks: which physical processes (|dT_k|) spatially predict the residual?

Method:
  - For each pressure level: Pearson r(|residual|, |dT_k|) across lat-lon
  - Also test cross-product term |dT_q| * |dT_cloud| (WV×Cloud nonlinearity)
  - Vertical profile of correlation for each term
  - Summary: domain-mean |residual| at each level, top-3 correlated terms

Usage:
    python3 scripts/analyze_residual_spatial.py --case eh13 [--plot]
"""

import os, sys, argparse
import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case, get_plev, PROJECT_ROOT

FORTRAN_PLEV = get_plev()   # TOA → surface, hPa
NLEV = len(FORTRAN_PLEV)

RAD_TERMS  = ['co2', 'q', 'ts', 'o3', 'solar', 'albedo', 'cloud', 'aerosol']
DYN_TERMS  = ['atmdyn', 'sfcdyn']
ALL_TERMS  = RAD_TERMS + DYN_TERMS


def load_nc(path, varname):
    nc = Dataset(path, 'r')
    if varname not in nc.variables:
        nc.close()
        return None
    data = np.array(nc.variables[varname][:], dtype=np.float64)
    nc.close()
    return data


def pearson_r_spatial(x, y):
    """Pearson r between 2-D arrays x and y over valid (non-nan) points."""
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 5:
        return np.nan
    xv, yv = x[mask], y[mask]
    xv -= xv.mean(); yv -= yv.mean()
    denom = np.sqrt((xv**2).sum() * (yv**2).sum())
    return float((xv * yv).sum() / denom) if denom > 0 else np.nan


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--case', required=True, help='Case name (e.g. eh13)')
    parser.add_argument('--plot', action='store_true', help='Save diagnostic plots')
    args = parser.parse_args()

    cfg     = load_case(args.case)
    OUTDIR  = cfg['_output_dir']
    FIGDIR  = cfg['_figures_dir']
    outfile = os.path.join(OUTDIR, 'cfram_result.nc')

    if not os.path.exists(outfile):
        print(f'ERROR: {outfile} not found. Run run_parallel_python.py first.')
        sys.exit(1)

    print(f'Loading: {outfile}')

    # ── Load all dT terms ────────────────────────────────────────────────────
    CANDIDATE = RAD_TERMS + ['atmdyn', 'sfcdyn', 'lhflx', 'shflx']
    dT = {}
    for t in CANDIDATE + ['observed']:
        d = load_nc(outfile, 'dT_' + t)
        if d is not None:
            dT[t] = d

    has_ocndyn = load_nc(outfile, 'dT_ocndyn') is not None
    if 'sfcdyn' in dT and has_ocndyn:
        SUM_TERMS = RAD_TERMS + ['atmdyn', 'sfcdyn']
    else:
        SUM_TERMS = [t for t in CANDIDATE if t in dT]

    dT_sum   = sum(dT[t] for t in SUM_TERMS if t in dT)
    residual = dT_sum - dT['observed']

    # Valid mask: exclude fill values
    valid = (np.abs(dT['observed']) < 900) & (np.abs(dT_sum) < 900)
    res   = np.where(valid, residual, np.nan)          # (nlev+1, nlat, nlon)
    nlev_total = res.shape[0]

    # Pressure labels: file is stored surface→TOA (index 0 = surface? or TOA?)
    # Based on validate_additivity: plev_labels = FORTRAN_PLEV[::-1] + ['sfc']
    # i.e. index 0 = surface pressure, indices 1..N = bottom→top, last = sfc
    # Actually FORTRAN_PLEV is TOA→surface, so reversed it is surface→TOA.
    # Let's just use indices and label the key ones.
    plev_arr = np.concatenate([FORTRAN_PLEV[::-1], [np.nan]])  # surface→TOA + sfc

    # ── Per-level domain-mean |residual| ─────────────────────────────────────
    lev_mean = np.nanmean(np.abs(res), axis=(1, 2))
    lev_max  = np.nanmax( np.abs(res), axis=(1, 2))

    # ── Spatial correlations per level ───────────────────────────────────────
    # Predictors: |dT_k| for each k, plus cross-product |dT_q|*|dT_cloud|
    CORR_TERMS = [t for t in ALL_TERMS if t in dT]
    cross_key  = 'q×cloud'   # WV×Cloud cross-term proxy

    corr_profiles = {t: np.full(nlev_total, np.nan) for t in CORR_TERMS}
    corr_profiles[cross_key] = np.full(nlev_total, np.nan)

    for k in range(nlev_total):
        res_k = np.abs(res[k])                        # 2-D (nlat, nlon)
        if np.all(np.isnan(res_k)):
            continue
        for t in CORR_TERMS:
            pred = np.abs(dT[t][k])
            pred = np.where(valid[k], pred, np.nan)
            corr_profiles[t][k] = pearson_r_spatial(res_k, pred)

        # Cross-product term
        if 'q' in dT and 'cloud' in dT:
            cross = np.abs(dT['q'][k]) * np.abs(dT['cloud'][k])
            cross = np.where(valid[k], cross, np.nan)
            corr_profiles[cross_key][k] = pearson_r_spatial(res_k, cross)

    # ── Summary table ────────────────────────────────────────────────────────
    print(f'\n{"="*65}')
    print(f'Spatial Correlation Analysis: |residual| vs |dT_k|  [{args.case}]')
    print(f'{"="*65}')
    print(f'\nGlobal stats:')
    print(f'  Overall mean|res|  : {np.nanmean(np.abs(res)):.4f} K')
    print(f'  Overall P95 |res|  : {np.nanpercentile(np.abs(res), 95):.4f} K')
    print(f'  Worst level (mean) : k={np.nanargmax(lev_mean)}  '
          f'p={plev_arr[np.nanargmax(lev_mean)]:.0f} hPa  '
          f'mean={np.nanmax(lev_mean):.4f} K')

    # Column-mean correlation for each predictor
    all_keys = CORR_TERMS + [cross_key]
    col_mean_corr = {k: np.nanmean(corr_profiles[k]) for k in all_keys}
    ranked = sorted(col_mean_corr.items(), key=lambda x: -abs(x[1]))

    print(f'\nPredictor ranking (mean |r| across all levels):')
    print(f'  {"Predictor":<14}  {"mean r":<10}')
    print(f'  {"-"*26}')
    for name, r in ranked:
        print(f'  {name:<14}  {r:+.4f}')

    # Per-level top-3 predictors at worst level
    worst_k = int(np.nanargmax(lev_mean))
    worst_p = plev_arr[worst_k]
    level_corr = {k: corr_profiles[k][worst_k] for k in all_keys}
    top3 = sorted(level_corr.items(), key=lambda x: -abs(x[1]))[:3]
    print(f'\nAt worst level (k={worst_k}, ~{worst_p:.0f} hPa):')
    print(f'  mean|res|={lev_mean[worst_k]:.4f} K,  max|res|={lev_max[worst_k]:.4f} K')
    print(f'  Top-3 correlated predictors:')
    for name, r in top3:
        print(f'    {name:<14}  r={r:+.4f}')

    # Cross-term uplift: does q×cloud outperform individual q or cloud?
    if 'q' in col_mean_corr and 'cloud' in col_mean_corr:
        r_q   = col_mean_corr['q']
        r_cld = col_mean_corr['cloud']
        r_x   = col_mean_corr[cross_key]
        print(f'\nWV×Cloud cross-term test (nonlinearity proxy):')
        print(f'  r(|res|, |dT_q|)        = {r_q:+.4f}')
        print(f'  r(|res|, |dT_cloud|)    = {r_cld:+.4f}')
        print(f'  r(|res|, |dT_q×cloud|)  = {r_x:+.4f}')
        if abs(r_x) > max(abs(r_q), abs(r_cld)):
            print('  → Cross-term WINS: WV×Cloud interaction dominates residual')
        else:
            best_ind = 'q' if abs(r_q) > abs(r_cld) else 'cloud'
            print(f'  → Individual term wins: {best_ind}')

    # ── Optional plots ───────────────────────────────────────────────────────
    if args.plot:
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            os.makedirs(FIGDIR, exist_ok=True)

            # Fig 1: Per-level mean|residual| profile
            fig, ax = plt.subplots(figsize=(5, 8))
            ax.plot(lev_mean, range(nlev_total), 'k-', lw=2, label='mean|res|')
            ax.plot(lev_max,  range(nlev_total), 'r--', lw=1, label='max|res|')
            ax.axvline(0.1, color='gray', lw=0.8, linestyle=':', label='0.1K')
            ax.axvline(0.5, color='gray', lw=0.8, linestyle='--', label='0.5K')
            ax.set_xlabel('|residual| (K)')
            ax.set_ylabel('level index (0=surface)')
            ax.set_title(f'Per-level residual  [{args.case}]')
            ax.invert_yaxis()
            ax.legend(fontsize=8)
            plt.tight_layout()
            p = os.path.join(FIGDIR, f'residual_profile_{args.case}.png')
            plt.savefig(p, dpi=120); plt.close()
            print(f'\nSaved: {p}')

            # Fig 2: Correlation profiles for top-5 predictors
            top5_keys = [k for k, _ in ranked[:5]]
            fig, ax = plt.subplots(figsize=(7, 8))
            colors = ['#e41a1c','#377eb8','#4daf4a','#984ea3','#ff7f00']
            for i, key in enumerate(top5_keys):
                ax.plot(corr_profiles[key], range(nlev_total),
                        color=colors[i], lw=1.5, label=key)
            ax.axvline(0, color='k', lw=0.5)
            ax.set_xlabel('r(|residual|, |predictor|)')
            ax.set_ylabel('level index (0=surface)')
            ax.set_title(f'Spatial correlation profiles  [{args.case}]')
            ax.invert_yaxis()
            ax.legend(fontsize=8)
            plt.tight_layout()
            p = os.path.join(FIGDIR, f'residual_corr_profiles_{args.case}.png')
            plt.savefig(p, dpi=120); plt.close()
            print(f'Saved: {p}')

            # Fig 3: Spatial map of max|residual| across all levels
            res_maxmap = np.nanmax(np.abs(res), axis=0)
            fig, ax = plt.subplots(figsize=(9, 5))
            im = ax.imshow(res_maxmap, origin='lower', aspect='auto',
                           cmap='YlOrRd', vmin=0, vmax=5)
            plt.colorbar(im, ax=ax, label='K')
            ax.set_title(f'max_lev |residual| (K)  [{args.case}]')
            ax.set_xlabel('lon index'); ax.set_ylabel('lat index')
            plt.tight_layout()
            p = os.path.join(FIGDIR, f'residual_maxmap_{args.case}.png')
            plt.savefig(p, dpi=120); plt.close()
            print(f'Saved: {p}')

        except ImportError:
            print('matplotlib not available; skipping plots.')

    print('\nDone.')


if __name__ == '__main__':
    main()
