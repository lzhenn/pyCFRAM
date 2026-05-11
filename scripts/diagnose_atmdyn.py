#!/usr/bin/env python3
"""Diagnose atmdyn at heating-rate level (Cai-Tung 2012 / colleague's framework).

Two-layer CFRAM check:
  Layer 1 — heating rate (W/m^2 per layer):
      Q_atmdyn[atm] = -frc_warm[atm]
      → if magnitude is unreasonable (e.g. tropical mean |Q|>10 W/m^2/layer),
        the bug is at *physical input* level (missing column heating, etc.)

  Layer 2 — Planck inversion (K):
      dT_atmdyn = -(∂R/∂T)^-1 · Q_atmdyn
      → if heating rate is reasonable but dT saturates, the bug is Planck
        conditioning (stratospheric ill-posedness).

Usage
-----
    python3 scripts/diagnose_atmdyn.py <case_name>

Outputs three figures into cases/<case>/figures/:
  diag_zonalmean_frc_warm.png   — frc_warm and frc_atmdyn (=-frc_warm) at atm
  diag_zonalmean_dT_radiative.png — sum of radiative dT terms vs observed
  diag_zonalmean_dT_atmdyn.png  — dT_atmdyn (saturated panel) for context
"""
import os, sys
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case


def zonal_mean(arr3d):
    """Zonal mean of (lev, lat, lon) -> (lev, lat). Drop fills."""
    arr = np.where(np.abs(arr3d) > 900, np.nan, arr3d)
    return np.nanmean(arr, axis=2)


def panel(ax, x, y, field, title, levels, cmap, units, hline_at=None):
    cf = ax.contourf(x, y, field, levels=levels, cmap=cmap,
                     norm=mcolors.BoundaryNorm(levels, ncolors=cmap.N, extend='both'),
                     extend='both')
    cb = plt.colorbar(cf, ax=ax, label=units, ticks=levels)
    ax.set_xlabel('Latitude')
    ax.set_ylabel('Pressure (hPa)')
    # Match OLD CFRAM y-axis: 1000 hPa (bottom) to ~1 hPa (top), log scale to expose strat.
    ax.set_yscale('log')
    ax.set_ylim(1000, 1)
    ax.set_yticks([1000, 800, 600, 500, 400, 300, 200, 100, 50, 20, 10, 5, 2, 1])
    ax.get_yaxis().set_major_formatter(plt.matplotlib.ticker.ScalarFormatter())
    ax.set_title(title, fontsize=11)
    if hline_at is not None:
        ax.axhline(hline_at, color='k', linewidth=0.4, linestyle='--', alpha=0.5)


def main(case_name):
    cfg = load_case(case_name)
    nc_file = os.path.join(cfg['_output_dir'], 'cfram_result.nc')
    if not os.path.exists(nc_file):
        sys.exit('Missing: ' + nc_file)

    nc = Dataset(nc_file)
    lats = np.array(nc.variables['lat'][:])
    lev_full = np.array(nc.variables['lev'][:])
    nlev = len(lev_full)
    print('lev range: %.1f -> %.1f hPa, n=%d' % (lev_full.min(), lev_full.max(), nlev))

    # Helper: load and zonal-average
    def load_zm(varname):
        if varname not in nc.variables:
            print('  missing var:', varname)
            return None
        return zonal_mean(np.array(nc.variables[varname][:], dtype=np.float64))

    fig_dir = cfg['_figures_dir']
    os.makedirs(fig_dir, exist_ok=True)

    # ====== Figure 1: frc_warm + frc_atmdyn (heating-rate level) ======
    frc_warm_zm = load_zm('frc_warm')
    if frc_warm_zm is None:
        sys.exit('frc_warm not in output — solver may be older version')

    # frc_atmdyn[atm] = -frc_warm[atm];  surface row = 0 by construction
    frc_atmdyn_zm = -frc_warm_zm.copy()
    sfc_idx = int(np.argmax(lev_full))   # surface = highest pressure
    frc_atmdyn_zm[sfc_idx, :] = 0.0       # surface row zeroed in algorithm

    levels_w = [-50, -20, -10, -5, -2, -1, 0, 1, 2, 5, 10, 20, 50]
    cmap = plt.cm.RdBu_r

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    panel(axes[0], lats, lev_full, frc_warm_zm,
          'frc_warm  (= ΔR all-non-T-perturbed)',
          levels_w, cmap, 'W/m² per layer', hline_at=200)
    panel(axes[1], lats, lev_full, frc_atmdyn_zm,
          'frc_atmdyn  (= −frc_warm at atm)',
          levels_w, cmap, 'W/m² per layer', hline_at=200)

    fig.suptitle('Heating-rate level diagnostic — %s' % case_name, fontsize=13)
    plt.tight_layout()
    out1 = os.path.join(fig_dir, 'diag_zonalmean_frc_warm.png')
    plt.savefig(out1, dpi=130, bbox_inches='tight')
    print('Saved:', out1)
    plt.close()

    # Quick magnitude print: tropical band 20S-20N atm mean
    trop = (lats >= -20) & (lats <= 20)
    atm_mask = lev_full > 100      # troposphere only (above 100 hPa = stratosphere excluded)
    atm_mask[sfc_idx] = False
    trop_mean_q = np.nanmean(np.abs(frc_atmdyn_zm[atm_mask][:, trop]))
    strat_mean_q = np.nanmean(np.abs(frc_atmdyn_zm[~atm_mask & (lev_full < 100)][:, trop]))
    print()
    print('===== HEATING-RATE LEVEL DIAGNOSTIC (tropical 20S-20N) =====')
    print('Mean |frc_atmdyn|  in TROPOSPHERE (trop+lev>100hPa):  %.2f W/m²/layer' % trop_mean_q)
    print('Mean |frc_atmdyn|  in STRATOSPHERE (lev<100hPa):       %.2f W/m²/layer' % strat_mean_q)
    print()
    print('Reference scale: real atm dynamics typically  2-5 W/m²/layer')
    print('  > 10 W/m²/layer in trop -> missing column-heating physics likely')
    print('  > 30 W/m²/layer in strat -> stratospheric Planck issue likely')

    # ====== Figure 2: radiative dT additivity check ======
    rad_terms = ['co2', 'q', 'cloud', 'aerosol', 'o3', 'solar', 'albedo', 'ts']
    rad_sum = np.zeros_like(frc_warm_zm)
    rad_each = {}
    for t in rad_terms:
        zm = load_zm('dT_' + t)
        if zm is not None:
            rad_each[t] = zm
            rad_sum += np.nan_to_num(zm)

    dT_obs_zm = load_zm('dT_observed')
    dT_warm_zm = load_zm('dT_warm')

    # Match OLD CFRAM colorbar exactly: [-15, -7, -2, -0.5, 0, 0.5, 2, 7, 15]
    levels_t = [-15, -7, -2, -0.5, 0, 0.5, 2, 7, 15]
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    panel(axes[0, 0], lats, lev_full, dT_obs_zm,
          'dT_observed (T_warm − T_base)', levels_t, cmap, 'K', hline_at=200)
    panel(axes[0, 1], lats, lev_full, dT_warm_zm,
          'dT_warm  (Planck-inverse of frc_warm)', levels_t, cmap, 'K', hline_at=200)
    panel(axes[1, 0], lats, lev_full, rad_sum,
          'Σ dT_radiative  (8 partial perturbations)', levels_t, cmap, 'K', hline_at=200)
    panel(axes[1, 1], lats, lev_full, dT_obs_zm - rad_sum,
          'dT_observed − Σ dT_radiative  (= dynamics + nonrad)',
          levels_t, cmap, 'K', hline_at=200)
    fig.suptitle('Radiative additivity — %s' % case_name, fontsize=13)
    plt.tight_layout()
    out2 = os.path.join(fig_dir, 'diag_zonalmean_dT_radiative.png')
    plt.savefig(out2, dpi=130, bbox_inches='tight')
    print('Saved:', out2)
    plt.close()

    # ====== Figure 3: dT_atmdyn (the saturated one) for reference ======
    dT_atmdyn_zm = load_zm('dT_atmdyn')
    dT_lhflx_zm = load_zm('dT_lhflx')
    dT_shflx_zm = load_zm('dT_shflx')
    dT_sfcdyn_zm = load_zm('dT_sfcdyn')
    dT_ocndyn_zm = load_zm('dT_ocndyn')

    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    panel(axes[0, 0], lats, lev_full, dT_atmdyn_zm,
          'dT_atmdyn (the saturated panel)', levels_t, cmap, 'K', hline_at=200)
    panel(axes[0, 1], lats, lev_full, dT_lhflx_zm,
          'dT_lhflx', levels_t, cmap, 'K', hline_at=200)
    panel(axes[0, 2], lats, lev_full, dT_shflx_zm,
          'dT_shflx', levels_t, cmap, 'K', hline_at=200)
    panel(axes[1, 0], lats, lev_full, dT_sfcdyn_zm,
          'dT_sfcdyn', levels_t, cmap, 'K', hline_at=200)
    panel(axes[1, 1], lats, lev_full, dT_ocndyn_zm,
          'dT_ocndyn', levels_t, cmap, 'K', hline_at=200)
    panel(axes[1, 2], lats, lev_full,
          dT_atmdyn_zm + np.nan_to_num(dT_lhflx_zm) + np.nan_to_num(dT_shflx_zm)
          + np.nan_to_num(dT_sfcdyn_zm) + np.nan_to_num(dT_ocndyn_zm),
          'Σ dT_(atmdyn+lh+sh+sfcdyn+ocndyn)', levels_t, cmap, 'K', hline_at=200)
    fig.suptitle('Non-radiative dT components — %s' % case_name, fontsize=13)
    plt.tight_layout()
    out3 = os.path.join(fig_dir, 'diag_zonalmean_dT_atmdyn.png')
    plt.savefig(out3, dpi=130, bbox_inches='tight')
    print('Saved:', out3)
    plt.close()

    nc.close()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit('usage: diagnose_atmdyn.py <case_name>')
    main(sys.argv[1])
