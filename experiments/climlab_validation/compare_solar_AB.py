#!/usr/bin/env python3
"""Compare path A (pure-python CFRAM) vs path B (production pyCFRAM Fortran)
vs climlab observed RCE on the +1% TOA solar experiment.

Inputs:
- experiments/climlab_validation/output/cfram_partials_solar.nc  (path A)
- cases/climlab_solar/output/cfram_result.nc                      (path B)
- experiments/climlab_validation/output/rce_solar_{base,warm}.nc  (obs)
"""
import os
import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))  # pyCFRAM/
OUTDIR = os.path.join(HERE, 'output')

PLEV_TOA2SFC = np.array([
    1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150,
    175, 200, 225, 250, 300, 350, 400, 450, 500, 550,
    600, 650, 700, 750, 775, 800, 825, 850, 875, 900,
    925, 950, 975, 1000
], dtype=float)
LEV_PLOT = np.concatenate([PLEV_TOA2SFC, [1013.0]])  # add ground


def load_path_A():
    """Path A: pure-python decomposition. Lev order TOA→sfc + 1013 ground."""
    nc = Dataset(os.path.join(OUTDIR, 'cfram_partials_solar.nc'))
    out = {}
    for v in ('dT_obs', 'dT_sum_lucai', 'dT_solar', 'dT_q', 'dT_co2',
             'dT_dry', 'dT_warm', 'dT_atmdyn_lump',
             'frc_solar', 'frc_q', 'dS_solar', 'dS_q'):
        out[v] = np.array(nc.variables[v][:])
    nc.close()
    return out


def load_path_B():
    """Path B: production pyCFRAM. NetCDF stores lev sfc→TOA + ground (1013).

    Reorder atm rows to TOA→sfc + ground (matching path A's convention).
    """
    nc = Dataset(os.path.join(ROOT, 'cases/climlab_solar/output/cfram_result.nc'))
    NLEV = nc.dimensions['lev'].size - 1
    out = {}
    for v in ('dT_observed', 'dT_solar', 'dT_q', 'dT_co2', 'dT_dry',
             'dT_warm', 'dT_atmdyn', 'dT_sfcdyn', 'dT_ocndyn',
             'frc_solar', 'frc_q', 'frc_warm', 'frc_full'):
        if v in nc.variables:
            arr = np.array(nc.variables[v][:, 0, 0])  # sfc→TOA atm + ground
            atm_rev = arr[:NLEV][::-1]                 # → TOA→sfc atm
            out[v.replace('observed', 'obs')] = np.concatenate([atm_rev, arr[NLEV:]])
    nc.close()
    return out


def load_obs():
    """Climlab native ΔT (TOA→sfc + Ts at end)."""
    base = Dataset(os.path.join(OUTDIR, 'rce_solar_base.nc'))
    warm = Dataset(os.path.join(OUTDIR, 'rce_solar_warm.nc'))
    dT_atm = np.array(warm.variables['Tatm'][:]) - np.array(base.variables['Tatm'][:])
    dTs = float(warm.variables['Ts'][...]) - float(base.variables['Ts'][...])
    base.close(); warm.close()
    return np.concatenate([dT_atm, [dTs]])


def report(A, B, obs):
    print('=' * 60)
    print(' Path A (pure-python) vs Path B (production pyCFRAM Fortran)')
    print('  +1% TOA solar perturbation, climlab RCE 1×1 column')
    print('=' * 60)

    print('\n--- Surface values (lev[-1]) ---')
    print('  %-14s | %-9s | %-9s | %-9s' %
          ('var', 'PathA', 'PathB', 'A-B'))
    print('  ' + '-' * 50)
    rows = [
        ('ΔT_obs', A['dT_obs'][-1], B['dT_obs'][-1]),
        ('ΔT_solar',  A['dT_solar'][-1], B['dT_solar'][-1]),
        ('ΔT_q',      A['dT_q'][-1],     B['dT_q'][-1]),
        ('ΔT_co2',    A['dT_co2'][-1],   B['dT_co2'][-1]),
        ('ΔT_dry',    A['dT_dry'][-1],   B['dT_dry'][-1]),
        ('ΔT_warm',   A['dT_warm'][-1],  B['dT_warm'][-1]),
    ]
    for name, av, bv in rows:
        print('  %-14s | %+8.4f | %+8.4f | %+8.4f' % (name, av, bv, av - bv))

    print('\n--- Closure (Σ dT_X vs ΔT_obs) at surface ---')
    sum_A = A['dT_solar'][-1] + A['dT_q'][-1] + A['dT_dry'][-1]
    sum_B = B['dT_solar'][-1] + B['dT_q'][-1] + B['dT_dry'][-1]
    print('  Path A:  Σ = %+.4f K, obs = %+.4f K, residual = %+.4f K (%.2f%%)' %
          (sum_A, obs[-1], obs[-1] - sum_A,
           100. * (obs[-1] - sum_A) / obs[-1]))
    print('  Path B:  Σ = %+.4f K, obs = %+.4f K, residual = %+.4f K (%.2f%%)' %
          (sum_B, obs[-1], obs[-1] - sum_B,
           100. * (obs[-1] - sum_B) / obs[-1]))
    print('  A vs B:  diff in sum = %+.4f K  (target <0.001 K)' %
          (sum_A - sum_B))

    print('\n--- Profile rmsd (atm only, exclude ground) ---')
    n_atm = len(A['dT_obs']) - 1   # exclude ground row
    obs_atm = obs[:n_atm]
    sum_A_atm = (A['dT_solar'] + A['dT_q'] + A['dT_dry'])[:n_atm]
    sum_B_atm = (B['dT_solar'] + B['dT_q'] + B['dT_dry'])[:n_atm]
    print('  Path A:  rmsd vs obs = %.4f K, max|diff| = %.4f K' %
          (np.sqrt(np.mean((sum_A_atm - obs_atm)**2)),
           np.max(np.abs(sum_A_atm - obs_atm))))
    print('  Path B:  rmsd vs obs = %.4f K, max|diff| = %.4f K' %
          (np.sqrt(np.mean((sum_B_atm - obs_atm)**2)),
           np.max(np.abs(sum_B_atm - obs_atm))))

    print('\n--- Forcing comparison (W/m² at surface) ---')
    print('  frc_solar:    A=%+8.4f, B=%+8.4f  (Path A is LW-only;' %
          (A['frc_solar'][-1], B['frc_solar'][-1]))
    print('                                         Path B is LW+SW net)')
    print('  dS_solar (A): %+8.4f' % A['dS_solar'][-1])
    print('  -> A puts SW absorption in dS, B puts it in frc; both yield same total')

    print('\n--- Validation criteria ---')
    crits = [
        ('A closure < 1%',
         abs(100. * (obs[-1] - sum_A) / obs[-1]) < 1.0),
        ('B closure < 1%',
         abs(100. * (obs[-1] - sum_B) / obs[-1]) < 1.0),
        ('A-B sum identical (<1e-3 K)',
         abs(sum_A - sum_B) < 1e-3),
        ('A solar dominates (>50%)',
         100. * A['dT_solar'][-1] / obs[-1] > 50),
        ('B solar dominates (>50%)',
         100. * B['dT_solar'][-1] / obs[-1] > 50),
        ('A & B dT_q match (<1e-3 K)',
         abs(A['dT_q'][-1] - B['dT_q'][-1]) < 1e-3),
        ('A & B dT_co2 ≈ 0 (<5e-3)',
         abs(A['dT_co2'][-1]) < 5e-3 and abs(B['dT_co2'][-1]) < 5e-3),
    ]
    for name, ok in crits:
        print('  [%s] %s' % ('PASS' if ok else 'FAIL', name))


def plot_compare(A, B, obs, path):
    plev = LEV_PLOT
    fig, axes = plt.subplots(1, 5, figsize=(22, 7), sharey=True)

    sum_A = A['dT_solar'] + A['dT_q'] + A['dT_dry']
    sum_B = B['dT_solar'] + B['dT_q'] + B['dT_dry']

    panels = [
        ('Total ΔT (closure)', [
            (obs,     'climlab obs',    'black',  '-',  'o', 1.0),
            (sum_A,   'Path A Σ(LuCai)','red',    '--', 's', 0.7),
            (sum_B,   'Path B Σ(LuCai)','blue',   ':',  '^', 0.7),
        ]),
        ('ΔT_solar', [
            (A['dT_solar'], 'A', 'red',  '-', 's', 0.7),
            (B['dT_solar'], 'B', 'blue', ':', '^', 0.7),
        ]),
        ('ΔT_q (water vapor)', [
            (A['dT_q'], 'A', 'red',  '-', 's', 0.7),
            (B['dT_q'], 'B', 'blue', ':', '^', 0.7),
        ]),
        ('ΔT_co2 (sanity ≈0)', [
            (A['dT_co2'], 'A', 'red',  '-', 's', 0.7),
            (B['dT_co2'], 'B', 'blue', ':', '^', 0.7),
        ]),
        ('ΔT_dry (Lu/Cai ΔQ_dry)', [
            (A['dT_dry'], 'A', 'red',  '-', 's', 0.7),
            (B['dT_dry'], 'B', 'blue', ':', '^', 0.7),
        ]),
    ]
    for ax, (title, lines) in zip(axes, panels):
        for arr, lbl, c, ls, m, alpha in lines:
            ax.plot(arr, plev, color=c, linestyle=ls, marker=m,
                    markersize=4, label=lbl, alpha=alpha)
        ax.set_yscale('log')
        ax.set_ylim(1013, plev[0])
        ax.set_yticks([1000, 500, 200, 100, 50, 20, 10])
        ax.get_yaxis().set_major_formatter(plt.matplotlib.ticker.ScalarFormatter())
        ax.axvline(0, color='k', linewidth=0.4, alpha=0.5)
        ax.axhline(100, color='k', linewidth=0.3, linestyle=':', alpha=0.4)
        ax.set_xlabel('ΔT (K)')
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=8, loc='upper right')
    axes[0].set_ylabel('Pressure (hPa)')
    fig.suptitle('Path A (pure-python) vs Path B (production pyCFRAM Fortran)\n'
                 '+1%% TOA solar climlab RCE  |  ΔTs_obs = %+.3f K' % obs[-1],
                 fontsize=12)
    plt.tight_layout()
    plt.savefig(path, dpi=130, bbox_inches='tight')
    print('Saved:', path)


def main():
    A = load_path_A()
    B = load_path_B()
    obs = load_obs()
    report(A, B, obs)
    plot_compare(A, B, obs, os.path.join(OUTDIR, 'compare_solar_AB.png'))


if __name__ == '__main__':
    main()
