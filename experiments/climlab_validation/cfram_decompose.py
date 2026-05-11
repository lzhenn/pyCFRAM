#!/usr/bin/env python3
"""CFRAM decomposition driver: climlab RCE → Python CFRAM → validate vs obs.

Uses cfram_core for all CFRAM math. This script is just glue: it loads the
two RCE NetCDFs, builds a closure that calls climlab's RRTMG offline, and
hands the closure to cfram_core. Then it plots and prints validation stats.
"""
import os
import sys
import numpy as np
from netCDF4 import Dataset
import climlab
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from cfram_core import (R_per_layer, S_per_layer, compute_drdt, solve_dT,
                        lump_atm_residual)

OUTDIR = os.path.join(HERE, 'output')


def load_rce(label):
    """Load equilibrium state saved by run_rce.py."""
    nc = Dataset(os.path.join(OUTDIR, 'rce_%s.nc' % label))
    s = {
        'lev': np.array(nc.variables['lev'][:]),
        'Tatm': np.array(nc.variables['Tatm'][:]),
        'q': np.array(nc.variables['q'][:]),
        'Ts': float(nc.variables['Ts'][...]),
        'co2': float(nc.co2_ppm) * 1e-6,
    }
    nc.close()
    return s


def build_evaluator(s_base):
    """Build an offline RRTMG evaluator pinned to base column structure.

    Returns
    -------
    eval_RS : callable f(Tatm, Ts, q=None, co2=None) -> (R, S)
    eval_R  : callable f(Tatm, Ts) -> R only (for compute_drdt)
    """
    cs = climlab.column_state(num_lev=len(s_base['lev']), water_depth=1.0)
    cs['Tatm'][:] = s_base['Tatm']
    cs['Ts'][:] = s_base['Ts']
    h2o = climlab.radiation.ManabeWaterVapor(state=cs, name='h2o')
    rad = climlab.radiation.RRTMG(
        state=cs, specific_humidity=s_base['q'].copy(),
        albedo=0.3, insolation=341.5, icld=0, name='rad',
    )
    rad.absorber_vmr['CO2'] = s_base['co2']
    # store baseline q so caller can reset
    base_q = s_base['q'].copy()
    base_co2 = s_base['co2']

    def eval_RS(Tatm, Ts, q=None, co2=None):
        rad.Tatm[:] = Tatm
        rad.Ts[:] = Ts
        rad.specific_humidity[:] = q if q is not None else base_q
        rad.absorber_vmr['CO2'] = co2 if co2 is not None else base_co2
        rad.compute_diagnostics()
        R = R_per_layer(np.array(rad.LW_flux_up).flatten(),
                        np.array(rad.LW_flux_down).flatten())
        S = S_per_layer(np.array(rad.SW_flux_up).flatten(),
                        np.array(rad.SW_flux_down).flatten())
        return R, S

    def eval_R(Tatm, Ts):
        return eval_RS(Tatm, Ts)[0]

    return eval_RS, eval_R


def decompose_2xco2(s1, s2):
    """Lu & Cai 2009 Eq.4 decomposition (no atmdyn lumping).

        ΔT = drdt⁻¹ · [ΔF_ext + Δ(S−R)_w + ΔQ_dry]

    where
        ΔF_ext  = -ΔR_co2|_T1            (linearized at base, Eq.5)
        Δ(S−R)_w = ΔS_q|_T1 − ΔR_q|_T1   (linearized at base, Eq.6)
        ΔQ_dry  = (R−S)_2 − (R−S)_1       (full state difference, no linearization)

    For our climlab setup with no explicit LH/SH processes, ΔQ_lh = ΔQ_sh = 0.
    Q_dry per layer at equilibrium = R - S (Lu/Cai Eq.2 with Q_lh = Q_sh = 0).
    """
    nlev = len(s1['Tatm'])
    n = nlev + 1

    eval_RS, eval_R = build_evaluator(s1)

    # ===== Planck matrix at base T1 =====
    print('Computing ∂R/∂T (%dx%d perturbations)...' % (n, n))
    drdt, R0 = compute_drdt(eval_R, s1['Tatm'], s1['Ts'], dT=1.0)
    drdt_inv = np.linalg.inv(drdt)
    print('  drdt diag: top=%.3f mid=%.3f sfc=%.3f W/m²/K  (cond=%.1e)' %
          (drdt[0, 0], drdt[nlev//2, nlev//2], drdt[nlev, nlev],
           np.linalg.cond(drdt)))

    # ===== Reference (R, S) at base state 1 =====
    R_base, S_base = eval_RS(s1['Tatm'], s1['Ts'], q=s1['q'], co2=s1['co2'])

    # ===== Linearized rad perturbations (T held at T1) =====
    R_co2,  S_co2  = eval_RS(s1['Tatm'], s1['Ts'], q=s1['q'], co2=s2['co2'])
    R_q,    S_q    = eval_RS(s1['Tatm'], s1['Ts'], q=s2['q'], co2=s1['co2'])
    R_warm, S_warm = eval_RS(s1['Tatm'], s1['Ts'], q=s2['q'], co2=s2['co2'])

    frc = {  # ΔR per layer, T held at T1
        'co2':  R_co2  - R_base,
        'q':    R_q    - R_base,
        'warm': R_warm - R_base,
    }
    dS = {
        'co2':  S_co2  - S_base,
        'q':    S_q    - S_base,
        'warm': S_warm - S_base,
    }

    # ===== ΔQ_dry from FULL state energy balance =====
    # Q_dry per layer = R - S (Lu/Cai Eq.2 at equilibrium with Q_lh=Q_sh=0).
    # Take difference between FULL state 2 and FULL state 1 (NOT at base T).
    R_full2, S_full2 = eval_RS(s2['Tatm'], s2['Ts'], q=s2['q'], co2=s2['co2'])
    Q_dry_1 = R_base   - S_base
    Q_dry_2 = R_full2 - S_full2
    dQ_dry = Q_dry_2 - Q_dry_1

    # ===== Solve partial dT_X = drdt⁻¹ · forcing_X =====
    # Lu/Cai Eq.4 partial form: ΔT_X = drdt⁻¹ · ΔF_X (positive ΔF = warming).
    # Our solve_dT(frc, delta_S) returns drdt⁻¹ · (delta_S - frc), so:
    #   - For rad: pass frc=ΔR, delta_S=ΔS  →  drdt⁻¹·(ΔS-ΔR) = drdt⁻¹·ΔF_LuCai ✓
    #   - For ΔQ_dry: pass frc=0, delta_S=ΔQ_dry  →  drdt⁻¹·ΔQ_dry ✓
    dT = {
        'ext':   solve_dT(drdt_inv, frc['co2'],  dS['co2']),   # ΔT_ext
        'w':     solve_dT(drdt_inv, frc['q'],    dS['q']),     # ΔT_w
        'dry':   solve_dT(drdt_inv, np.zeros_like(dQ_dry), dQ_dry),  # ΔT_dry
        'warm':  solve_dT(drdt_inv, frc['warm'], dS['warm']),  # ΔT_co2 + ΔT_q (combined rad)
    }

    # Lu/Cai 5-term sum (here only 3 nonzero since Q_lh = Q_sh = 0 in our setup):
    # ΔT_obs ≈ ΔT_ext + ΔT_w + ΔT_dry (+ ΔT_lh + ΔT_sh)
    dT_obs = np.concatenate([s2['Tatm'] - s1['Tatm'], [s2['Ts'] - s1['Ts']]])
    dT_sum = dT['ext'] + dT['w'] + dT['dry']

    # Also compute the OLD pyCFRAM-style atmdyn for comparison (lump shortcut)
    frc_atm, dS_atm = lump_atm_residual(frc['warm'], dS['warm'])
    dT['atmdyn_lump'] = solve_dT(drdt_inv, frc_atm, dS_atm)
    dT_sum_lump = dT['ext'] + dT['w'] + dT['atmdyn_lump']

    return {
        'lev': np.concatenate([s1['lev'], [1013.0]]),
        'frc': frc, 'dS': dS, 'dT': dT,
        'dQ_dry': dQ_dry, 'Q_dry_1': Q_dry_1, 'Q_dry_2': Q_dry_2,
        'dT_obs': dT_obs, 'dT_sum_lucai': dT_sum, 'dT_sum_lump': dT_sum_lump,
        'drdt': drdt, 'R_base': R_base, 'S_base': S_base,
    }


def report(res):
    obs = res['dT_obs']
    sum_lc = res['dT_sum_lucai']
    sum_lp = res['dT_sum_lump']

    print('\n===== Surface dT (Lu/Cai Eq.4 decomposition) =====')
    for k in ['ext', 'w', 'dry', 'warm']:
        print('  ΔT_%-7s = %+.3f K' % (k, res['dT'][k][-1]))
    print('  Σ Lu/Cai      = %+.3f K  (ext + w + dry)' % sum_lc[-1])
    print('  Σ pyCFRAM lump = %+.3f K  (ext + w + atmdyn_lump)' % sum_lp[-1])
    print('  observed       = %+.3f K' % obs[-1])
    resid_lc = obs[-1] - sum_lc[-1]
    resid_lp = obs[-1] - sum_lp[-1]
    print('  residual Lu/Cai  = %+.3f K  (%.1f%% of obs)' %
          (resid_lc, 100. * resid_lc / obs[-1] if obs[-1] != 0 else 0))
    print('  residual lumping = %+.3f K  (%.1f%% of obs)' %
          (resid_lp, 100. * resid_lp / obs[-1] if obs[-1] != 0 else 0))

    print('\nProfile rmsd (Σ partials vs observed):')
    print('  Lu/Cai:  rmsd=%.3f K, max|diff|=%.3f K' %
          (np.sqrt(np.mean((sum_lc - obs)**2)), np.max(np.abs(sum_lc - obs))))
    print('  Lumping: rmsd=%.3f K, max|diff|=%.3f K' %
          (np.sqrt(np.mean((sum_lp - obs)**2)), np.max(np.abs(sum_lp - obs))))

    # ΔQ_dry sanity: should be 0 in stratosphere (no convection up there)
    nlev = len(res['lev']) - 1
    strat_mask = res['lev'][:-1] < 100  # atm rows only
    if strat_mask.any():
        print('\nΔQ_dry sanity (W/m² per layer):')
        print('  stratosphere (lev<100): max|ΔQ_dry| = %.3f W/m²' %
              np.max(np.abs(res['dQ_dry'][:nlev][strat_mask])))
        print('  troposphere  (lev>=100): max|ΔQ_dry| = %.3f W/m²' %
              np.max(np.abs(res['dQ_dry'][:nlev][~strat_mask])))


def save_nc(res, path):
    nc = Dataset(path, 'w')
    n = len(res['lev'])
    nc.createDimension('lev', n)
    nc.createVariable('lev', 'f8', ('lev',))[:] = res['lev']
    for name, arr in [
        ('dT_obs', res['dT_obs']),
        ('dT_sum_lucai', res['dT_sum_lucai']),
        ('dT_sum_lump',  res['dT_sum_lump']),
        ('dT_ext',    res['dT']['ext']),
        ('dT_w',      res['dT']['w']),
        ('dT_dry',    res['dT']['dry']),
        ('dT_warm',   res['dT']['warm']),
        ('dT_atmdyn_lump', res['dT']['atmdyn_lump']),
        ('frc_co2', res['frc']['co2']),  ('frc_q', res['frc']['q']),
        ('frc_warm', res['frc']['warm']),
        ('dS_co2',  res['dS']['co2']),   ('dS_q', res['dS']['q']),
        ('dS_warm', res['dS']['warm']),
        ('dQ_dry',  res['dQ_dry']),
        ('Q_dry_1', res['Q_dry_1']), ('Q_dry_2', res['Q_dry_2']),
        ('drdt_diag', np.diag(res['drdt'])),
        ('R_base', res['R_base']), ('S_base', res['S_base']),
    ]:
        nc.createVariable(name, 'f8', ('lev',))[:] = arr
    nc.close()
    print('Saved:', path)


def plot_fig(res, path):
    plev = res['lev']
    fig, axes = plt.subplots(1, 5, figsize=(20, 6), sharey=True)

    panels = [
        ('Total ΔT', [
            (res['dT_obs'],         'climlab obs',     'black',  '-',  'o'),
            (res['dT_sum_lucai'],   'CFRAM Lu/Cai Σ',  'red',    '--', 's'),
            (res['dT_sum_lump'],    'pyCFRAM lump Σ',  'orange', ':',  '^'),
        ]),
        ('ΔT_ext  (CO2)',  [(res['dT']['ext'],   None, 'red',    '-', 'o')]),
        ('ΔT_w  (water vapor)', [(res['dT']['w'],     None, 'blue',   '-', 'o')]),
        ('ΔT_dry  (Lu/Cai ΔQ_dry)', [(res['dT']['dry'],   None, 'purple', '-', 'o')]),
        ('ΔT_atmdyn (pyCFRAM lump)', [(res['dT']['atmdyn_lump'], None, 'orange', '-', 'o')]),
    ]
    for ax, (title, lines) in zip(axes, panels):
        for arr, lbl, c, ls, m in lines:
            ax.plot(arr, plev, color=c, linestyle=ls, marker=m,
                    markersize=3, label=lbl)
        ax.set_yscale('log')
        ax.set_ylim(1013, plev[0])
        ax.set_yticks([1000, 500, 200, 100, 50, 20, 10])
        ax.get_yaxis().set_major_formatter(plt.matplotlib.ticker.ScalarFormatter())
        ax.axvline(0, color='k', linewidth=0.4, alpha=0.5)
        ax.axhline(100, color='k', linewidth=0.3, linestyle=':', alpha=0.4)
        ax.set_xlabel('ΔT (K)')
        ax.set_title(title, fontsize=10)
        if any(l[1] for l in lines):
            ax.legend(fontsize=8, loc='upper right')
    axes[0].set_ylabel('Pressure (hPa)')
    fig.suptitle('Lu/Cai vs pyCFRAM-lump CFRAM decomposition on climlab RCE (2×CO2)',
                 fontsize=12)
    plt.tight_layout()
    plt.savefig(path, dpi=130, bbox_inches='tight')
    print('Saved:', path)


def main():
    s1 = load_rce('1xco2')
    s2 = load_rce('2xco2')
    print('Base 1×CO2: Ts=%.3f K  →  Pert 2×CO2: Ts=%.3f K  (ΔTs=%+.3f K)' %
          (s1['Ts'], s2['Ts'], s2['Ts'] - s1['Ts']))

    res = decompose_2xco2(s1, s2)
    report(res)
    save_nc(res, os.path.join(OUTDIR, 'cfram_partials.nc'))
    plot_fig(res, os.path.join(OUTDIR, 'cfram_partials.png'))


if __name__ == '__main__':
    main()
