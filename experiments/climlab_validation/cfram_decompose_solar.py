#!/usr/bin/env python3
"""Path A: pure-python CFRAM decomposition of climlab +1% solar RCE pair.

Setup
-----
Base:  insolation = 341.5 W/m², CO2 = 348 ppm
Warm:  insolation = 341.5 × 1.01 = 345.015 W/m² (+3.515 W/m² globally), same CO2
No clouds, no aerosols.

Decomposition (Lu/Cai 2009 Eq.4, no LH/SH since RCE has no LH/SH)
-----------------------------------------------------------------
ΔT_obs ≈ ΔT_solar + ΔT_w + ΔT_dry

where each partial solves drdt · ΔT_X = ΔS_X − ΔR_X (T held at base):
- ΔT_solar = drdt⁻¹ · (ΔS_solar − ΔR_solar)        # ΔR_solar ≈ 0 (LW unchanged)
- ΔT_w     = drdt⁻¹ · (ΔS_q − ΔR_q)
- ΔT_dry   = drdt⁻¹ · ΔQ_dry                       # ΔQ_dry from full-state diff

Validation criteria
-------------------
1. Surface closure: |ΔT_obs − Σ ΔT_X|_sfc < 1% of ΔT_obs
2. ΔT_solar dominates (since solar is the only external forcing)
3. ΔT_co2/o3/cld/aerosol ≈ 0 (none of these changed)
4. ΔT_w is positive (Manabe RH-fixed: q rises with T → positive H2O feedback)
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


def load_rce_solar(label):
    """Load equilibrium state saved by run_rce_solar.py."""
    nc = Dataset(os.path.join(OUTDIR, 'rce_solar_%s.nc' % label))
    s = {
        'lev':         np.array(nc.variables['lev'][:]),
        'Tatm':        np.array(nc.variables['Tatm'][:]),
        'q':           np.array(nc.variables['q'][:]),
        'o3':          np.array(nc.variables['o3'][:]),
        'Ts':          float(nc.variables['Ts'][...]),
        'co2':         float(nc.co2_ppm) * 1e-6,
        'insolation':  float(nc.insolation_Wm2),
        'albedo':      float(nc.albedo),
    }
    nc.close()
    return s


def build_evaluator(s_base):
    """Offline RRTMG evaluator pinned to base column structure.

    Allows perturbing T, Ts, q, co2, and **insolation** independently.
    """
    cs = climlab.column_state(lev=s_base['lev'], water_depth=1.0)
    cs['Tatm'][:] = s_base['Tatm']
    cs['Ts'][:] = s_base['Ts']
    h2o = climlab.radiation.ManabeWaterVapor(state=cs, name='h2o')
    rad = climlab.radiation.RRTMG(
        state=cs, specific_humidity=s_base['q'].copy(),
        albedo=s_base['albedo'], insolation=s_base['insolation'], icld=0,
        name='rad',
    )
    rad.absorber_vmr['CO2'] = s_base['co2']
    rad.absorber_vmr['O3'][:] = s_base['o3']

    base_q = s_base['q'].copy()
    base_co2 = s_base['co2']
    base_insol = s_base['insolation']

    def eval_RS(Tatm, Ts, q=None, co2=None, insolation=None):
        rad.Tatm[:] = Tatm
        rad.Ts[:] = Ts
        rad.specific_humidity[:] = q if q is not None else base_q
        rad.absorber_vmr['CO2'] = co2 if co2 is not None else base_co2
        rad.insolation = insolation if insolation is not None else base_insol
        rad.compute_diagnostics()
        R = R_per_layer(np.array(rad.LW_flux_up).flatten(),
                        np.array(rad.LW_flux_down).flatten())
        S = S_per_layer(np.array(rad.SW_flux_up).flatten(),
                        np.array(rad.SW_flux_down).flatten())
        return R, S

    def eval_R(Tatm, Ts):
        return eval_RS(Tatm, Ts)[0]

    return eval_RS, eval_R


def decompose_solar(s1, s2):
    """Lu/Cai-style decomposition for +1% solar perturbation."""
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
    R_base, S_base = eval_RS(s1['Tatm'], s1['Ts'])

    # ===== Linearized partial perturbations (T held at T1) =====
    # solar: only insolation → warm; everything else (q, co2, T) base
    R_solar, S_solar = eval_RS(s1['Tatm'], s1['Ts'],
                               insolation=s2['insolation'])
    # q: only q → warm; T, insolation, co2 stay base
    R_q, S_q = eval_RS(s1['Tatm'], s1['Ts'], q=s2['q'])
    # co2: in this experiment co2 unchanged, but we evaluate to confirm ≈ 0
    R_co2, S_co2 = eval_RS(s1['Tatm'], s1['Ts'], co2=s2['co2'])
    # warm (all together except T at T1)
    R_warm, S_warm = eval_RS(s1['Tatm'], s1['Ts'],
                             q=s2['q'], co2=s2['co2'],
                             insolation=s2['insolation'])

    frc = {
        'solar': R_solar - R_base,
        'q':     R_q     - R_base,
        'co2':   R_co2   - R_base,
        'warm':  R_warm  - R_base,
    }
    dS = {
        'solar': S_solar - S_base,
        'q':     S_q     - S_base,
        'co2':   S_co2   - S_base,
        'warm':  S_warm  - S_base,
    }

    # ===== ΔQ_dry from FULL state energy balance (no LH/SH in RCE) =====
    R_full2, S_full2 = eval_RS(s2['Tatm'], s2['Ts'],
                               q=s2['q'], co2=s2['co2'],
                               insolation=s2['insolation'])
    Q_dry_1 = R_base   - S_base
    Q_dry_2 = R_full2 - S_full2
    dQ_dry = Q_dry_2 - Q_dry_1

    # ===== Solve partial dT_X =====
    dT = {
        'solar': solve_dT(drdt_inv, frc['solar'], dS['solar']),
        'q':     solve_dT(drdt_inv, frc['q'],     dS['q']),
        'co2':   solve_dT(drdt_inv, frc['co2'],   dS['co2']),
        'dry':   solve_dT(drdt_inv, np.zeros_like(dQ_dry), dQ_dry),
        'warm':  solve_dT(drdt_inv, frc['warm'],  dS['warm']),
    }

    dT_obs = np.concatenate([s2['Tatm'] - s1['Tatm'], [s2['Ts'] - s1['Ts']]])
    dT_sum = dT['solar'] + dT['q'] + dT['dry']

    # Lump-style atmdyn for historical comparison (also include co2 should be ~0)
    frc_atm, dS_atm = lump_atm_residual(frc['warm'], dS['warm'])
    dT['atmdyn_lump'] = solve_dT(drdt_inv, frc_atm, dS_atm)

    return {
        'lev': np.concatenate([s1['lev'], [1013.0]]),
        'frc': frc, 'dS': dS, 'dT': dT,
        'dQ_dry': dQ_dry, 'Q_dry_1': Q_dry_1, 'Q_dry_2': Q_dry_2,
        'dT_obs': dT_obs, 'dT_sum_lucai': dT_sum,
        'drdt': drdt, 'R_base': R_base, 'S_base': S_base,
    }


def report(res):
    obs = res['dT_obs']
    sum_lc = res['dT_sum_lucai']

    print('\n===== Surface dT (Lu/Cai Eq.4 decomposition) =====')
    for k in ['solar', 'q', 'co2', 'dry', 'warm']:
        print('  ΔT_%-6s = %+.4f K' % (k, res['dT'][k][-1]))
    print('  Σ Lu/Cai     = %+.4f K  (solar + q + dry)' % sum_lc[-1])
    print('  observed     = %+.4f K' % obs[-1])
    resid = obs[-1] - sum_lc[-1]
    pct = 100. * resid / obs[-1] if obs[-1] != 0 else 0
    print('  residual     = %+.4f K  (%.2f%% of obs)' % (resid, pct))

    print('\n===== Profile rmsd (Σ partials vs observed) =====')
    print('  Lu/Cai:  rmsd=%.4f K, max|diff|=%.4f K' %
          (np.sqrt(np.mean((sum_lc - obs)**2)),
           np.max(np.abs(sum_lc - obs))))

    print('\n===== Validation criteria =====')
    obs_sfc = obs[-1]
    sol_sfc = res['dT']['solar'][-1]
    q_sfc = res['dT']['q'][-1]
    co2_sfc = res['dT']['co2'][-1]
    dry_sfc = res['dT']['dry'][-1]
    print('  [1] Closure: residual %.2f%% of obs (target <1%%) -> %s' %
          (abs(pct), 'PASS' if abs(pct) < 1.0 else 'FAIL'))
    sol_share = 100. * sol_sfc / obs_sfc if obs_sfc != 0 else 0
    print('  [2] Solar dominates: ΔT_solar = %.1f%% of ΔT_obs (target >50%%) -> %s' %
          (sol_share, 'PASS' if sol_share > 50 else 'FAIL'))
    print('  [3] CO2 absent: ΔT_co2 = %+.4f K (target |.|<0.005 K) -> %s' %
          (co2_sfc, 'PASS' if abs(co2_sfc) < 0.005 else 'FAIL'))
    print('  [4] H2O positive feedback: ΔT_q = %+.4f K (target >0) -> %s' %
          (q_sfc, 'PASS' if q_sfc > 0 else 'FAIL'))


def save_nc(res, path):
    nc = Dataset(path, 'w')
    nc.createDimension('lev', len(res['lev']))
    nc.createVariable('lev', 'f8', ('lev',))[:] = res['lev']
    for name, arr in [
        ('dT_obs', res['dT_obs']),
        ('dT_sum_lucai', res['dT_sum_lucai']),
        ('dT_solar',  res['dT']['solar']),
        ('dT_q',      res['dT']['q']),
        ('dT_co2',    res['dT']['co2']),
        ('dT_dry',    res['dT']['dry']),
        ('dT_warm',   res['dT']['warm']),
        ('dT_atmdyn_lump', res['dT']['atmdyn_lump']),
        ('frc_solar', res['frc']['solar']), ('dS_solar', res['dS']['solar']),
        ('frc_q',     res['frc']['q']),     ('dS_q',     res['dS']['q']),
        ('frc_co2',   res['frc']['co2']),   ('dS_co2',   res['dS']['co2']),
        ('frc_warm',  res['frc']['warm']),  ('dS_warm',  res['dS']['warm']),
        ('dQ_dry',    res['dQ_dry']),
        ('drdt_diag', np.diag(res['drdt'])),
        ('R_base',    res['R_base']),       ('S_base',   res['S_base']),
    ]:
        nc.createVariable(name, 'f8', ('lev',))[:] = arr
    nc.close()
    print('Saved:', path)


def plot_fig(res, path):
    plev = res['lev']
    fig, axes = plt.subplots(1, 5, figsize=(20, 6), sharey=True)

    panels = [
        ('Total ΔT', [
            (res['dT_obs'],       'climlab obs',     'black',  '-',  'o'),
            (res['dT_sum_lucai'], 'CFRAM Lu/Cai Σ',  'red',    '--', 's'),
        ]),
        ('ΔT_solar (TOA solar)', [(res['dT']['solar'], None, 'red',    '-', 'o')]),
        ('ΔT_q (water vapor fb)', [(res['dT']['q'],     None, 'blue',   '-', 'o')]),
        ('ΔT_co2 (sanity ≈0)',    [(res['dT']['co2'],   None, 'gray',   '-', 'o')]),
        ('ΔT_dry (Lu/Cai ΔQ_dry)', [(res['dT']['dry'],   None, 'purple', '-', 'o')]),
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
    fig.suptitle('Lu/Cai CFRAM decomposition: climlab RCE +1% TOA solar', fontsize=12)
    plt.tight_layout()
    plt.savefig(path, dpi=130, bbox_inches='tight')
    print('Saved:', path)


def main():
    s1 = load_rce_solar('base')
    s2 = load_rce_solar('warm')
    print('Base solar=%.3f W/m², Ts=%.4f K' % (s1['insolation'], s1['Ts']))
    print('Warm solar=%.3f W/m², Ts=%.4f K  (ΔTs=%+.4f K, Δsolar=%+.3f W/m²)' %
          (s2['insolation'], s2['Ts'], s2['Ts'] - s1['Ts'],
           s2['insolation'] - s1['insolation']))

    res = decompose_solar(s1, s2)
    report(res)
    save_nc(res, os.path.join(OUTDIR, 'cfram_partials_solar.nc'))
    plot_fig(res, os.path.join(OUTDIR, 'cfram_partials_solar.png'))


if __name__ == '__main__':
    main()
