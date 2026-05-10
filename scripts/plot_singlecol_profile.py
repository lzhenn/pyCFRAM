#!/usr/bin/env python3
"""5-subplot vertical-profile decomposition figure for 1×1 single-column cases.

Identity (clear-sky climlab RCE):
    dT_obs ≈ dT_co2 + dT_q + dT_ts + dT_dyn_proper
where dT_dyn_proper = dT_obs − Σ(radiative) is the genuine non-radiative
residual (convective adjustment + transport — by construction closes exactly).

Layout
------
[1] Total ΔT  : dT_observed (RT-independent) + Σ(co2+q+ts) per engine
                (lines should overlap closely; gap = `dT_dyn_proper`)
[2] CO2       : dT_co2 — RRTMG vs Fu
[3] WV        : dT_q   — RRTMG vs Fu
[4] Surf-T    : dT_ts  — RRTMG vs Fu  (surface-emission radiative response)
[5] Dynamics  : dT_dyn_proper = dT_obs − (co2+q+ts) — RRTMG vs Fu
                (convective heat redistribution + LH/SH residual; column
                integral ≈ 0 in pure RCE)

Usage
-----
    python3 scripts/plot_singlecol_profile.py climlab_4xco2 climlab_4xco2_fu

The first case is "primary" (its figures/ dir gets the output).
"""
import os
import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case


# Display label for radiation.scheme key
SCHEME_LABEL = {'fu': 'Fu', 'rrtmg': 'RRTMG'}
SCHEME_COLOR = {'fu': '#d62728', 'rrtmg': '#1f77b4'}    # red, blue


def _load(case_name):
    cfg = load_case(case_name)
    out = os.path.join(cfg['_output_dir'], 'cfram_result.nc')
    if not os.path.exists(out):
        raise SystemExit("Missing %s — run `python3 run_case.py %s --step run` first."
                         % (out, case_name))
    nc = Dataset(out)
    lev = np.array(nc.variables['lev'][:])           # sfc→TOA atm + ground at end
    nlev = len(lev) - 1                              # ground is last
    scheme = (cfg.get('radiation', {}).get('scheme')
              or 'rrtmg')                            # default

    def get(name):
        if name not in nc.variables:
            return None
        return np.array(nc.variables[name][:, 0, 0])

    data = {k: get('dT_' + k) for k in
            ['observed', 'co2', 'q', 'ts', 'dry', 'atmdyn']}
    nc.close()
    return cfg, lev, nlev, scheme, data


def _plot_panel(ax, title, lev, atm_idx, sfc_idx, curves, xlim=None):
    """One panel: plot each (label, color, dT) over (atm profile + sfc marker).

    atm_idx  : slice(0, nlev)   — atmospheric layers (lev sfc→TOA)
    sfc_idx  : nlev             — ground (single point)
    """
    for label, color, dT, linestyle in curves:
        if dT is None:
            continue
        ax.plot(dT[atm_idx], lev[atm_idx], color=color, lw=1.6,
                ls=linestyle, label=label, zorder=3)
        ax.scatter(dT[sfc_idx], lev[sfc_idx], color=color, marker='*',
                   s=80, zorder=4, edgecolor='k', linewidth=0.4)

    ax.axvline(0, color='gray', lw=0.5, ls='--', alpha=0.7)
    ax.set_title(title, fontsize=11)
    ax.set_xlabel(r'$\Delta T$ (K)', fontsize=10)
    ax.set_ylim(lev.max() + 5, lev.min() - 10)       # invert: TOA at top
    if xlim is not None:
        ax.set_xlim(xlim)
    ax.grid(True, alpha=0.3)


def main(case_names):
    if not case_names:
        sys.exit("Usage: plot_singlecol_profile.py <case1> [<case2> ...]")

    cases = []
    for cn in case_names:
        cfg, lev, nlev, scheme, data = _load(cn)
        cases.append({'name': cn, 'cfg': cfg, 'lev': lev, 'nlev': nlev,
                      'scheme': scheme, 'data': data})

    # All cases must share the same vertical grid (climlab pair => identical
    # grid by construction).
    lev_ref = cases[0]['lev']
    nlev = cases[0]['nlev']
    for c in cases[1:]:
        if not np.allclose(c['lev'], lev_ref):
            raise SystemExit("Cases have mismatched lev grids — cannot overlay.")

    atm = slice(0, nlev)
    sfc = nlev
    primary = cases[0]
    desc = primary['cfg'].get('description', primary['name'])

    fig, axes = plt.subplots(1, 5, figsize=(19, 5), sharey=True)

    # --- Panel 1: Total — dT_observed + Σ(co2+q+ts) per engine ---
    # If dT_ts is missing (legacy build), Σ degenerates to (co2+q+dry).
    obs = primary['data']['observed']
    panel1 = [('observed (input)', 'black', obs, '-')]
    for c in cases:
        sch = c['scheme']
        ts = c['data'].get('ts')
        if ts is not None:
            s = c['data']['co2'] + c['data']['q'] + ts
            label = 'Σ %s = co2+q+ts' % SCHEME_LABEL.get(sch, sch)
        else:
            s = c['data']['co2'] + c['data']['q'] + c['data']['dry']
            label = 'Σ %s = co2+q+dry (legacy)' % SCHEME_LABEL.get(sch, sch)
        panel1.append((label, SCHEME_COLOR.get(sch, 'gray'), s, '--'))
    _plot_panel(axes[0], '(a) Total ΔT  (radiative closure)',
                lev_ref, atm, sfc, panel1)
    axes[0].set_ylabel('Pressure (hPa)', fontsize=10)
    axes[0].legend(loc='best', fontsize=8, framealpha=0.9)

    # --- Panels 2–4: dT_co2, dT_q, dT_ts per engine ---
    for col, (key, label) in enumerate([
        ('co2', '(b) ΔT$_{CO_2}$'),
        ('q',   '(c) ΔT$_{WV}$'),
        ('ts',  '(d) ΔT$_{T_s}$  (surface T radiative response)'),
    ], start=1):
        curves = []
        for c in cases:
            sch = c['scheme']
            dT = c['data'].get(key)
            curves.append((SCHEME_LABEL.get(sch, sch),
                           SCHEME_COLOR.get(sch, 'gray'), dT, '-'))
        _plot_panel(axes[col], label, lev_ref, atm, sfc, curves)
        axes[col].legend(loc='best', fontsize=8, framealpha=0.9)

    # --- Panel 5: True non-radiative residual ΔT_dyn = dT_obs − (co2+q+ts) ---
    panel5 = []
    for c in cases:
        sch = c['scheme']
        ts = c['data'].get('ts')
        if ts is None:
            continue
        dyn = c['data']['observed'] - c['data']['co2'] - c['data']['q'] - ts
        panel5.append((SCHEME_LABEL.get(sch, sch),
                       SCHEME_COLOR.get(sch, 'gray'), dyn, '-'))
    _plot_panel(axes[4], '(e) ΔT$_{dyn}$  (convective + transport)',
                lev_ref, atm, sfc, panel5)
    axes[4].legend(loc='best', fontsize=8, framealpha=0.9)

    # --- Title ---
    title = ('pyCFRAM single-column decomposition: %s' %
             ', '.join(c['name'] for c in cases))
    fig.suptitle(title, fontsize=12)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Save into the primary case's figures/ dir
    fig_dir = primary['cfg']['_figures_dir']
    os.makedirs(fig_dir, exist_ok=True)
    out_path = os.path.join(fig_dir, 'fig_singlecol_profile.png')
    plt.savefig(out_path, dpi=130, bbox_inches='tight')
    print('Saved: %s' % out_path)

    # If two engines were given, mirror into the second case's figures/ too
    # (convenient when both engines are runnable).
    if len(cases) > 1:
        fig_dir2 = cases[1]['cfg']['_figures_dir']
        os.makedirs(fig_dir2, exist_ok=True)
        out_path2 = os.path.join(fig_dir2, 'fig_singlecol_profile.png')
        plt.savefig(out_path2, dpi=130, bbox_inches='tight')
        print('Saved: %s' % out_path2)


if __name__ == '__main__':
    main(sys.argv[1:])
