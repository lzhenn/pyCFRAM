"""CFRAM-A full decomposition workflow.

Orchestrates the complete CFRAM analysis:
1. Define base and anomalous states
2. Run radiation for base state
3. Run radiation for each perturbed state (one variable at a time)
4. Build Planck feedback matrix
5. Compute forcing = perturbed - base
6. Solve for partial temperature changes

This module is backend-agnostic: the actual radiation calls are done
through a callback function, which can be backed by subprocess (Fortran)
or eventually pure f2py/Python.
"""

import numpy as np
from .radiation import (compute_interface_levels, compute_flux_convergence,
                        prepare_rrtmg_inputs)
from .planck_matrix import build_planck_matrix, solve_cfram_multi
from .constants import CH4_PPMV, N2O_PPMV, ICLD_ON, IAER_ON


class CFRAMDecomposition:
    """Single-column CFRAM-A decomposition."""

    # Perturbation terms: name -> which variables change from base to warm
    TERMS = {
        'co2':     ['co2'],
        'q':       ['q'],
        'ts':      ['ts'],
        'o3':      ['o3'],
        'solar':   ['zenith'],
        'albedo':  ['albedo_sw'],
        'cloud':   ['cldfrac', 'cldlwc', 'cldiwc'],
        'aerosol': ['tauaer_sw', 'ssaaer_sw', 'asmaer_sw', 'tauaer_lw'],
    }

    def __init__(self, rad_func, rad_lw_func):
        """Initialize decomposition.

        Args:
            rad_func: callable(state_dict) -> (rad, lw, sw, fd, fu, fdl, ful)
                Full radiation (LW+SW). Returns flux convergences and fluxes
                in ERA order. state_dict has same keys as prepare_rrtmg_inputs output.

                Returns:
                - rad: (nlayer,) total flux convergence
                - lw: (nlayer,) LW convergence
                - sw: (nlayer,) SW convergence
                - fd: (nlayer+1,) total downward flux
                - fu: (nlayer+1,) total upward flux
                - fdl: (nlayer+1,) LW downward flux
                - ful: (nlayer+1,) LW upward flux

            rad_lw_func: callable(state_dict) -> (lw, fdl, ful)
                LW-only radiation for Planck matrix construction.
        """
        self.rad_func = rad_func
        self.rad_lw_func = rad_lw_func

    def run(self, base_state, warm_state):
        """Run full CFRAM decomposition.

        Args:
            base_state: dict with keys matching Column fields:
                plev, t, q, o3, co2, ts, ps, zenith, albedo_sw, albedo_lw,
                cldfrac, cldlwc, cldiwc, tauaer_sw, ssaaer_sw, asmaer_sw,
                tauaer_lw
            warm_state: same structure, warm/anomalous period values

        Returns:
            results: dict with:
                'nlayer': number of active layers
                'plev': pressure levels used
                'frc_{term}': (nlayer+1,) forcing for each term
                'dT_{term}': (nlayer,) partial temperature change
                'frc_warm': (nlayer+1,) total forcing
                'dT_warm': (nlayer,) total temperature change
                'drdt': (nlayer+1, nlayer+1) Planck matrix
        """
        nlayer = self._determine_nlayer(base_state['plev'], base_state['ps'])

        # Prepare RRTMG inputs for base state
        base_inputs = self._prepare_inputs(base_state, nlayer)

        # 1. Run base state radiation
        rad_base, lw_base, sw_base, fd_base, fu_base, fdl_base, ful_base = \
            self.rad_func(base_inputs)

        # 2. Run warm state radiation (all variables changed)
        warm_inputs = self._prepare_inputs(warm_state, nlayer)
        rad_warm, lw_warm, sw_warm, fd_warm, fu_warm, fdl_warm, ful_warm = \
            self.rad_func(warm_inputs)

        # 3. Build Planck matrix (LW only, base state)
        drdt = build_planck_matrix(
            self.rad_lw_func, base_inputs,
            lw_base, fdl_base, ful_base, nlayer)

        # 4. Compute forcings for each perturbation term
        forcings = {}

        # Total (warm) forcing
        frc_warm = self._compute_forcing(
            rad_warm, fd_warm, fu_warm, rad_base, fd_base, fu_base, nlayer)
        forcings['warm'] = frc_warm

        # Per-term forcings
        for term, changed_vars in self.TERMS.items():
            perturbed_state = self._make_perturbed_state(
                base_state, warm_state, changed_vars)
            pert_inputs = self._prepare_inputs(perturbed_state, nlayer)
            rad_pert, _, _, fd_pert, fu_pert, _, _ = self.rad_func(pert_inputs)

            frc = self._compute_forcing(
                rad_pert, fd_pert, fu_pert, rad_base, fd_base, fu_base, nlayer)
            forcings[term] = frc

        # 5. Solve for partial temperature changes
        dT_dict = solve_cfram_multi(drdt, forcings, nlayer)

        # Package results
        results = {
            'nlayer': nlayer,
            'plev': base_state['plev'][:nlayer],
            'drdt': drdt,
        }
        for term in list(self.TERMS.keys()) + ['warm']:
            results[f'frc_{term}'] = forcings[term]
            results[f'dT_{term}'] = dT_dict[term]

        return results

    def _determine_nlayer(self, plev, ps):
        """Find number of active layers (pressure levels above surface)."""
        ps_hpa = ps / 100.0 if ps > 2000 else ps  # handle Pa or hPa
        nlayer = 0
        for i, p in enumerate(plev):
            if p < ps_hpa:
                nlayer = i + 1
            else:
                break
        return nlayer

    def _prepare_inputs(self, state, nlayer):
        """Convert column state dict to RRTMG input dict."""
        return prepare_rrtmg_inputs(
            nlayer=nlayer,
            plev=state['plev'],
            t=state['t'],
            q=state['q'],
            o3=state['o3'],
            co2_ppmv=state['co2'],
            ch4_ppmv=state.get('ch4', CH4_PPMV),
            n2o_ppmv=state.get('n2o', N2O_PPMV),
            ps=state['ps'],
            ts=state['ts'],
            zenith=state['zenith'],
            albedo_sw=state['albedo_sw'],
            albedo_lw=state.get('albedo_lw', 0.0),
            cldfrac=state['cldfrac'],
            cldlwc=state['cldlwc'],
            cldiwc=state['cldiwc'],
            tauaer_sw=state.get('tauaer_sw'),
            ssaaer_sw=state.get('ssaaer_sw'),
            asmaer_sw=state.get('asmaer_sw'),
            tauaer_lw=state.get('tauaer_lw'),
            icld=ICLD_ON,
            iaer=IAER_ON,
        )

    def _compute_forcing(self, rad_pert, fd_pert, fu_pert,
                         rad_base, fd_base, fu_base, nlayer):
        """Compute radiative forcing vector (nlayer+1).

        Atmospheric layers: frc = rad_pert - rad_base
        Surface: frc = net_flux_pert(sfc) - net_flux_base(sfc)
        """
        frc = np.zeros(nlayer + 1)
        frc[:nlayer] = rad_pert[:nlayer] - rad_base[:nlayer]
        frc[nlayer] = ((fd_pert[nlayer] - fu_pert[nlayer]) -
                       (fd_base[nlayer] - fu_base[nlayer]))
        return frc

    def _make_perturbed_state(self, base, warm, changed_vars):
        """Create a perturbed state: base values except for changed_vars from warm."""
        state = {}
        for key in base:
            if key in changed_vars:
                state[key] = warm[key]
            else:
                state[key] = base[key]
        return state

    def print_summary(self, results):
        """Print decomposition summary."""
        nlayer = results['nlayer']
        print(f"CFRAM Decomposition: {nlayer} active layers")
        print(f"{'Term':<10} {'frc(atm) range':>24}  {'dT range':>24}")
        print("-" * 62)
        for term in list(self.TERMS.keys()) + ['warm']:
            frc = results[f'frc_{term}']
            dT = results[f'dT_{term}']
            print(f"{term:<10} [{frc[:nlayer].min():>10.4f}, {frc[:nlayer].max():>10.4f}]  "
                  f"[{dT.min():>10.4f}, {dT.max():>10.4f}]")

        # Additivity check
        dT_sum = sum(results[f'dT_{t}'] for t in self.TERMS)
        dT_warm = results['dT_warm']
        residual = np.max(np.abs(dT_sum - dT_warm))
        print(f"\nAdditivity: max|sum - warm| = {residual:.4f} K")
