#!/usr/bin/env python3
"""Run radiative-convective equilibrium with climlab at 1×CO2 and 2×CO2.

Output: rce_1xco2.nc, rce_2xco2.nc — clean equilibrium states for CFRAM
        decomposition validation (route B).

Setup is climlab-standard; we are NOT trying to match Lu/Cai 2009 paper
parameters. The point is to provide an internally-consistent atmosphere column
where pyCFRAM's CFRAM math can be tested against a known-correct reference.

Components
----------
- 30-layer column (climlab default grid, sigma-based)
- Slab surface with 1 m water depth (fast convergence)
- ManabeWaterVapor: q follows fixed RH profile
- RRTMG_LW + RRTMG_SW: full radiation (no clouds)
- ConvectiveAdjustment: dry, lapse_rate ≤ 6.5 K/km
- Insolation 341.5 W/m² (global mean), albedo 0.3

CO2 perturbation
----------------
- 1×CO2: 348 ppm (RRTMG default; close to 1990s level)
- 2×CO2: 696 ppm

Convergence
-----------
Integrate until |ASR - OLR| < 0.02 W/m² (typically 200-400 days).
"""
import os
import sys
import numpy as np
from netCDF4 import Dataset
import climlab


HERE = os.path.dirname(os.path.abspath(__file__))


def build_model(co2_ppm, name):
    """Build a climlab RCE single-column model."""
    state = climlab.column_state(num_lev=30, water_depth=1.0)
    nlev = state['Tatm'].shape[0]
    state['Tatm'][:] = np.linspace(200., 290., nlev)
    state['Ts'][:] = 290.

    h2o = climlab.radiation.ManabeWaterVapor(state=state, name='h2o')

    # Build with full default absorber set, then override CO2 only
    rad = climlab.radiation.RRTMG(
        state=state, specific_humidity=h2o.q,
        albedo=0.3, insolation=341.5, icld=0,
        name='rad',
    )
    rad.absorber_vmr['CO2'] = co2_ppm * 1e-6

    conv = climlab.convection.ConvectiveAdjustment(
        state=state, adj_lapse_rate=6.5, name='dry_conv',
    )

    model = climlab.couple([rad, h2o, conv], name=name)
    return model, h2o, rad, conv


def integrate_to_equilibrium(model, rad, max_days=2000, tol=0.02, chunk=20):
    """Integrate until TOA imbalance < tol W/m² or max_days reached."""
    elapsed = 0
    while elapsed < max_days:
        model.integrate_days(chunk)
        elapsed += chunk
        asr = float(np.array(rad.ASR).flat[0])
        olr = float(np.array(rad.OLR).flat[0])
        imbalance = asr - olr
        ts = float(np.array(model.Ts).flat[0])
        print('  day %d: Ts=%.3f K, ASR-OLR=%+.4f W/m² (target |.|<%.2f)'
              % (elapsed, ts, imbalance, tol))
        if abs(imbalance) < tol:
            return elapsed, imbalance
    return elapsed, imbalance


def save_state(model, h2o, rad, path, co2_ppm):
    """Save equilibrium state + diagnostics to NetCDF."""
    lev = np.array(model.lev)               # (nlev,) hPa
    Tatm = np.array(model.Tatm).flatten()   # (nlev,)
    Ts = float(np.array(model.Ts).flat[0])
    q = np.array(h2o.q).flatten()           # specific humidity kg/kg

    # Radiation diagnostics — flux profiles at half-levels (nlev+1 elements,
    # TOA at index 0, surface at index nlev)
    LW_flux_up = np.array(rad.LW_flux_up).flatten()
    LW_flux_down = np.array(rad.LW_flux_down).flatten()
    SW_flux_up = np.array(rad.SW_flux_up).flatten()
    SW_flux_down = np.array(rad.SW_flux_down).flatten()
    asr = float(np.array(rad.ASR).flat[0])
    olr = float(np.array(rad.OLR).flat[0])

    nc = Dataset(path, 'w')
    nc.createDimension('lev', len(lev))
    nc.createDimension('half', len(LW_flux_up))
    nc.co2_ppm = co2_ppm
    nc.insolation_Wm2 = 341.5
    nc.albedo = 0.3
    nc.toa_imbalance_Wm2 = asr - olr
    nc.olr_Wm2 = olr
    nc.asr_Wm2 = asr

    v = nc.createVariable('lev', 'f8', ('lev',));  v.units = 'hPa';  v[:] = lev
    v = nc.createVariable('Tatm', 'f8', ('lev',)); v.units = 'K';    v[:] = Tatm
    v = nc.createVariable('q', 'f8', ('lev',));    v.units = 'kg/kg';v[:] = q
    v = nc.createVariable('Ts', 'f8', ());          v.units = 'K';    v[...] = Ts
    v = nc.createVariable('LW_flux_up', 'f8', ('half',));   v.units = 'W/m^2'; v[:] = LW_flux_up
    v = nc.createVariable('LW_flux_down', 'f8', ('half',)); v.units = 'W/m^2'; v[:] = LW_flux_down
    v = nc.createVariable('SW_flux_up', 'f8', ('half',));   v.units = 'W/m^2'; v[:] = SW_flux_up
    v = nc.createVariable('SW_flux_down', 'f8', ('half',)); v.units = 'W/m^2'; v[:] = SW_flux_down
    nc.close()
    print('Saved %s  (Ts=%.3f K, OLR=%.2f, ASR=%.2f)' % (path, Ts, olr, asr))


def main():
    out_dir = os.path.join(HERE, 'output')
    os.makedirs(out_dir, exist_ok=True)

    for co2_ppm, label in [(348., '1xco2'), (696., '2xco2')]:
        print('\n===== Building model: CO2 = %.1f ppm =====' % co2_ppm)
        model, h2o, rad, conv = build_model(co2_ppm, name='RCE_' + label)
        elapsed, imbalance = integrate_to_equilibrium(model, rad,
                                                      max_days=2000, tol=0.02)
        print('  -> equilibrium reached after %d days (imbalance %.4f W/m²)' %
              (elapsed, imbalance))
        save_state(model, h2o, rad,
                   os.path.join(out_dir, 'rce_%s.nc' % label), co2_ppm)


if __name__ == '__main__':
    main()
