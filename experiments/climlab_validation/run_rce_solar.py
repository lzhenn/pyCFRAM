#!/usr/bin/env python3
"""Run RCE base + 1.01×TOA-solar perturbation, on pyCFRAM's 37-level plev grid.

This produces a PURE solar-perturbation pair (CO2 fixed at 348 ppm, no cloud,
no aerosol, ManabeWaterVapor RH-fixed), suitable for validating both:

- Path A: cfram_decompose_solar.py (pure-python CFRAM)
- Path B: production pyCFRAM via 1×1 case (NetCDF saved in pyCFRAM input format)

Output (in `output/`)
---------------------
- rce_solar_base.nc  / rce_solar_warm.nc       : climlab native (lev TOA→sfc)
- pycfram_input/base_pres.nc + base_surf.nc    : pyCFRAM 1×1 case input (sfc→TOA)
- pycfram_input/perturbed_pres.nc + perturbed_surf.nc

Setup
-----
- 37 layers matching pyCFRAM defaults.yaml plev
- Slab surface 1 m water depth
- ManabeWaterVapor (q follows fixed RH), no cloud, no aerosol
- CO2 = 348 ppm (fixed in both runs)
- Albedo = 0.3
- Insolation: base=341.5 W/m², warm=341.5×1.01=345.015 W/m²
- ConvectiveAdjustment (dry, lapse_rate ≤ 6.5 K/km)

Convergence: |ASR-OLR| < 0.02 W/m²
"""
import os
import sys
import numpy as np
from netCDF4 import Dataset
import climlab


HERE = os.path.dirname(os.path.abspath(__file__))

# pyCFRAM 37-level plev (TOA → surface)
PLEV = np.array([
    1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150,
    175, 200, 225, 250, 300, 350, 400, 450, 500, 550,
    600, 650, 700, 750, 775, 800, 825, 850, 875, 900,
    925, 950, 975, 1000
], dtype=float)
NLEV = len(PLEV)

CO2_PPM = 348.0
ALBEDO = 0.3
INSOLATION_BASE = 341.5
SOLAR_FACTOR = 1.01


def build_model(insolation, name):
    """climlab RCE on pyCFRAM 37-level grid."""
    state = climlab.column_state(lev=PLEV, water_depth=1.0)
    state['Tatm'][:] = np.linspace(200., 290., NLEV)
    state['Ts'][:] = 290.

    h2o = climlab.radiation.ManabeWaterVapor(state=state, name='h2o')
    rad = climlab.radiation.RRTMG(
        state=state, specific_humidity=h2o.q,
        albedo=ALBEDO, insolation=insolation, icld=0, name='rad',
    )
    rad.absorber_vmr['CO2'] = CO2_PPM * 1e-6

    conv = climlab.convection.ConvectiveAdjustment(
        state=state, adj_lapse_rate=6.5, name='dry_conv',
    )
    model = climlab.couple([rad, h2o, conv], name=name)
    return model, h2o, rad


def integrate_to_equilibrium(model, rad, max_days=2000, tol=0.02, chunk=20):
    elapsed = 0
    while elapsed < max_days:
        model.integrate_days(chunk)
        elapsed += chunk
        asr = float(np.array(rad.ASR).flat[0])
        olr = float(np.array(rad.OLR).flat[0])
        imbalance = asr - olr
        ts = float(np.array(model.Ts).flat[0])
        print('  day %d: Ts=%.4f K, ASR-OLR=%+.4f W/m²' %
              (elapsed, ts, imbalance))
        if abs(imbalance) < tol:
            return elapsed, imbalance
    return elapsed, imbalance


def save_native(model, h2o, rad, path, insolation):
    """Save climlab native state (lev TOA→sfc) for path A."""
    Tatm = np.array(model.Tatm).flatten()       # (NLEV,) TOA→sfc
    Ts = float(np.array(model.Ts).flat[0])
    q = np.array(h2o.q).flatten()
    o3 = np.array(rad.absorber_vmr['O3']).flatten()  # vmr (mol/mol)

    nc = Dataset(path, 'w')
    nc.createDimension('lev', NLEV)
    nc.createDimension('half', NLEV + 1)
    nc.co2_ppm = CO2_PPM
    nc.insolation_Wm2 = insolation
    nc.albedo = ALBEDO
    olr = float(np.array(rad.OLR).flat[0])
    asr = float(np.array(rad.ASR).flat[0])
    nc.toa_imbalance_Wm2 = asr - olr
    nc.olr_Wm2 = olr
    nc.asr_Wm2 = asr

    nc.createVariable('lev', 'f8', ('lev',))[:] = PLEV
    nc.createVariable('Tatm', 'f8', ('lev',))[:] = Tatm
    nc.createVariable('q', 'f8', ('lev',))[:] = q
    nc.createVariable('o3', 'f8', ('lev',))[:] = o3
    nc.createVariable('Ts', 'f8', ())[...] = Ts
    nc.createVariable('LW_flux_up',   'f8', ('half',))[:] = np.array(rad.LW_flux_up).flatten()
    nc.createVariable('LW_flux_down', 'f8', ('half',))[:] = np.array(rad.LW_flux_down).flatten()
    nc.createVariable('SW_flux_up',   'f8', ('half',))[:] = np.array(rad.SW_flux_up).flatten()
    nc.createVariable('SW_flux_down', 'f8', ('half',))[:] = np.array(rad.SW_flux_down).flatten()
    nc.close()
    print('Saved %s  (Ts=%.4f K, OLR=%.3f, ASR=%.3f)' % (path, Ts, olr, asr))


def save_pycfram_input(model, h2o, rad, insolation, out_pres, out_surf):
    """Save in pyCFRAM input NetCDF format (1×1 grid, lev sfc→TOA convention).

    pyCFRAM run_parallel_python.py reads with `[0, ::-1, :, :]`, so we store
    arrays in surface→TOA order; after [::-1] they become TOA→surface internally.
    """
    Tatm = np.array(model.Tatm).flatten()       # TOA→sfc
    Ts = float(np.array(model.Ts).flat[0])
    q = np.array(h2o.q).flatten()
    o3 = np.array(rad.absorber_vmr['O3']).flatten()

    # Surface→TOA in NetCDF (so pyCFRAM's [::-1] gets back to TOA→sfc)
    plev_sfc2top = PLEV[::-1]
    Tatm_sfc2top = Tatm[::-1]
    q_sfc2top = q[::-1]
    o3_sfc2top = o3[::-1]

    nc = Dataset(out_pres, 'w')
    nc.createDimension('time', 1)
    nc.createDimension('lev', NLEV)
    nc.createDimension('lat', 1)
    nc.createDimension('lon', 1)
    nc.createVariable('time', 'f8', ('time',))[:] = [0.0]
    nc.createVariable('lev',  'f8', ('lev',))[:]  = plev_sfc2top
    nc.createVariable('lat',  'f8', ('lat',))[:]  = [0.0]
    nc.createVariable('lon',  'f8', ('lon',))[:]  = [0.0]

    def write3d(name, arr_sfc2top):
        v = nc.createVariable(name, 'f8', ('time', 'lev', 'lat', 'lon'))
        v[:] = arr_sfc2top[None, :, None, None]

    write3d('ta_lay', Tatm_sfc2top)
    write3d('q', q_sfc2top)
    write3d('o3', o3_sfc2top)
    # CO2 as full 3D field of constant value (volume mixing ratio mol/mol)
    write3d('co2', np.full(NLEV, CO2_PPM * 1e-6))
    # Cloud + aerosols all zero (no cloud, no aerosol RCE)
    for name in ['camt', 'cliq', 'cice', 'bc', 'ocphi', 'ocpho', 'sulf', 'ss', 'dust']:
        write3d(name, np.zeros(NLEV))
    nc.close()

    nc = Dataset(out_surf, 'w')
    nc.createDimension('time', 1)
    nc.createDimension('lat', 1)
    nc.createDimension('lon', 1)
    nc.createVariable('time', 'f8', ('time',))[:] = [0.0]
    nc.createVariable('lat',  'f8', ('lat',))[:]  = [0.0]
    nc.createVariable('lon',  'f8', ('lon',))[:]  = [0.0]

    def write2d(name, val):
        v = nc.createVariable(name, 'f8', ('time', 'lat', 'lon'))
        v[:] = np.full((1, 1, 1), val)

    write2d('ts', Ts)
    # ps slightly above lowest atm plev (1000 hPa) to avoid zero-thickness
    # bottom layer in pyCFRAM Fortran (nlayer would otherwise be 36 not 37,
    # leaving uninitialized -999 at the lowest atm row).
    write2d('ps', 101325.0)              # Pa (= 1013.25 hPa, standard sea level)
    write2d('solar', insolation)         # W/m²
    write2d('albedo', ALBEDO)
    nc.close()
    print('Saved pyCFRAM input: %s  +  %s' % (out_pres, out_surf))


def main():
    out_dir = os.path.join(HERE, 'output')
    pycfram_dir = os.path.join(out_dir, 'pycfram_input')
    os.makedirs(pycfram_dir, exist_ok=True)

    cases = [
        ('base', INSOLATION_BASE),
        ('warm', INSOLATION_BASE * SOLAR_FACTOR),
    ]
    for label, insolation in cases:
        print('\n===== %s : insolation = %.4f W/m² =====' %
              (label.upper(), insolation))
        model, h2o, rad = build_model(insolation, name='RCE_solar_' + label)
        elapsed, imbalance = integrate_to_equilibrium(model, rad)
        print('  -> equilibrium reached after %d days (imbalance %.4f W/m²)' %
              (elapsed, imbalance))

        save_native(model, h2o, rad,
                    os.path.join(out_dir, 'rce_solar_%s.nc' % label),
                    insolation)
        # pyCFRAM input naming: base→base_pres/base_surf, warm→perturbed_*
        prefix = 'base' if label == 'base' else 'perturbed'
        save_pycfram_input(
            model, h2o, rad, insolation,
            os.path.join(pycfram_dir, '%s_pres.nc' % prefix),
            os.path.join(pycfram_dir, '%s_surf.nc' % prefix),
        )


if __name__ == '__main__':
    main()
