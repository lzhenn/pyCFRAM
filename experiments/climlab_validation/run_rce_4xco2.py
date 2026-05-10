#!/usr/bin/env python3
"""Run climlab RCE at 1×CO2 and 4×CO2, save pyCFRAM 1×1 input directly.

Output
------
- output/rce_1xco2.nc / rce_4xco2.nc — climlab native equilibrium states
  (Tatm, q, fluxes; for inspection / pure-Python CFRAM validation)
- cases/climlab_4xco2/input/{base,perturbed}_{pres,surf}.nc — pyCFRAM
  standard 1×1 input format consumed by `run_case.py climlab_4xco2`
- cases/climlab_4xco2_fu/input/* — symlinks to the same files (Fu engine
  variant)

Setup
-----
- 30-layer column (climlab default sigma grid; `cfram_*_1col` infers
  nlev at runtime from data_prep/plev.dat size, so no Fortran rebuild)
- Slab surface, 1 m water depth
- ManabeWaterVapor: q follows fixed RH (water-vapor feedback in radiation only)
- RRTMG_LW + RRTMG_SW: clear-sky radiation (icld=0)
- ConvectiveAdjustment: dry, lapse_rate ≤ 6.5 K/km
- Insolation 341.5 W/m², albedo 0.3
- 1×CO2: 348 ppm; 4×CO2: 1392 ppm

Convergence
-----------
Integrate until |ASR-OLR| < 0.02 W/m². 4×CO2 typically needs ~600-1500 days.
"""
import os
import sys
import numpy as np
from netCDF4 import Dataset

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# Reuse build_model / integrate_to_equilibrium / save_state from run_rce.py
# (build_model uses climlab.column_state(num_lev=30) → native 30-level sigma grid).
from run_rce import build_model, integrate_to_equilibrium, save_state

# Match run_rce.py constants (kept private; if those change, update here).
ALBEDO = 0.3
INSOLATION = 341.5

CASE_DIR_PRIMARY = os.path.join(
    os.path.dirname(os.path.dirname(HERE)), 'cases', 'climlab_4xco2'
)
CASE_DIR_FU = os.path.join(
    os.path.dirname(os.path.dirname(HERE)), 'cases', 'climlab_4xco2_fu'
)


def save_pycfram_input(model, h2o, rad, co2_ppm, out_pres, out_surf):
    """Save equilibrium state in pyCFRAM 1×1 input format.

    pyCFRAM `run_parallel_python.py` reads NetCDFs with `[0, ::-1, :, :]`,
    so we store fields in surface→TOA order; after the slice they become
    TOA→surface internally.
    """
    nlev = model.lev.size
    plev_top2sfc = np.array(model.lev)            # climlab native: TOA→sfc
    Tatm_top2sfc = np.array(model.Tatm).flatten()
    Ts = float(np.array(model.Ts).flat[0])
    q_top2sfc = np.array(h2o.q).flatten()
    # O3 is part of RRTMG default absorber set. In pure RCE we don't perturb
    # O3, but pyCFRAM still expects an O3 field — write the climlab default.
    o3_top2sfc = np.array(rad.absorber_vmr['O3']).flatten()

    # Reverse to surface→TOA for NetCDF storage
    plev_sfc2top = plev_top2sfc[::-1]
    Tatm_sfc2top = Tatm_top2sfc[::-1]
    q_sfc2top = q_top2sfc[::-1]
    o3_sfc2top = o3_top2sfc[::-1]

    nc = Dataset(out_pres, 'w')
    nc.createDimension('time', 1)
    nc.createDimension('lev', nlev)
    nc.createDimension('lat', 1)
    nc.createDimension('lon', 1)
    nc.createVariable('time', 'f8', ('time',))[:] = [0.0]
    nc.createVariable('lev',  'f8', ('lev',))[:]  = plev_sfc2top
    nc.createVariable('lat',  'f8', ('lat',))[:]  = [0.0]
    nc.createVariable('lon',  'f8', ('lon',))[:]  = [0.0]

    def write3d(name, arr_sfc2top):
        v = nc.createVariable(name, 'f8', ('time', 'lev', 'lat', 'lon'))
        v[:] = np.asarray(arr_sfc2top, dtype=np.float64)[None, :, None, None]

    write3d('ta_lay', Tatm_sfc2top)
    write3d('q',      q_sfc2top)
    write3d('o3',     o3_sfc2top)
    # CO2: pyCFRAM expects 3D vmr field; constant column is fine.
    write3d('co2',    np.full(nlev, co2_ppm * 1e-6))
    # Cloud + aerosols all zero (clear-sky RCE, no aerosol).
    for name in ['camt', 'cliq', 'cice', 'bc', 'ocphi', 'ocpho', 'sulf', 'ss', 'dust']:
        write3d(name, np.zeros(nlev))
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
        v[:] = np.full((1, 1, 1), float(val))

    write2d('ts',     Ts)
    # ps just below the lowest atm plev so the Fortran sees nlayer == nlev.
    # climlab's default 30-level grid puts the bottom layer near 983 hPa, so
    # we set ps = 1000 hPa = 100000 Pa.
    write2d('ps',     100000.0)
    write2d('solar',  INSOLATION)
    write2d('albedo', ALBEDO)
    nc.close()
    print('  -> wrote pyCFRAM input: %s  +  %s' % (out_pres, out_surf))


def main():
    out_dir = os.path.join(HERE, 'output')
    os.makedirs(out_dir, exist_ok=True)

    # Ensure case input dirs exist
    for d in [CASE_DIR_PRIMARY, CASE_DIR_FU]:
        os.makedirs(os.path.join(d, 'input'), exist_ok=True)

    cases = [
        # (label, co2_ppm, max_days, pyCFRAM input filenames)
        ('1xco2', 348.,  1500, ('base_pres.nc',      'base_surf.nc')),
        ('4xco2', 1392., 3000, ('perturbed_pres.nc', 'perturbed_surf.nc')),
    ]
    for label, co2_ppm, max_days, (pres_name, surf_name) in cases:
        print('\n===== Building model: CO2 = %.1f ppm (%s) =====' % (co2_ppm, label))
        model, h2o, rad, _conv = build_model(co2_ppm, name='RCE_' + label)
        elapsed, imbalance = integrate_to_equilibrium(
            model, rad, max_days=max_days, tol=0.02, chunk=20,
        )
        print('  -> equilibrium after %d days (imbalance %.4f W/m²)' %
              (elapsed, imbalance))

        # Native climlab snapshot (for inspection)
        save_state(model, h2o, rad,
                   os.path.join(out_dir, 'rce_%s.nc' % label), co2_ppm)

        # pyCFRAM input — primary case directory
        save_pycfram_input(
            model, h2o, rad, co2_ppm,
            out_pres=os.path.join(CASE_DIR_PRIMARY, 'input', pres_name),
            out_surf=os.path.join(CASE_DIR_PRIMARY, 'input', surf_name),
        )

    # Mirror primary input into the Fu case dir via symlinks (avoid duplication;
    # both cases consume identical climlab output, only the radiation engine
    # differs).
    for fname in ['base_pres.nc', 'base_surf.nc',
                  'perturbed_pres.nc', 'perturbed_surf.nc']:
        src = os.path.relpath(
            os.path.join(CASE_DIR_PRIMARY, 'input', fname),
            os.path.join(CASE_DIR_FU, 'input'),
        )
        dst = os.path.join(CASE_DIR_FU, 'input', fname)
        if os.path.lexists(dst):
            os.remove(dst)
        os.symlink(src, dst)
    print('\nSymlinked input/ → %s/input/' % CASE_DIR_FU)


if __name__ == '__main__':
    main()
