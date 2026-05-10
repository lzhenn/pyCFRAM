#!/usr/bin/env python3
"""Build pyCFRAM input for cesm2_4xco2_official from CMIP6 raw output.

Pipeline:
1. Read piControl (year_start..year_end) climo + abrupt-4xCO2 climo
2. Re-project cl/clw/cli from 32-layer hybrid sigma → 19-layer plev
3. Compute albedo (rsus / rsds, polar-night safe)
4. Write base_pres.nc, base_surf.nc, perturbed_pres.nc, perturbed_surf.nc
5. Write nonrad_forcing.nc (hfls + hfss difference, warm − base)
6. (External) inject CESM 1850 O3 climo via scripts/inject_cesm_o3.py
7. (External) mask subsurface via scripts/mask_subsurface_layers.py

Output NetCDF lev convention: sfc→TOA (matches CMIP6 plev order, pyCFRAM
flips internally with [::-1]). 19 lev: [1000, 925, 850, ..., 5, 1] hPa.
"""
import os
import sys
import argparse
import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case, PROJECT_ROOT
from data.cesm2_cmip6_source import (
    PLEV19_PA, load_climo_pres, hybrid_to_plev_mass_conserving, compute_albedo
)


def build_pres_nc(out_path, climo, plev_pa, co2_ppm, lat, lon):
    """Write pyCFRAM pressure-level input NetCDF.

    All variables stored in NetCDF lev = sfc→TOA order (PLEV19_PA convention).
    """
    nlev = len(plev_pa)
    nlat = len(lat)
    nlon = len(lon)

    nc = Dataset(out_path, 'w')
    nc.createDimension('time', 1)
    nc.createDimension('lev', nlev)
    nc.createDimension('lat', nlat)
    nc.createDimension('lon', nlon)
    nc.createVariable('time', 'f8', ('time',))[:] = [0.0]
    nc.createVariable('lev',  'f8', ('lev',))[:]  = plev_pa / 100.0   # hPa
    nc.createVariable('lat',  'f8', ('lat',))[:]  = lat
    nc.createVariable('lon',  'f8', ('lon',))[:]  = lon

    def write3d(name, arr_sfc2top):
        v = nc.createVariable(name, 'f8', ('time', 'lev', 'lat', 'lon'))
        v[:] = arr_sfc2top[None, :, :, :]

    # ta/hus already in CMIP6 NetCDF order (sfc→TOA); use directly
    write3d('ta_lay', climo['ta'])
    write3d('q',      climo['hus'])
    # cl/clw/cli have been re-projected to plev_pa order (sfc→TOA, matching plev_pa)
    write3d('camt', climo['cl_plev'])
    write3d('cliq', climo['clw_plev'])
    write3d('cice', climo['cli_plev'])
    # O3 placeholder: zeros, to be filled by inject_cesm_o3.py later
    write3d('o3', np.zeros((nlev, nlat, nlon), dtype=np.float64))
    # CO2 as 3D constant, units mol/mol
    write3d('co2', np.full((nlev, nlat, nlon), co2_ppm * 1e-6, dtype=np.float64))
    nc.close()
    print('  Saved %s' % out_path)


def build_surf_nc(out_path, climo, lat, lon):
    """Write pyCFRAM surface input NetCDF."""
    nc = Dataset(out_path, 'w')
    nc.createDimension('time', 1)
    nc.createDimension('lat', len(lat))
    nc.createDimension('lon', len(lon))
    nc.createVariable('time', 'f8', ('time',))[:] = [0.0]
    nc.createVariable('lat',  'f8', ('lat',))[:]  = lat
    nc.createVariable('lon',  'f8', ('lon',))[:]  = lon

    def write2d(name, arr):
        v = nc.createVariable(name, 'f8', ('time', 'lat', 'lon'))
        v[:] = arr[None, :, :]

    write2d('ts',     climo['ts'])
    write2d('ps',     climo['ps'])
    write2d('solar',  climo['rsdt'])
    write2d('albedo', climo['albedo'])
    # huss = 2m specific humidity. Apple-to-apple OLD CFRAM raw/CFRAM.zip
    # GW-base.f L322-330: ph(nv1) = huss(I,J) — surface row of /atmosp/.
    write2d('huss',   climo['huss'])
    nc.close()
    print('  Saved %s' % out_path)


def build_nonrad_forcing(out_path, base, warm, lat, lon, plev_pa):
    """Write nonrad_forcing.nc with hfls/hfss differences.

    Convention (per pyCFRAM run_parallel_python.py):
    - frc_lhflx and frc_shflx are 4D (time, lev, lat, lon)
    - Atm rows (lev != surface): zero
    - Surface row (last lev index in NetCDF order): warm − base difference
    Sign: pyCFRAM uses dT_X = -drdt⁻¹·frc_X. For hfls/hfss to give correct
    dT sign, frc = +(warm − base) (sign flips happen in worker code).
    """
    nlev = len(plev_pa)
    nlat = len(lat)
    nlon = len(lon)

    nc = Dataset(out_path, 'w')
    nc.createDimension('time', 1)
    nc.createDimension('lev', nlev)
    nc.createDimension('lat', nlat)
    nc.createDimension('lon', nlon)
    nc.createVariable('time', 'f8', ('time',))[:] = [0.0]
    nc.createVariable('lev',  'f8', ('lev',))[:]  = plev_pa / 100.0
    nc.createVariable('lat',  'f8', ('lat',))[:]  = lat
    nc.createVariable('lon',  'f8', ('lon',))[:]  = lon

    # Zero atm rows; only surface row (lev[0]=1000 hPa, index 0 in sfc→TOA NetCDF) gets diff
    # NB: pyCFRAM reads with [0, ::-1, :, :] so internal index NLEV (= surface row equivalent)
    # corresponds to NetCDF index 0. Our lev[0]=1000 is sfc-most → NetCDF index 0 holds the
    # surface forcing.
    def write_sfc_only(name, sfc_field):
        v = nc.createVariable(name, 'f8', ('time', 'lev', 'lat', 'lon'))
        arr = np.zeros((nlev, nlat, nlon), dtype=np.float64)
        arr[0, :, :] = sfc_field   # NetCDF index 0 = sfc-most level
        v[:] = arr[None, :, :, :]

    # Sign convention for pyCFRAM nonrad_forcing.nc:
    # OLD CFRAM Fortran: fc_lhflx = +Δhfls, dt_lhflx = drdt⁻¹·fc.
    # pyCFRAM Python: dT = -drdt⁻¹·frc. To match OLD dt_lhflx, store frc = -Δhfls.
    # Same for shflx. (Collaborator's cesm2_4xco2/input/nonrad_forcing.nc follows
    # this convention: lhflx mean = -9.78 W/m² for 4xCO2.)
    delta_hfls = warm['hfls'] - base['hfls']
    delta_hfss = warm['hfss'] - base['hfss']
    write_sfc_only('lhflx', -delta_hfls)
    write_sfc_only('shflx', -delta_hfss)
    nc.close()
    print('  Saved %s  (frc_lhflx=-Δhfls mean=%+.3f, frc_shflx=-Δhfss mean=%+.3f W/m²)' %
          (out_path, -np.mean(delta_hfls), -np.mean(delta_hfss)))


def main():
    parser = argparse.ArgumentParser(description='Build pyCFRAM input from CESM2 CMIP6 raw')
    parser.add_argument('--case', required=True, help='Case name')
    args = parser.parse_args()

    cfg = load_case(args.case)
    src = cfg.get('source', {})
    if src.get('type') != 'cesm2_cmip6':
        sys.exit('ERROR: case %s source.type != cesm2_cmip6' % args.case)

    raw_dir = src['raw_dir']
    pictrl_y = src['pictrl_years']     # [start, end]
    warm_y = src['warm_years']
    co2_pi = src['co2_pictrl_ppm']
    co2_warm = src['co2_warm_ppm']

    case_dir = cfg['_case_dir']
    input_dir = os.path.join(case_dir, 'input')
    os.makedirs(input_dir, exist_ok=True)

    print('=' * 60)
    print('Building pyCFRAM input: %s' % args.case)
    print('  raw_dir: %s' % raw_dir)
    print('  piControl years: %d-%d (CO2=%.1f ppm)' % (pictrl_y[0], pictrl_y[1], co2_pi))
    print('  warm years:      %d-%d (CO2=%.1f ppm)' % (warm_y[0], warm_y[1], co2_warm))
    print('=' * 60)

    # Load both states' annual climos
    print('\n1. Loading piControl climo')
    base = load_climo_pres(raw_dir, src.get('pictrl_subdir', 'piControl'),
                            pictrl_y[0], pictrl_y[1])
    print('\n2. Loading abrupt-4xCO2 climo')
    warm = load_climo_pres(raw_dir, src.get('warm_subdir', 'abrupt-4XCO2'),
                            warm_y[0], warm_y[1])

    lat = base['lat']; lon = base['lon']

    # Re-project cl/clw/cli from 32 hybrid → PLEV19_PA. Use the
    # mass-conserving variant: cumulative-integration → layer means. The old
    # log-linear sampler (`hybrid_to_plev`) lost ~5.5% column liquid water
    # and up to ~24% per-cell in BL clouds (verified by
    # scripts/diag_cloud_column.py).
    print('\n3. Re-projecting clouds hybrid 32 → plev 19 (mass-conserving)')
    a, b, p0 = base['hybrid_a'], base['hybrid_b'], base['hybrid_p0']
    for state, label in [(base, 'piControl'), (warm, '4xCO2')]:
        for var in ('cl', 'clw', 'cli'):
            field_hyb = state[var]   # (32, nlat, nlon) TOA→sfc
            ps_2d = state['ps']
            field_plev_top_down = hybrid_to_plev_mass_conserving(
                field_hyb, a, b, p0, ps_2d, PLEV19_PA[::-1])
            # PLEV19_PA[::-1] is TOA→sfc in Pa; result top_down is TOA→sfc.
            # Flip back to sfc→TOA to match PLEV19_PA NetCDF order:
            state[var + '_plev'] = field_plev_top_down[::-1]
        print('  %s: cl_plev shape=%s, range [%.3f, %.3f]' %
              (label, state['cl_plev'].shape, state['cl_plev'].min(), state['cl_plev'].max()))

    # Albedo
    print('\n4. Computing albedo')
    base['albedo'] = compute_albedo(base['rsus'], base['rsds'])
    warm['albedo'] = compute_albedo(warm['rsus'], warm['rsds'])
    print('  base albedo mean=%.3f, warm=%.3f' %
          (np.mean(base['albedo']), np.mean(warm['albedo'])))

    # Write 4 input NCs
    print('\n5. Writing pyCFRAM input NetCDFs')
    build_pres_nc(os.path.join(input_dir, 'base_pres.nc'),      base, PLEV19_PA, co2_pi, lat, lon)
    build_pres_nc(os.path.join(input_dir, 'perturbed_pres.nc'), warm, PLEV19_PA, co2_warm, lat, lon)
    build_surf_nc(os.path.join(input_dir, 'base_surf.nc'),      base, lat, lon)
    build_surf_nc(os.path.join(input_dir, 'perturbed_surf.nc'), warm, lat, lon)
    build_nonrad_forcing(os.path.join(input_dir, 'nonrad_forcing.nc'),
                         base, warm, lat, lon, PLEV19_PA)

    print('\nDone. Next steps:')
    print('  python3 scripts/inject_cesm_o3.py %s     # inject CESM 1850 O3 climo' % args.case)
    print('  python3 scripts/mask_subsurface_layers.py %s   # mask sub-ps levels' % args.case)
    print('  python3 run_case.py %s --step run --nproc 384' % args.case)


if __name__ == '__main__':
    main()
