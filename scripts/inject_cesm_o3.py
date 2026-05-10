#!/usr/bin/env python3
"""Replace cesm2_4xco2 input O3 with CESM2 1850 climatology O3.

Source: /tmp/ozone_1.9x2.5_L26_1850clim_c090420.nc (from hqlx86 metctm1).
- Annual + zonal mean → (26 lev, 96 lat) profile.
- Convert mol/mol → kg/kg (×48/29).
- Interpolate vertically to 37-level pyCFRAM grid (1000-1 hPa).
- Interpolate latitudinally to 192-point cesm2 grid.
- Broadcast over 288 longitudes.
- Apply identically to base_pres.nc AND perturbed_pres.nc (no frc_o3 perturbation).

Backup: base_pres.nc.bak_o3eh13 / perturbed_pres.nc.bak_o3eh13 (the eh13-derived O3).
"""
import os
import sys
import shutil
import numpy as np
from netCDF4 import Dataset

# Default O3 source path on hqlx220/hqlx74 (CESM 1850 climo from hqlx86 metctm1)
CESM_O3_CANDIDATES = [
    '/tmp/ozone_1.9x2.5_L26_1850clim_c090420.nc',
    '/home/lzhenn/work/ust-jumper/pyCFRAM/raw_data/ozone_1.9x2.5_L26_1850clim_c090420.nc',
    '/Users/zhenningli/work/ust-jumper/pyCFRAM/raw_data/ozone_1.9x2.5_L26_1850clim_c090420.nc',
]


def find_o3_source():
    for p in CESM_O3_CANDIDATES:
        if os.path.exists(p):
            return p
    raise FileNotFoundError(
        'CESM 1850 O3 climo not found. Place at one of: %s' %
        ', '.join(CESM_O3_CANDIDATES))


CESM_O3 = None  # set in main() once case CWD known

# mol/mol → kg/kg (M_O3 / M_air = 48 / 29)
VMR_TO_MMR = 48.0 / 29.0


sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case


def load_cesm_o3():
    """Load CESM 1850 climatology O3, return zonal+annual mean profile.

    Returns
    -------
    o3_zm : (nlev_src, nlat_src) array, kg/kg
    lev_src : (nlev_src,) hPa, ascending
    lat_src : (nlat_src,) degrees, ascending
    """
    nc = Dataset(CESM_O3)
    o3_vmr = np.array(nc.variables['O3'][:])    # (12, 26, 96, 144) mol/mol
    lev = np.array(nc.variables['lev'][:])       # 26 levels, ascending pressure
    lat = np.array(nc.variables['lat'][:])       # 96 lat, ascending
    nc.close()

    # Annual mean + zonal mean
    o3_zm_vmr = np.nanmean(o3_vmr, axis=(0, 3))  # (26, 96)
    o3_zm_mmr = o3_zm_vmr * VMR_TO_MMR
    return o3_zm_mmr, lev, lat


def interp_to_pycfram_grid(o3_zm_src, lev_src, lat_src, lev_tgt, lat_tgt):
    """Interpolate (lev_src, lat_src) → (lev_tgt, lat_tgt) by 1D linear in log(p) and lat."""
    # Vertical interpolation in log(p)
    logp_src = np.log(lev_src)
    logp_tgt = np.log(lev_tgt)

    # Step 1: vertical interp at each src latitude
    o3_v = np.empty((len(lev_tgt), len(lat_src)))
    for j in range(len(lat_src)):
        # extrapolate via clip to ends (ozone is small at top, smaller at bottom)
        o3_v[:, j] = np.interp(logp_tgt, logp_src, o3_zm_src[:, j],
                               left=o3_zm_src[0, j], right=o3_zm_src[-1, j])

    # Step 2: latitudinal interp at each tgt level
    o3_out = np.empty((len(lev_tgt), len(lat_tgt)))
    for k in range(len(lev_tgt)):
        o3_out[k, :] = np.interp(lat_tgt, lat_src, o3_v[k, :])
    return o3_out


def inject(target_path, o3_3d):
    """Write o3_3d (lev, lat, lon) into target_path's o3 variable in place."""
    nc = Dataset(target_path, 'r+')
    nc.variables['o3'][0, :, :, :] = o3_3d
    if hasattr(nc, 'o3_source'):
        del nc.o3_source
    nc.o3_source = 'CESM2 1850 climatology (zonal+annual mean), inject_cesm_o3.py'
    nc.close()


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', required=True)
    args = ap.parse_args()
    cfg = load_case(args.case)
    case_input = os.path.join(cfg['_case_dir'], 'input')

    global CESM_O3
    CESM_O3 = find_o3_source()
    print('Using O3 source: %s' % CESM_O3)
    print('Loading CESM2 1850 climatology O3...')
    o3_src, lev_src, lat_src = load_cesm_o3()
    print('  src: lev %.1f-%.1f hPa (%d lev), lat %.1f-%.1f (%d lat)' %
          (lev_src.min(), lev_src.max(), len(lev_src),
           lat_src.min(), lat_src.max(), len(lat_src)))
    print('  global O3 range: %.2e - %.2e kg/kg' % (o3_src.min(), o3_src.max()))

    # Read pyCFRAM target grid
    nc = Dataset(os.path.join(case_input, 'base_pres.nc'))
    lev_tgt = np.array(nc.variables['lev'][:])
    lat_tgt = np.array(nc.variables['lat'][:])
    lon_tgt = np.array(nc.variables['lon'][:])
    nc.close()
    nlev, nlat, nlon = len(lev_tgt), len(lat_tgt), len(lon_tgt)
    print('  tgt: lev %.1f-%.1f hPa (%d lev), lat %.1f-%.1f (%d lat), %d lon' %
          (lev_tgt.min(), lev_tgt.max(), nlev,
           lat_tgt.min(), lat_tgt.max(), nlat, nlon))

    # Interpolate to (37, 192)
    o3_zm = interp_to_pycfram_grid(o3_src, lev_src, lat_src, lev_tgt, lat_tgt)
    print('  interpolated zonal-mean O3 range: %.2e - %.2e kg/kg' %
          (o3_zm.min(), o3_zm.max()))

    # Sample at 30 hPa (peak): equator vs poles
    k30 = np.abs(lev_tgt - 30).argmin()
    print('  O3 at %d hPa: equator=%.2e, 60N=%.2e, 60S=%.2e kg/kg' %
          (lev_tgt[k30], o3_zm[k30, np.abs(lat_tgt-0).argmin()],
           o3_zm[k30, np.abs(lat_tgt-60).argmin()],
           o3_zm[k30, np.abs(lat_tgt+60).argmin()]))

    # Broadcast over lon
    o3_3d = np.broadcast_to(o3_zm[:, :, None], (nlev, nlat, nlon)).copy()

    # Backup + write to both base and perturbed (no frc_o3 perturbation)
    for fname in ['base_pres.nc', 'perturbed_pres.nc']:
        path = os.path.join(case_input, fname)
        bak = path + '.bak_o3eh13'
        if not os.path.exists(bak):
            shutil.copy(path, bak)
            print('  backup: %s' % bak)
        inject(path, o3_3d)
        print('  injected CESM2 1850 climatology O3 into %s' % fname)


if __name__ == '__main__':
    main()
