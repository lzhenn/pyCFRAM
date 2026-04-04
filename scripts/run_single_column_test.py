#!/usr/bin/env python3
"""Phase 2: Single-column test with real ERA5 data.

Extracts one grid point (default: 115E, 32N in EH13 key area),
computes base/warm states from ERA5 daily data, writes Fortran
binary input, and runs cfram_rrtmg.

Run on hqlx204 (has netCDF4):
    python3 scripts/run_single_column_test.py
"""

import os
import sys
import numpy as np
from netCDF4 import Dataset

# Config
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(PROJECT_ROOT, "era5_data", "daily")
FORTRAN_DIR = os.path.join(PROJECT_ROOT, "fortran")
YEARS = range(2003, 2023)

# EH13 warm period: Aug 5-17 (0-based: days 4-16)
EH13_WARM_DAYS = list(range(4, 17))  # 13 days
# Target grid point: ~115E, 32N (within EH13 key area)
TARGET_LAT = 32.0
TARGET_LON = 115.0


def load_pl_all_years(short_name, years, day_indices, lat_idx, lon_idx):
    """Load PL variable for all years, extract warm-period mean at one point."""
    all_means = []
    for yr in years:
        fname = os.path.join(DATA_DIR, "era5_pl_%s_%d08.nc" % (short_name, yr))
        nc = Dataset(fname, 'r')
        # Find data variable
        skip = {'number', 'valid_time', 'pressure_level', 'latitude', 'longitude', 'expver', 'time', 'lat', 'lon', 'level'}
        var_name = [v for v in nc.variables if v not in skip][0]
        # (ntimes, nlev, nlat, nlon) - 6-hourly
        raw = np.array(nc.variables[var_name][:], dtype=np.float64)
        nc.close()
        # 6-hourly -> daily mean: (31, 4, nlev, nlat, nlon) -> (31, nlev, nlat, nlon)
        ndays = 31
        ntpd = raw.shape[0] // ndays
        daily = raw.reshape(ndays, ntpd, *raw.shape[1:]).mean(axis=1)
        # Warm period mean at target point
        warm_mean = daily[day_indices, :, lat_idx, lon_idx].mean(axis=0)
        all_means.append(warm_mean)
    return np.array(all_means)  # (nyears, nlev)


def load_sl_all_years(var_name, step_type, years, day_indices, lat_idx, lon_idx):
    """Load SL variable for all years, extract warm-period mean at one point."""
    all_means = []
    for yr in years:
        sl_dir = os.path.join(DATA_DIR, "era5_sl_%d08" % yr)
        fname = os.path.join(sl_dir, "data_stream-oper_stepType-%s.nc" % step_type)
        nc = Dataset(fname, 'r')
        raw = np.array(nc.variables[var_name][:], dtype=np.float64)
        nc.close()
        ndays = 31
        ntpd = raw.shape[0] // ndays
        daily = raw.reshape(ndays, ntpd, *raw.shape[1:]).mean(axis=1)
        warm_mean = daily[day_indices, lat_idx, lon_idx].mean()
        all_means.append(warm_mean)
    return np.array(all_means)  # (nyears,)


def find_nearest_idx(arr, val):
    return int(np.argmin(np.abs(arr - val)))


def write_bin(filepath, data):
    np.asarray(data, dtype=np.float64).tofile(filepath)


def main():
    print("=== Phase 2: Single-column test with real ERA5 ===")

    # Get lat/lon/level arrays from first file
    nc = Dataset(os.path.join(DATA_DIR, "era5_pl_t_200308.nc"), 'r')
    lats = np.array(nc.variables['latitude'][:])
    lons = np.array(nc.variables['longitude'][:])
    levels = np.array(nc.variables['pressure_level'][:])
    nc.close()

    lat_idx = find_nearest_idx(lats, TARGET_LAT)
    lon_idx = find_nearest_idx(lons, TARGET_LON)
    print("Target: %.1fE, %.1fN -> idx (%d, %d) = (%.2fN, %.2fE)" %
          (TARGET_LON, TARGET_LAT, lat_idx, lon_idx, lats[lat_idx], lons[lon_idx]))
    print("Levels: %d, range %s-%s hPa" % (len(levels), levels.min(), levels.max()))
    print("Warm days: Aug %d-%d (0-based: %s)" %
          (EH13_WARM_DAYS[0]+1, EH13_WARM_DAYS[-1]+1, EH13_WARM_DAYS))

    # Load all PL variables
    print("\nLoading PL variables (20 years x 6 vars)...")
    pl_vars = {'t': 't', 'q': 'q', 'o3': 'o3', 'cc': 'cc', 'clwc': 'clwc', 'ciwc': 'ciwc'}
    pl_data = {}
    for short in pl_vars:
        print("  %s..." % short, end=' ', flush=True)
        data = load_pl_all_years(short, YEARS, EH13_WARM_DAYS, lat_idx, lon_idx)
        pl_data[short] = data
        print("shape=%s, range=[%.4e, %.4e]" % (data.shape, data.min(), data.max()))

    # Climatology (mean over years) and event year (2013, index 10)
    event_idx = list(YEARS).index(2013)
    print("\nEvent year index: %d (year %d)" % (event_idx, 2013))

    # Load SL variables
    print("\nLoading SL variables...")
    sl_spec = {
        'skt': ('skt', 'instant'),
        'sp': ('sp', 'instant'),
        'ssrd': ('ssrd', 'accum'),
        'ssr': ('ssr', 'accum'),
        'tisr': ('tisr', 'accum'),
    }
    sl_data = {}
    for key, (var, step) in sl_spec.items():
        print("  %s..." % key, end=' ', flush=True)
        data = load_sl_all_years(var, step, YEARS, EH13_WARM_DAYS, lat_idx, lon_idx)
        sl_data[key] = data
        print("shape=%s, range=[%.4f, %.4f]" % (data.shape, data.min(), data.max()))

    # CO2 from co2.txt
    co2_file = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                            "OA/202604020900/CFRAM_RRTMG/data_prep/co2.txt")
    if os.path.exists(co2_file):
        co2_data = np.loadtxt(co2_file)
        co2_years = co2_data[:, 0].astype(int)
        co2_ppmv = co2_data[:, 1]
    else:
        # Default approximate values
        co2_years = np.arange(2003, 2023)
        co2_ppmv = 380 + (co2_years - 2003) * 2.2  # rough approximation
        print("  Warning: co2.txt not found, using approximation")

    # Compute base and warm profiles
    print("\n=== Computing base/warm states ===")
    outdir = os.path.join(FORTRAN_DIR, "data_prep")
    os.makedirs(outdir, exist_ok=True)

    # PL variables
    for short, fname_prefix in [('t', 't'), ('q', 'hus'), ('o3', 'O3'),
                                 ('cc', 'cc'), ('clwc', 'clwc'), ('ciwc', 'ciwc')]:
        base = pl_data[short].mean(axis=0)  # 20-year mean
        warm = pl_data[short][event_idx]     # 2013
        write_bin(os.path.join(outdir, "%s_base.dat" % fname_prefix), base)
        write_bin(os.path.join(outdir, "%s_warm.dat" % fname_prefix), warm)
        print("  %s: base=[%.4f,%.4f], warm=[%.4f,%.4f]" %
              (fname_prefix, base.min(), base.max(), warm.min(), warm.max()))

    # SL variables
    for key, fname in [('skt', 'skt'), ('sp', 'sp')]:
        base = sl_data[key].mean()
        warm = sl_data[key][event_idx]
        write_bin(os.path.join(outdir, "%s_base.dat" % fname), np.array([base]))
        write_bin(os.path.join(outdir, "%s_warm.dat" % fname), np.array([warm]))
        print("  %s: base=%.4f, warm=%.4f" % (fname, base, warm))

    # Solar radiation (needs unit conversion if accumulated)
    for key, fname in [('tisr', 'solarin'), ('ssrd', 'ssrd')]:
        base = sl_data[key].mean()
        warm = sl_data[key][event_idx]
        # ERA5 accumulated fields: J/m² per time step → W/m²
        # For 6-hourly data averaged to daily, divide by seconds in accumulation period
        # Actually for daily mean of 6h accumulations, the values may already represent rates
        # Check magnitude: if > 10000, likely J/m² needing /86400
        if abs(base) > 10000:
            base /= 86400.0
            warm /= 86400.0
        write_bin(os.path.join(outdir, "%s_base.dat" % fname), np.array([base]))
        write_bin(os.path.join(outdir, "%s_warm.dat" % fname), np.array([warm]))
        print("  %s: base=%.4f W/m2, warm=%.4f W/m2" % (fname, base, warm))

    # ssru = ssrd - ssr
    for period, idx_or_mean in [('base', None), ('warm', event_idx)]:
        if idx_or_mean is None:
            ssrd_val = sl_data['ssrd'].mean()
            ssr_val = sl_data['ssr'].mean()
        else:
            ssrd_val = sl_data['ssrd'][idx_or_mean]
            ssr_val = sl_data['ssr'][idx_or_mean]
        if abs(ssrd_val) > 10000:
            ssrd_val /= 86400.0
            ssr_val /= 86400.0
        ssru = ssrd_val - ssr_val
        write_bin(os.path.join(outdir, "ssru_%s.dat" % period), np.array([ssru]))
        print("  ssru_%s: %.4f W/m2" % (period, ssru))

    # CO2
    co2_mask = np.isin(co2_years, list(YEARS))
    if co2_mask.any():
        co2_base = co2_ppmv[co2_mask].mean()
    else:
        co2_base = 390.0
    co2_yr_mask = co2_years == 2013
    co2_warm = co2_ppmv[co2_yr_mask][0] if co2_yr_mask.any() else 396.5
    write_bin(os.path.join(outdir, "co2_b.dat"), np.array([co2_base]))
    write_bin(os.path.join(outdir, "co2_w.dat"), np.array([co2_warm]))
    print("  co2: base=%.2f ppmv, warm=%.2f ppmv" % (co2_base, co2_warm))

    # Aerosol (zero for now - no MERRA-2)
    nlev = len(levels)
    zeros_lw = np.zeros((nlev, 1, 1, 16))
    zeros_sw = np.zeros((nlev, 1, 1, 14))
    for suffix in ('base', 'warm'):
        write_bin(os.path.join(outdir, "aerosol_aod_lw_%s.dat" % suffix), zeros_lw)
        write_bin(os.path.join(outdir, "aerosol_aod_sw_%s.dat" % suffix), zeros_sw)
        write_bin(os.path.join(outdir, "aerosol_ssa_sw_%s.dat" % suffix), zeros_sw)
        write_bin(os.path.join(outdir, "aerosol_g_sw_%s.dat" % suffix), zeros_sw)
    print("  aerosol: zeros (no MERRA-2 yet)")

    print("\n=== All input files written to %s ===" % outdir)
    print("Files: %d" % len([f for f in os.listdir(outdir) if f.endswith('.dat')]))
    print("\nTo run Fortran CFRAM on hqlx74:")
    print("  cd %s && ./cfram_rrtmg" % FORTRAN_DIR)


if __name__ == '__main__':
    main()
