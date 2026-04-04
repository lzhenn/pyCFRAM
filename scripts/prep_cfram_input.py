#!/usr/bin/env python3
"""Generate CFRAM-A Fortran binary input files from ERA5 NetCDF data.

Replaces all NCL preprocessing scripts (t_3d.ncl, q_3d.ncl, etc.).
Reads ERA5 monthly data, computes base/warm period averages, and writes
binary files compatible with cfram_rrtmg.f90.

Usage:
    python scripts/prep_cfram_input.py \
        --era5_pl  era5_pressure_monthly.nc \
        --era5_sl  era5_single_monthly.nc \
        --outdir   CFRAM_RRTMG/data_prep \
        --base_years 1986-1995 \
        --warm_years 2004-2013 \
        --lat_idx 0 --lon_idx 0
"""

import os
import argparse
import numpy as np
from netCDF4 import Dataset, num2date


def get_year_month(nc, time_var='time'):
    """Extract year and month arrays from NetCDF time variable."""
    t = nc.variables[time_var]
    dates = num2date(t[:], t.units, t.calendar if hasattr(t, 'calendar') else 'standard')
    years = np.array([d.year for d in dates])
    months = np.array([d.month for d in dates])
    return years, months


def monthly_to_annual_mean(data, years, months):
    """Convert monthly data to annual means.

    Args:
        data: array with time as first dimension
        years: year for each time step
        months: month for each time step

    Returns:
        annual_means: array with year as first dimension
        unique_years: array of years
    """
    unique_years = np.unique(years)
    shape = (len(unique_years),) + data.shape[1:]
    annual = np.zeros(shape, dtype=np.float64)

    for i, yr in enumerate(unique_years):
        mask = years == yr
        annual[i] = np.mean(data[mask], axis=0)

    return annual, unique_years


def period_mean(annual_data, annual_years, yr_start, yr_end):
    """Compute mean over a period of years."""
    mask = (annual_years >= yr_start) & (annual_years <= yr_end)
    return np.mean(annual_data[mask], axis=0)


def write_bin(filepath, data):
    """Write numpy array to Fortran-compatible binary (raw float64)."""
    np.asarray(data, dtype=np.float64).tofile(filepath)
    print(f"  Wrote {filepath} ({data.size} values, {data.nbytes} bytes)")


def prep_from_era5(era5_pl_file, era5_sl_file, outdir,
                   base_years, warm_years, lat_idx=0, lon_idx=0,
                   co2_file=None):
    """Generate all CFRAM binary input files from ERA5 NetCDF.

    Args:
        era5_pl_file: ERA5 pressure level monthly means NetCDF
        era5_sl_file: ERA5 single level monthly means NetCDF
        outdir: output directory for .dat files
        base_years: (start, end) tuple for base period
        warm_years: (start, end) tuple for warm period
        lat_idx: latitude index to extract (single column)
        lon_idx: longitude index to extract (single column)
        co2_file: path to co2.txt (year, ppmv); if None, uses defaults
    """
    os.makedirs(outdir, exist_ok=True)

    # --- Pressure level variables ---
    print(f"Reading pressure level data from {era5_pl_file}")
    nc_pl = Dataset(era5_pl_file, 'r')
    years_pl, months_pl = get_year_month(nc_pl)

    # Variable name mapping (ERA5 NetCDF name -> CFRAM name)
    # Adjust variable names based on actual NetCDF content
    var_names_3d = {}
    for vname in nc_pl.variables:
        vlow = vname.lower()
        if vlow in ('t', 'temperature'):
            var_names_3d['t'] = vname
        elif vlow in ('q', 'specific_humidity'):
            var_names_3d['q'] = vname
        elif vlow in ('o3', 'ozone_mass_mixing_ratio', 'ozone'):
            var_names_3d['o3'] = vname
        elif vlow in ('cc', 'fraction_of_cloud_cover'):
            var_names_3d['cc'] = vname
        elif vlow in ('clwc', 'specific_cloud_liquid_water_content'):
            var_names_3d['clwc'] = vname
        elif vlow in ('ciwc', 'specific_cloud_ice_water_content'):
            var_names_3d['ciwc'] = vname

    print(f"  Found 3D variables: {list(var_names_3d.keys())}")
    print(f"  Extracting lat_idx={lat_idx}, lon_idx={lon_idx}")

    for short_name, nc_name in var_names_3d.items():
        print(f"  Processing {short_name} ({nc_name})...")
        # Read: (time, level, lat, lon) -> extract single column
        raw = np.array(nc_pl.variables[nc_name][:, :, lat_idx, lon_idx],
                       dtype=np.float64)
        # Monthly -> annual mean
        annual, ann_years = monthly_to_annual_mean(raw, years_pl, months_pl)
        # Period means
        base_mean = period_mean(annual, ann_years, *base_years)
        warm_mean = period_mean(annual, ann_years, *warm_years)

        # Output file names matching Fortran expectations
        name_map = {
            't': ('t_base.dat', 't_warm.dat'),
            'q': ('hus_base.dat', 'hus_warm.dat'),
            'o3': ('O3_base.dat', 'O3_warm.dat'),
            'cc': ('cc_base.dat', 'cc_warm.dat'),
            'clwc': ('clwc_base.dat', 'clwc_warm.dat'),
            'ciwc': ('ciwc_base.dat', 'ciwc_warm.dat'),
        }
        base_f, warm_f = name_map[short_name]
        write_bin(os.path.join(outdir, base_f), base_mean)
        write_bin(os.path.join(outdir, warm_f), warm_mean)

    nc_pl.close()

    # --- Single level variables ---
    print(f"\nReading single level data from {era5_sl_file}")
    nc_sl = Dataset(era5_sl_file, 'r')
    years_sl, months_sl = get_year_month(nc_sl)

    # Detect variable names
    sl_vars = {}
    for vname in nc_sl.variables:
        vlow = vname.lower()
        if vlow in ('skt', 'skin_temperature'):
            sl_vars['skt'] = vname
        elif vlow in ('sp', 'surface_pressure'):
            sl_vars['sp'] = vname
        elif vlow in ('ssrd', 'surface_solar_radiation_downwards'):
            sl_vars['ssrd'] = vname
        elif vlow in ('ssr', 'surface_net_solar_radiation'):
            sl_vars['ssr'] = vname
        elif vlow in ('tisr', 'toa_incident_solar_radiation'):
            sl_vars['tisr'] = vname

    print(f"  Found SL variables: {list(sl_vars.keys())}")

    for short_name, nc_name in sl_vars.items():
        print(f"  Processing {short_name} ({nc_name})...")
        raw = np.array(nc_sl.variables[nc_name][:, lat_idx, lon_idx],
                       dtype=np.float64)

        # Unit conversion: ERA5 accumulations J/m2 -> W/m2 (divide by seconds in a day)
        if short_name in ('ssrd', 'ssr', 'tisr'):
            raw = raw / 86400.0

        annual, ann_years = monthly_to_annual_mean(raw, years_sl, months_sl)
        base_val = period_mean(annual, ann_years, *base_years)
        warm_val = period_mean(annual, ann_years, *warm_years)

        name_map = {
            'skt': ('skt_base.dat', 'skt_warm.dat'),
            'sp': ('sp_base.dat', 'sp_warm.dat'),
            'ssrd': ('ssrd_base.dat', 'ssrd_warm.dat'),
            'tisr': ('solarin_base.dat', 'solarin_warm.dat'),
        }

        if short_name == 'ssr':
            # Compute ssru = ssrd - ssr (upwelling = downwelling - net)
            # Need ssrd to be processed first
            pass  # handled below
        elif short_name in name_map:
            base_f, warm_f = name_map[short_name]
            # Scalar for single column
            write_bin(os.path.join(outdir, base_f), np.array([base_val]))
            write_bin(os.path.join(outdir, warm_f), np.array([warm_val]))

    # Compute ssru = ssrd - ssr
    if 'ssrd' in sl_vars and 'ssr' in sl_vars:
        print("  Computing ssru = ssrd - ssr...")
        ssrd_raw = np.array(nc_sl.variables[sl_vars['ssrd']][:, lat_idx, lon_idx],
                            dtype=np.float64) / 86400.0
        ssr_raw = np.array(nc_sl.variables[sl_vars['ssr']][:, lat_idx, lon_idx],
                           dtype=np.float64) / 86400.0
        ssru_raw = ssrd_raw - ssr_raw

        annual_ssru, ann_years = monthly_to_annual_mean(ssru_raw, years_sl, months_sl)
        ssru_base = period_mean(annual_ssru, ann_years, *base_years)
        ssru_warm = period_mean(annual_ssru, ann_years, *warm_years)
        write_bin(os.path.join(outdir, 'ssru_base.dat'), np.array([ssru_base]))
        write_bin(os.path.join(outdir, 'ssru_warm.dat'), np.array([ssru_warm]))

    nc_sl.close()

    # --- CO2 ---
    print("\nProcessing CO2...")
    if co2_file and os.path.exists(co2_file):
        co2_data = np.loadtxt(co2_file)  # year, ppmv
        co2_years = co2_data[:, 0].astype(int)
        co2_ppmv = co2_data[:, 1]
        mask_base = (co2_years >= base_years[0]) & (co2_years <= base_years[1])
        mask_warm = (co2_years >= warm_years[0]) & (co2_years <= warm_years[1])
        co2_base = np.mean(co2_ppmv[mask_base])
        co2_warm = np.mean(co2_ppmv[mask_warm])
    else:
        # Default values (approximate)
        co2_base = 352.0  # ~1986-1995 mean
        co2_warm = 385.0  # ~2004-2013 mean
        print(f"  Using default CO2: base={co2_base}, warm={co2_warm} ppmv")

    write_bin(os.path.join(outdir, 'co2_b.dat'), np.array([co2_base]))
    write_bin(os.path.join(outdir, 'co2_w.dat'), np.array([co2_warm]))

    # --- Aerosol (placeholder zeros if no MERRA-2 data) ---
    print("\nAerosol optical properties...")
    nlev = 37
    nbndlw = 16
    jpband = 14
    zeros_lw = np.zeros((nlev, 1, 1, nbndlw), dtype=np.float64)
    zeros_sw = np.zeros((nlev, 1, 1, jpband), dtype=np.float64)

    for suffix in ('base', 'warm'):
        write_bin(os.path.join(outdir, f'aerosol_aod_lw_{suffix}.dat'), zeros_lw)
        write_bin(os.path.join(outdir, f'aerosol_aod_sw_{suffix}.dat'), zeros_sw)
        write_bin(os.path.join(outdir, f'aerosol_ssa_sw_{suffix}.dat'), zeros_sw)
        write_bin(os.path.join(outdir, f'aerosol_g_sw_{suffix}.dat'), zeros_sw)
    print("  (Using zero aerosol - replace with MERRA-2 data later)")

    print(f"\n=== Done. Output files in {outdir} ===")
    print(f"Total files: {len([f for f in os.listdir(outdir) if f.endswith('.dat')])}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Generate CFRAM binary input from ERA5 NetCDF')
    parser.add_argument('--era5_pl', required=True,
                        help='ERA5 pressure level monthly means NetCDF')
    parser.add_argument('--era5_sl', required=True,
                        help='ERA5 single level monthly means NetCDF')
    parser.add_argument('--outdir', default='data_prep',
                        help='Output directory for binary files')
    parser.add_argument('--base_years', default='1986-1995',
                        help='Base period (e.g., 1986-1995)')
    parser.add_argument('--warm_years', default='2004-2013',
                        help='Warm period (e.g., 2004-2013)')
    parser.add_argument('--lat_idx', type=int, default=0,
                        help='Latitude index for single-column extraction')
    parser.add_argument('--lon_idx', type=int, default=0,
                        help='Longitude index for single-column extraction')
    parser.add_argument('--co2_file', default=None,
                        help='Path to co2.txt (year ppmv)')
    args = parser.parse_args()

    base = tuple(map(int, args.base_years.split('-')))
    warm = tuple(map(int, args.warm_years.split('-')))

    prep_from_era5(
        era5_pl_file=args.era5_pl,
        era5_sl_file=args.era5_sl,
        outdir=args.outdir,
        base_years=base,
        warm_years=warm,
        lat_idx=args.lat_idx,
        lon_idx=args.lon_idx,
        co2_file=args.co2_file,
    )
