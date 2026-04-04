#!/usr/bin/env python3
"""Download ERA5 monthly mean data for CFRAM-A from CDS API.

Downloads pressure-level and single-level monthly means needed for CFRAM.
Requires ~/.cdsapirc to be configured.

Usage:
    python scripts/download_era5.py --outdir /path/to/era5_data [--years 1979-2020]
"""

import os
import argparse


def download_pressure_levels(outdir, years):
    """Download ERA5 monthly pressure level data."""
    import cdsapi
    c = cdsapi.Client()

    variables_pl = [
        'temperature',
        'specific_humidity',
        'ozone_mass_mixing_ratio',
        'fraction_of_cloud_cover',
        'specific_cloud_liquid_water_content',
        'specific_cloud_ice_water_content',
    ]

    pressure_levels = [
        '1', '2', '3', '5', '7', '10', '20', '30', '50', '70',
        '100', '125', '150', '175', '200', '225', '250', '300',
        '350', '400', '450', '500', '550', '600', '650', '700',
        '750', '775', '800', '825', '850', '875', '900', '925',
        '950', '975', '1000',
    ]

    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, f'era5_pressure_monthly_{years[0]}_{years[-1]}.nc')

    if os.path.exists(outfile):
        print(f"  Already exists: {outfile}")
        return outfile

    print(f"  Downloading pressure level data to {outfile}...")
    c.retrieve(
        'reanalysis-era5-pressure-levels-monthly-means',
        {
            'product_type': 'monthly_averaged_reanalysis',
            'variable': variables_pl,
            'pressure_level': pressure_levels,
            'year': [str(y) for y in years],
            'month': [f'{m:02d}' for m in range(1, 13)],
            'time': '00:00',
            'data_format': 'netcdf',
            'download_format': 'unarchived',
            'area': [42, 95, 18, 125],  # N, W, S, E
        },
        outfile
    )
    return outfile


def download_single_levels(outdir, years):
    """Download ERA5 monthly single level data."""
    import cdsapi
    c = cdsapi.Client()

    variables_sl = [
        'skin_temperature',
        'surface_pressure',
        'surface_solar_radiation_downwards',
        'surface_net_solar_radiation',
        'toa_incident_solar_radiation',
    ]

    os.makedirs(outdir, exist_ok=True)
    outfile = os.path.join(outdir, f'era5_single_monthly_{years[0]}_{years[-1]}.nc')

    if os.path.exists(outfile):
        print(f"  Already exists: {outfile}")
        return outfile

    print(f"  Downloading single level data to {outfile}...")
    c.retrieve(
        'reanalysis-era5-single-levels-monthly-means',
        {
            'product_type': 'monthly_averaged_reanalysis',
            'variable': variables_sl,
            'year': [str(y) for y in years],
            'month': [f'{m:02d}' for m in range(1, 13)],
            'time': '00:00',
            'data_format': 'netcdf',
            'download_format': 'unarchived',
            'area': [42, 95, 18, 125],
        },
        outfile
    )
    return outfile


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Download ERA5 data for CFRAM')
    parser.add_argument('--outdir', default='era5_data', help='Output directory')
    parser.add_argument('--years', default='1979-2020',
                        help='Year range (e.g., 1979-2020)')
    args = parser.parse_args()

    y0, y1 = map(int, args.years.split('-'))
    years = list(range(y0, y1 + 1))

    print("=== ERA5 Pressure Level Download ===")
    download_pressure_levels(args.outdir, years)

    print("=== ERA5 Single Level Download ===")
    download_single_levels(args.outdir, years)

    print("Done. Files saved to:", args.outdir)
