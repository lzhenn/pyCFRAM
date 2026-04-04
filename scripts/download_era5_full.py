#!/usr/bin/env python3
"""Download all ERA5 monthly data needed for CFRAM-A.

Run on hqlx204 (has working SSL + cdsapi):
    nohup python3 -u download_era5_full.py > download_era5.log 2>&1 &

Downloads pressure level + single level monthly means, 1979-2020,
domain [42N, 95E, 18S, 125E], all 37 pressure levels.
"""

import os
import cdsapi

OUTDIR = '/home/lzhenn/work/ust-jumper/pyCFRAM/era5_data'
os.makedirs(OUTDIR, exist_ok=True)

c = cdsapi.Client(
    url='https://cds.climate.copernicus.eu/api',
    key='410ef3ec-11bc-4f8c-b477-e389d5d21b48',
)

PRESSURE_LEVELS = [
    '1', '2', '3', '5', '7', '10', '20', '30', '50', '70',
    '100', '125', '150', '175', '200', '225', '250', '300',
    '350', '400', '450', '500', '550', '600', '650', '700',
    '750', '775', '800', '825', '850', '875', '900', '925',
    '950', '975', '1000',
]

PL_VARIABLES = [
    'temperature',
    'specific_humidity',
    'ozone_mass_mixing_ratio',
    'fraction_of_cloud_cover',
    'specific_cloud_liquid_water_content',
    'specific_cloud_ice_water_content',
]

SL_VARIABLES = [
    'skin_temperature',
    'surface_pressure',
    'surface_solar_radiation_downwards',
    'surface_net_solar_radiation',
    'toa_incident_solar_radiation',
]

MONTHS = [f'{m:02d}' for m in range(1, 13)]
AREA = [42, 95, 18, 125]  # N, W, S, E
YEARS = list(range(1979, 2021))


def download_pl(year):
    """Download pressure level data for one year."""
    outfile = os.path.join(OUTDIR, f'era5_pl_{year}.nc')
    if os.path.exists(outfile):
        print(f'  PL {year}: skip (exists)')
        return
    print(f'  PL {year}: downloading...')
    c.retrieve(
        'reanalysis-era5-pressure-levels-monthly-means',
        {
            'product_type': ['monthly_averaged_reanalysis'],
            'variable': PL_VARIABLES,
            'pressure_level': PRESSURE_LEVELS,
            'year': [str(year)],
            'month': MONTHS,
            'time': ['00:00'],
            'data_format': 'netcdf',
            'download_format': 'unarchived',
            'area': AREA,
        },
        outfile,
    )
    fsize = os.path.getsize(outfile) / 1e6
    print(f'  PL {year}: done ({fsize:.1f} MB)')


def download_sl(year):
    """Download single level data for one year."""
    outfile = os.path.join(OUTDIR, f'era5_sl_{year}.nc')
    if os.path.exists(outfile):
        print(f'  SL {year}: skip (exists)')
        return
    print(f'  SL {year}: downloading...')
    c.retrieve(
        'reanalysis-era5-single-levels-monthly-means',
        {
            'product_type': ['monthly_averaged_reanalysis'],
            'variable': SL_VARIABLES,
            'year': [str(year)],
            'month': MONTHS,
            'time': ['00:00'],
            'data_format': 'netcdf',
            'download_format': 'unarchived',
            'area': AREA,
        },
        outfile,
    )
    fsize = os.path.getsize(outfile) / 1e6
    print(f'  SL {year}: done ({fsize:.1f} MB)')


if __name__ == '__main__':
    print(f'=== ERA5 Download: {len(YEARS)} years, outdir={OUTDIR} ===')
    for year in YEARS:
        try:
            download_pl(year)
        except Exception as e:
            print(f'  PL {year}: ERROR {e}')
        try:
            download_sl(year)
        except Exception as e:
            print(f'  SL {year}: ERROR {e}')
    print('=== All done ===')
