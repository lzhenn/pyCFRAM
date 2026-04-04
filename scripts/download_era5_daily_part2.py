#!/usr/bin/env python3
"""Download ERA5 daily data - PART 2 (2020-2022).
Run on hqlx205 to parallelize with hqlx204 (doing 2017-2019).

    nohup python3 -u scripts/download_era5_daily_part2.py > era5_daily_download_p2.log 2>&1 &
"""

import os
import cdsapi

OUTDIR = '/home/lzhenn/work/ust-jumper/pyCFRAM/era5_data/daily'
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
    'maximum_2m_temperature',
]

PL_SHORT = {
    'temperature': 't',
    'specific_humidity': 'q',
    'ozone_mass_mixing_ratio': 'o3',
    'fraction_of_cloud_cover': 'cc',
    'specific_cloud_liquid_water_content': 'clwc',
    'specific_cloud_ice_water_content': 'ciwc',
}

DAYS = [f'{d:02d}' for d in range(1, 32)]
TIMES = ['00:00', '06:00', '12:00', '18:00']
AREA = [42, 95, 18, 125]
YEARS = [2020, 2021, 2022]  # Part 2 only


def download_pl_var(year, variable):
    short = PL_SHORT.get(variable, variable[:4])
    outfile = os.path.join(OUTDIR, f'era5_pl_{short}_{year}08.nc')
    if os.path.exists(outfile):
        print(f'  PL {short} {year}: skip', flush=True)
        return
    print(f'  PL {short} {year}08: downloading...', flush=True)
    c.retrieve(
        'reanalysis-era5-pressure-levels',
        {
            'product_type': ['reanalysis'],
            'variable': [variable],
            'pressure_level': PRESSURE_LEVELS,
            'year': [str(year)],
            'month': ['08'],
            'day': DAYS,
            'time': TIMES,
            'data_format': 'netcdf',
            'download_format': 'unarchived',
            'area': AREA,
        },
        outfile,
    )
    fsize = os.path.getsize(outfile) / 1e6
    print(f'  PL {short} {year}08: done ({fsize:.1f} MB)', flush=True)


def download_sl(year):
    outfile = os.path.join(OUTDIR, f'era5_sl_{year}08.nc')
    if os.path.exists(outfile):
        print(f'  SL {year}: skip', flush=True)
        return
    print(f'  SL {year}08: downloading...', flush=True)
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': ['reanalysis'],
            'variable': SL_VARIABLES,
            'year': [str(year)],
            'month': ['08'],
            'day': DAYS,
            'time': TIMES,
            'data_format': 'netcdf',
            'download_format': 'unarchived',
            'area': AREA,
        },
        outfile,
    )
    fsize = os.path.getsize(outfile) / 1e6
    print(f'  SL {year}08: done ({fsize:.1f} MB)', flush=True)


if __name__ == '__main__':
    print(f'=== ERA5 Daily Part2: Aug {YEARS} ===', flush=True)
    for year in YEARS:
        for var in PL_VARIABLES:
            try:
                download_pl_var(year, var)
            except Exception as e:
                print(f'  PL {PL_SHORT.get(var,var)} {year}: ERROR {e}', flush=True)
        try:
            download_sl(year)
        except Exception as e:
            print(f'  SL {year}: ERROR {e}', flush=True)
    print('=== Part2 done ===', flush=True)
