#!/usr/bin/env python3
"""Download ERA5 surface heat fluxes (slhf, sshf) for pyCFRAM nonrad forcing.

These are 6-hourly accumulated fields (J/m²), same stepType as ssrd/ssr/tisr.

Run on hqlx204:
    nohup python3 -u scripts/download_era5_flux.py > /tmp/download_flux.log 2>&1 &

Output:
    era5_data/daily/era5_sl_{YYYY}08/data_stream-oper_stepType-accum_flux.nc
"""

import os
import sys
import cdsapi

OUTDIR_BASE = 'era5_data/daily'
AREA = [42, 95, 18, 125]  # N, W, S, E

VARIABLES = [
    'surface_latent_heat_flux',     # slhf (J/m², accum)
    'surface_sensible_heat_flux',   # sshf (J/m², accum)
]

YEAR_START = 2003
YEAR_END = 2022
MONTH = 8


def main():
    c = cdsapi.Client()

    for year in range(YEAR_START, YEAR_END + 1):
        sl_dir = os.path.join(OUTDIR_BASE, f'era5_sl_{year}{MONTH:02d}')
        os.makedirs(sl_dir, exist_ok=True)
        outfile = os.path.join(sl_dir, 'data_stream-oper_stepType-accum_flux.nc')

        if os.path.exists(outfile) and os.path.getsize(outfile) > 1e6:
            print(f'{year}: skip (exists, {os.path.getsize(outfile)/1e6:.1f} MB)')
            continue

        print(f'{year}: downloading slhf+sshf...', flush=True)
        ndays = 31  # August
        days = [f'{d:02d}' for d in range(1, ndays + 1)]
        times = ['00:00', '06:00', '12:00', '18:00']

        try:
            c.retrieve(
                'reanalysis-era5-single-levels',
                {
                    'product_type': ['reanalysis'],
                    'variable': VARIABLES,
                    'year': [str(year)],
                    'month': [f'{MONTH:02d}'],
                    'day': days,
                    'time': times,
                    'data_format': 'netcdf',
                    'download_format': 'unarchived',
                    'area': AREA,
                },
                outfile,
            )
            fsize = os.path.getsize(outfile) / 1e6
            print(f'{year}: done ({fsize:.1f} MB)')
        except Exception as e:
            print(f'{year}: FAILED - {e}')
            if os.path.exists(outfile):
                os.remove(outfile)

    print('=== All done ===')


if __name__ == '__main__':
    main()
