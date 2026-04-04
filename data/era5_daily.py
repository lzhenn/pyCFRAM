"""ERA5 daily data loader and preprocessor for CFRAM-A.

Implements the Wu et al. (2025) methodology:
1. Read ERA5 daily data (6-hourly → daily mean)
2. Compute daily climatology (2003-2022 mean for each calendar day)
3. Identify warm periods (Tmax > 90th percentile for 3+ consecutive days)
4. Compute base state (climatological warm period mean)
5. Compute anomalous state (event year warm period mean)
"""

import os
import numpy as np
from netCDF4 import Dataset, num2date


class ERA5DailyLoader:
    """Load and preprocess ERA5 daily data for CFRAM analysis."""

    def __init__(self, data_dir, years=range(2003, 2023)):
        """
        Args:
            data_dir: directory containing era5_pl_{var}_{YYYY}08.nc and
                      era5_sl_{YYYY}08.nc files
            years: range of years for climatology
        """
        self.data_dir = data_dir
        self.years = list(years)

    def load_pl_variable(self, short_name, year):
        """Load one pressure-level variable for one August.

        Args:
            short_name: 't', 'q', 'o3', 'cc', 'clwc', 'ciwc'
            year: integer year

        Returns:
            data: (ntime, nlev, nlat, nlon) daily mean array
            lats, lons, levels: coordinate arrays
        """
        fname = os.path.join(self.data_dir, f'era5_pl_{short_name}_{year}08.nc')
        nc = Dataset(fname, 'r')

        # Auto-detect variable name
        skip = {'time', 'latitude', 'longitude', 'level', 'pressure_level',
                'lat', 'lon', 'valid_time', 'number', 'expver'}
        var_name = None
        for v in nc.variables:
            if v.lower() not in skip and nc.variables[v].ndim >= 3:
                var_name = v
                break

        data_raw = np.array(nc.variables[var_name][:], dtype=np.float64)

        # Coordinate arrays
        for name in ('latitude', 'lat'):
            if name in nc.variables:
                lats = np.array(nc.variables[name][:])
                break
        for name in ('longitude', 'lon'):
            if name in nc.variables:
                lons = np.array(nc.variables[name][:])
                break
        for name in ('level', 'pressure_level'):
            if name in nc.variables:
                levels = np.array(nc.variables[name][:])
                break

        nc.close()

        # If 6-hourly (4 times per day), compute daily mean
        ndays = 31  # August
        ntimes_per_day = len(data_raw) // ndays
        if ntimes_per_day > 1:
            shape = (ndays, ntimes_per_day) + data_raw.shape[1:]
            data = data_raw.reshape(shape).mean(axis=1)
        else:
            data = data_raw

        return data, lats, lons, levels

    def load_sl_variable(self, var_name, year):
        """Load one single-level variable for one August.

        Args:
            var_name: NetCDF variable name (auto-detected if ambiguous)
            year: integer year

        Returns:
            data: (ndays, nlat, nlon) daily mean array
            lats, lons: coordinate arrays
        """
        fname = os.path.join(self.data_dir, f'era5_sl_{year}08.nc')
        nc = Dataset(fname, 'r')

        # Try exact name, then search
        if var_name in nc.variables:
            data_raw = np.array(nc.variables[var_name][:], dtype=np.float64)
        else:
            # Search by standard name or long name
            data_raw = None
            for v in nc.variables:
                if var_name.lower() in v.lower():
                    data_raw = np.array(nc.variables[v][:], dtype=np.float64)
                    break
            if data_raw is None:
                nc.close()
                raise KeyError(f"Variable '{var_name}' not found in {fname}")

        for name in ('latitude', 'lat'):
            if name in nc.variables:
                lats = np.array(nc.variables[name][:])
                break
        for name in ('longitude', 'lon'):
            if name in nc.variables:
                lons = np.array(nc.variables[name][:])
                break

        nc.close()

        # 6-hourly → daily mean
        ndays = 31
        ntimes_per_day = len(data_raw) // ndays
        if ntimes_per_day > 1:
            shape = (ndays, ntimes_per_day) + data_raw.shape[1:]
            data = data_raw.reshape(shape).mean(axis=1)
        else:
            data = data_raw

        return data, lats, lons

    def compute_daily_climatology_pl(self, short_name):
        """Compute daily climatology for a PL variable over all years.

        Returns:
            clim: (31, nlev, nlat, nlon) climatological daily mean for August
            lats, lons, levels
        """
        all_data = []
        for yr in self.years:
            try:
                data, lats, lons, levels = self.load_pl_variable(short_name, yr)
                all_data.append(data[:31])  # ensure 31 days
            except Exception as e:
                print(f"  Warning: skip {short_name} {yr}: {e}")
        # (nyears, 31, nlev, nlat, nlon) → mean over years
        stacked = np.array(all_data)
        clim = stacked.mean(axis=0)
        return clim, lats, lons, levels

    def compute_daily_climatology_sl(self, var_name):
        """Compute daily climatology for a SL variable."""
        all_data = []
        for yr in self.years:
            try:
                data, lats, lons = self.load_sl_variable(var_name, yr)
                all_data.append(data[:31])
            except Exception as e:
                print(f"  Warning: skip {var_name} {yr}: {e}")
        stacked = np.array(all_data)
        clim = stacked.mean(axis=0)
        return clim, lats, lons

    def identify_warm_period(self, year, lat_range, lon_range, threshold_pct=90,
                             min_consecutive=3):
        """Identify warm period dates for an extreme heat event.

        Uses maximum 2m temperature (Tmax) exceeding the 90th percentile
        for 3+ consecutive days.

        Args:
            year: event year
            lat_range: (lat_min, lat_max) for key region
            lon_range: (lon_min, lon_max) for key region
            threshold_pct: percentile threshold (default 90)
            min_consecutive: minimum consecutive days (default 3)

        Returns:
            warm_days: list of day indices (0-based, August)
        """
        # Load Tmax for all years to compute threshold
        all_tmax = []
        for yr in self.years:
            try:
                data, lats, lons = self.load_sl_variable('mx2t', yr)
                # Select key region
                lat_mask = (lats >= lat_range[0]) & (lats <= lat_range[1])
                lon_mask = (lons >= lon_range[0]) & (lons <= lon_range[1])
                regional = data[:31, :, :][:, lat_mask, :][:, :, lon_mask]
                # Area mean
                all_tmax.append(regional.mean(axis=(1, 2)))  # (31,)
            except Exception:
                pass

        # Stack: (nyears, 31) → percentile for each calendar day
        tmax_stack = np.array(all_tmax)  # (nyears, 31)

        # Compute threshold: for each day, use ±7 day window across all years
        # Simplified: use all years' same-day values
        threshold = np.percentile(tmax_stack, threshold_pct, axis=0)  # (31,)

        # Event year Tmax
        event_data, lats, lons = self.load_sl_variable('mx2t', year)
        lat_mask = (lats >= lat_range[0]) & (lats <= lat_range[1])
        lon_mask = (lons >= lon_range[0]) & (lons <= lon_range[1])
        event_tmax = event_data[:31, :, :][:, lat_mask, :][:, :, lon_mask].mean(axis=(1, 2))

        # Find consecutive days exceeding threshold
        exceed = event_tmax > threshold
        warm_days = []
        consecutive = 0
        start = -1
        for d in range(31):
            if exceed[d]:
                if consecutive == 0:
                    start = d
                consecutive += 1
            else:
                if consecutive >= min_consecutive:
                    warm_days.extend(range(start, start + consecutive))
                consecutive = 0
                start = -1
        if consecutive >= min_consecutive:
            warm_days.extend(range(start, start + consecutive))

        return sorted(set(warm_days))

    def get_base_warm_states(self, event_year, warm_days,
                              lat_idx, lon_idx):
        """Extract base and warm state profiles for a single column.

        Args:
            event_year: year of the extreme event
            warm_days: list of day indices (0-based August)
            lat_idx: latitude index
            lon_idx: longitude index

        Returns:
            base_state: dict of climatological warm-period mean profiles
            warm_state: dict of event-year warm-period mean profiles
        """
        pl_vars = {
            't': 't', 'q': 'q', 'o3': 'o3',
            'cc': 'cldfrac', 'clwc': 'cldlwc', 'ciwc': 'cldiwc',
        }

        base_state = {}
        warm_state = {}

        # Pressure level variables
        for short, state_key in pl_vars.items():
            clim, lats, lons, levels = self.compute_daily_climatology_pl(short)
            # Climatological warm period mean at (lat_idx, lon_idx)
            base_profile = clim[warm_days, :, lat_idx, lon_idx].mean(axis=0)

            # Event year warm period mean
            event_data, _, _, _ = self.load_pl_variable(short, event_year)
            warm_profile = event_data[warm_days, :, lat_idx, lon_idx].mean(axis=0)

            base_state[state_key] = base_profile
            warm_state[state_key] = warm_profile

        base_state['plev'] = levels
        warm_state['plev'] = levels

        # Single level variables
        sl_mapping = {
            'skt': 'ts',
            'sp': 'ps',
            'tisr': 'solin',
            'ssrd': 'ssrd',
        }
        for sl_var, state_key in sl_mapping.items():
            clim_sl, _, _ = self.compute_daily_climatology_sl(sl_var)
            base_val = clim_sl[warm_days, lat_idx, lon_idx].mean()

            event_sl, _, _ = self.load_sl_variable(sl_var, event_year)
            warm_val = event_sl[warm_days, lat_idx, lon_idx].mean()

            base_state[state_key] = base_val
            warm_state[state_key] = warm_val

        # Derived quantities
        # ssru = ssrd - ssr
        for period, state in [('clim', base_state), ('event', warm_state)]:
            if period == 'clim':
                ssr_data, _, _ = self.compute_daily_climatology_sl('ssr')
                ssrd_data, _, _ = self.compute_daily_climatology_sl('ssrd')
                ssr_val = ssr_data[warm_days, lat_idx, lon_idx].mean()
                ssrd_val = ssrd_data[warm_days, lat_idx, lon_idx].mean()
            else:
                ssr_event, _, _ = self.load_sl_variable('ssr', event_year)
                ssrd_event, _, _ = self.load_sl_variable('ssrd', event_year)
                ssr_val = ssr_event[warm_days, lat_idx, lon_idx].mean()
                ssrd_val = ssrd_event[warm_days, lat_idx, lon_idx].mean()

            # Unit conversion: J/m² → W/m² for accumulated fields
            ssrd_wm2 = ssrd_val / 86400.0 if ssrd_val > 1000 else ssrd_val
            ssr_wm2 = ssr_val / 86400.0 if ssr_val > 1000 else ssr_val
            ssru = ssrd_wm2 - ssr_wm2

            state['ssrd'] = ssrd_wm2
            state['ssru'] = ssru
            solin = state.get('solin', 0)
            state['solin'] = solin / 86400.0 if solin > 1000 else solin

        # Compute derived fields
        for state in (base_state, warm_state):
            solin = state['solin']
            state['zenith'] = solin / 1360.98 if solin > 0 else 0.0
            ssrd = state['ssrd']
            ssru = state['ssru']
            state['albedo_sw'] = ssru / ssrd if ssrd > 0 else 0.2
            state['albedo_lw'] = 0.0

            # Placeholder aerosol (zero until MERRA-2 processing)
            nlev = len(state.get('t', state.get('plev', [])))
            state['tauaer_sw'] = np.zeros((nlev, 14))
            state['ssaaer_sw'] = np.zeros((nlev, 14))
            state['asmaer_sw'] = np.zeros((nlev, 14))
            state['tauaer_lw'] = np.zeros((nlev, 16))

        return base_state, warm_state
