"""Physical constants and configuration for CFRAM-A calculations."""

import numpy as np

# --- Physical constants ---
GRAVITY = 9.8066          # m/s^2
R_GAS = 8.31446           # J/(mol·K), universal gas constant
R_DRY = 287.058           # J/(kg·K), dry air gas constant
CP_DRY = 1004.0           # J/(kg·K), specific heat at constant pressure
MW_DRY = 28.966e-3        # kg/mol, molecular weight of dry air
MW_H2O = 18.016e-3        # kg/mol, molecular weight of water vapor
MW_O3 = 47.998e-3         # kg/mol, molecular weight of ozone
MW_CO2 = 44.010e-3        # kg/mol, molecular weight of CO2
MW_CH4 = 16.043e-3        # kg/mol, molecular weight of CH4
MW_N2O = 44.013e-3        # kg/mol, molecular weight of N2O
AVOGADRO = 6.02214e23     # molecules/mol
SOLAR_CONSTANT = 1360.98  # W/m^2 (from cfram_rrtmg.f90)

# --- RRTMG spectral bands ---
NBND_LW = 16   # longwave bands
NBND_SW = 14   # shortwave bands
NBND_TOTAL = 30  # SW(1-14) + LW(15-30)

# --- ERA5 standard pressure levels (37 levels, hPa, top to bottom) ---
ERA5_PLEVELS = np.array([
    1, 2, 3, 5, 7, 10, 20, 30, 50, 70,
    100, 125, 150, 175, 200, 225, 250, 300, 350, 400,
    450, 500, 550, 600, 650, 700, 750, 775, 800, 825,
    850, 875, 900, 925, 950, 975, 1000
], dtype=np.float64)

NLEV = 37  # number of vertical levels

# --- Default gas concentrations (from cfram_rrtmg.f90) ---
CH4_PPMV = 1.6    # methane
N2O_PPMV = 0.28   # nitrous oxide

# --- RRTMG configuration flags ---
ICLD_OFF = 0     # no clouds
ICLD_ON = 2      # clouds enabled (McICA)
IAER_OFF = 0     # no aerosols
IAER_ON = 10     # aerosols from input (band-by-band)

# --- Cloud effective sizes (from rad_driver.f90) ---
CLOUD_RE_ICE = 5.0    # microns, ice effective radius
CLOUD_RE_LIQ = 20.0   # microns, liquid effective radius

# --- Study domain (Wu et al. 2025) ---
DOMAIN_LON = (95.0, 125.0)   # degrees East
DOMAIN_LAT = (18.0, 42.0)    # degrees North

# Key areas for EH13 and EH22
KEY_AREA_EH13 = {'lon': (110.0, 122.0), 'lat': (28.0, 36.0)}
KEY_AREA_EH22 = {'lon': (103.0, 122.0), 'lat': (27.0, 35.0)}

# --- MERRA-2 aerosol species ---
MERRA2_AEROSOL_VARS = {
    'DU': ['DU001', 'DU002', 'DU003', 'DU004', 'DU005'],
    'SS': ['SS001', 'SS002', 'SS003', 'SS004', 'SS005'],
    'BC': ['BCPHOBIC', 'BCPHILIC'],
    'OC': ['OCPHOBIC', 'OCPHILIC'],
    'SU': ['SO4'],
}

MERRA2_AUX_VARS = ['AIRDENS', 'RH', 'DELP']

# --- CFRAM decomposition terms ---
CFRAM_TERMS = [
    'warm',     # total (all variables changed)
    'co2',      # CO2 forcing
    'q',        # water vapor
    'ts',       # surface temperature
    'o3',       # ozone
    'solar',    # solar irradiance / zenith angle
    'albedo',   # surface albedo
    'cloud',    # cloud properties
    'aerosol',  # aerosol properties
]
