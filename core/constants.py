"""Physical constants and RRTMG-intrinsic parameters.

These are universal constants that never change between cases.
Case-specific and configurable parameters are in configs/defaults.yaml.
"""

import numpy as np

# --- Physical constants ---
GRAVITY = 9.8066          # m/s^2
R_GAS = 8.31446           # J/(mol*K), universal gas constant
R_DRY = 287.058           # J/(kg*K), dry air gas constant
CP_DRY = 1004.0           # J/(kg*K), specific heat at constant pressure
MW_DRY = 28.966e-3        # kg/mol, molecular weight of dry air
MW_H2O = 18.016e-3        # kg/mol, molecular weight of water vapor
MW_O3 = 47.998e-3         # kg/mol, molecular weight of ozone
MW_CO2 = 44.010e-3        # kg/mol, molecular weight of CO2
MW_CH4 = 16.043e-3        # kg/mol, molecular weight of CH4
MW_N2O = 44.013e-3        # kg/mol, molecular weight of N2O
AVOGADRO = 6.02214e23     # molecules/mol

# --- RRTMG intrinsic spectral band counts (fixed by RRTMG code) ---
NBND_LW = 16     # longwave bands
NBND_SW = 14     # shortwave bands
NBND_TOTAL = 30  # SW(1-14) + LW(15-30)
