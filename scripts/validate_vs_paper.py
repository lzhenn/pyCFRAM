#!/usr/bin/env python3
"""Correct surface-level comparison: surface is at index -1, not 0."""

import os, sys, numpy as np
from netCDF4 import Dataset
from scipy.stats import pearsonr

BASE = '/disk/r074/lzhenn/workspace/ust-jumper/pyCFRAM'

# Accept case from command line, default to eh22
case = sys.argv[1] if len(sys.argv) > 1 else 'eh22'
case_upper = case.upper()

# Map case to paper_data directory
paper_dirs = {
    'eh22': 'case_eh22_c20250118',
    'eh13': 'case_eh13_c20250102',
}
if case not in paper_dirs:
    print(f"Unknown case: {case}. Available: {list(paper_dirs.keys())}")
    sys.exit(1)

ind_path = os.path.join(BASE, f'cases/{case}/output/cfram_result.nc')
paper_path = os.path.join(BASE, f'paper_data/cfram_out/{paper_dirs[case]}/merra2_{case}_partial_t.nc')

def get_surface(nc, v):
    """Get surface-level field. In cfram_result.nc, surface = last lev index."""
    d = np.array(nc.variables[v][:])
    if d.ndim == 4:
        return d[0, -1, :, :]
    elif d.ndim == 3:
        return d[-1, :, :]
    return d

def get_paper_surface(nc, v):
    """Paper partial_t: surface = first lev index (lev[0]=1013 hPa)."""
    d = np.array(nc.variables[v][:])
    if d.ndim == 4:
        sfc = d[0, 0, :, :]
        valid = sfc[sfc > -900]
        return sfc, float(valid.mean()) if len(valid) > 0 else float('nan')
    elif d.ndim == 3:
        sfc = d[0, :, :]
        valid = sfc[sfc > -900]
        return sfc, float(valid.mean()) if len(valid) > 0 else float('nan')
    return d, float(d.mean())

print(f"=== {case_upper}: Checking level structure of independent result ===")
inc = Dataset(ind_path)
levs = np.array(inc.variables['lev'][:])
print(f"lev: {levs[:5]}...{levs[-5:]}")
print(f"lev[0]={levs[0]:.0f} hPa, lev[-1]={levs[-1]:.0f} hPa")
print(f"Shape of dT_cloud: {inc.variables['dT_cloud'][:].shape}")
print(f"dT_cloud at lev[0] ({levs[0]:.0f} hPa): mean={inc.variables['dT_cloud'][0,:,:].mean():.4f}")
print(f"dT_cloud at lev[-1] ({levs[-1]:.0f} hPa): mean={inc.variables['dT_cloud'][-1,:,:].mean():.4f}")
inc.close()

print()
print(f"=== {case_upper}: Paper result level structure ===")
pnc = Dataset(paper_path)
plevs = np.array(pnc.variables['lev'][:])
print(f"lev: {plevs[:5]}...{plevs[-5:]}")
print(f"lev[0]={plevs[0]:.0f} hPa, lev[-1]={plevs[-1]:.0f} hPa")
pnc.close()

print()
print(f"=== {case_upper}: CORRECTED comparison: surface dT ===")
print(f"  {'term':<15} {'paper':>10} {'indep':>10} {'corr':>8}")

pnc = Dataset(paper_path)
plat = np.array(pnc.variables['lat'][:])

inc = Dataset(ind_path)
ilat = np.array(inc.variables['lat'][:])
lat_mask = (ilat >= 20) & (ilat <= 40)

term_map = [
    ('cld', 'dT_cloud'), ('wv', 'dT_q'), ('co2', 'dT_co2'),
    ('albedo', 'dT_albedo'), ('solar', 'dT_solar'),
    ('lhflx', 'dT_lhflx'), ('shflx', 'dT_shflx'),
]

for pvar, ivar in term_map:
    if pvar not in pnc.variables or ivar not in inc.variables:
        continue

    pfield_4d, _ = get_paper_surface(pnc, pvar)
    pfield = pfield_4d
    valid_mask = pfield > -900
    pmean = float(pfield[valid_mask].mean()) if valid_mask.any() else float('nan')

    ifield_3d = np.array(inc.variables[ivar][:])
    ifield_sfc = ifield_3d[-1, :, :]
    ifield_crop = ifield_sfc[lat_mask, :]
    imean = float(ifield_crop.mean())

    pf = pfield.flatten()
    if_crop = ifield_crop.flatten()
    if len(pf) == len(if_crop):
        valid = (pf > -900)
        r, _ = pearsonr(pf[valid], if_crop[valid])
    else:
        r = float('nan')

    print(f"  {pvar:<15} {pmean:>10.4f} {imean:>10.4f} {r:>8.3f}")

# Aerosol total
print()
print(f"=== {case_upper}: Aerosol comparison ===")
pnc2 = Dataset(paper_path)
paper_aer_terms = ['bc', 'oc', 'sulf', 'seas', 'dust']
paper_aer_total = None
for pvar in paper_aer_terms:
    if pvar in pnc2.variables:
        pf, pm = get_paper_surface(pnc2, pvar)
        print(f"  paper {pvar}: {pm:.4f} K")
        if paper_aer_total is None:
            paper_aer_total = np.where(pf > -900, pf, 0.0)
        else:
            paper_aer_total += np.where(pf > -900, pf, 0.0)
if paper_aer_total is not None:
    print(f"  paper total aerosol: {paper_aer_total.mean():.4f} K")
pnc2.close()

ifield_aer = np.array(inc.variables['dT_aerosol'][:])[-1, :, :]
ifield_aer_crop = ifield_aer[lat_mask, :]
print(f"  indep dT_aerosol: {ifield_aer_crop.mean():.4f} K")

inc.close()
pnc.close()
