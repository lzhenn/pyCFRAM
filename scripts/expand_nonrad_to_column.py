#!/usr/bin/env python3
"""Convert surface-mode nonrad_forcing.nc to column-mode.

Why this is needed
------------------
CFRAM's Planck-matrix inversion expects each non-radiative forcing term to be
provided as a *vertical profile* of layer-integrated heating (W/m^2 per layer).
GCM history files give this directly via tendency diagnostics (DTCOND, ZMDT,
DTV). Reanalyses and idealized RCM input usually give only a *surface* bulk
flux (slhf, sshf) — the column-internal heating is missing.

Putting a surface-only flux into the solver is physically valid only when the
column-internal heating is negligible (extreme heat events with no convection,
single-column RCE without clouds). For climate-mean / convection-active
scenarios, ignoring column heating dumps the whole tropical condensation
release into the atmospheric-dynamics residual — a known failure mode.

This tool detects that case and rebuilds nonrad_forcing.nc with a column-mode
distribution using the Lu & Cai (2009) Eq.(1) parameterization for latent
heat and an exponential PBL profile for sensible heat.

Tendency-direct path
--------------------
If the user has CAM (or other GCM) tendency diagnostics available, they should
*not* run this script — instead, generate nonrad_forcing.nc directly from
those tendencies (DTCOND/ZMDT/EVAPTZM → lhflx column; DTV → shflx column),
and tag the file with global attribute ``nonrad_source='tendency_direct'``.
This script refuses to overwrite tendency-direct files.

Usage
-----
    python3 scripts/expand_nonrad_to_column.py --case cesm2_4xco2
    python3 scripts/expand_nonrad_to_column.py --case eh22 --dry-run
    python3 scripts/expand_nonrad_to_column.py --in nonrad.nc --out nonrad_col.nc
"""
import os
import sys
import argparse
import shutil
import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case
from core.heating_profile import (
    distribute_to_column,
    expand_surface_to_column_energy_conserved,
)

FILL = -999.0
TOL = 900.0   # values larger than this in abs are treated as fill


def detect_mode(arr3d):
    """Return 'surface', 'column', or 'empty' for a (lev, lat, lon) variable.

    surface: only k=0 (sfc layer in the file convention) has valid data
    column:  more than one level has valid data
    empty:   no valid data anywhere
    """
    valid = np.abs(arr3d) < TOL
    valid_per_level = valid.any(axis=(1, 2))
    n_active = int(valid_per_level.sum())
    if n_active == 0:
        return 'empty'
    if n_active == 1:
        return 'surface'
    return 'column'


def read_var(nc, name):
    """Read a (time=1, lev, lat, lon) variable and squeeze time."""
    if name not in nc.variables:
        return None
    arr = np.array(nc.variables[name][0, :, :, :], dtype=np.float64)
    return arr


def write_var(nc, name, data4d, like=None):
    """Create or overwrite a (time, lev, lat, lon) variable."""
    if name in nc.variables:
        nc.variables[name][:] = data4d
        return
    dims = ('time', 'lev', 'lat', 'lon')
    v = nc.createVariable(name, 'f8', dims, fill_value=FILL)
    if like is not None and like in nc.variables:
        for attr in ('units', 'long_name'):
            try:
                setattr(v, attr, getattr(nc.variables[like], attr))
            except Exception:
                pass
    v[:] = data4d


def expand_file(in_path, out_path, lh_profile='lu_cai', sh_profile='pbl_exp',
                ps_hpa_2d=None, dry_run=False):
    """Expand surface-mode nonrad_forcing.nc to column-mode.

    Parameters
    ----------
    in_path, out_path : str
        Input/output NetCDF paths. Can be the same path (in-place overwrite).
    lh_profile, sh_profile : str
        Profile names for latent / sensible heat. See core.heating_profile.
    ps_hpa_2d : (nlat, nlon) or None
        Surface pressure for terrain-aware sh_profile. If None and pbl_exp
        is used, falls back to a uniform 1000 hPa.
    dry_run : bool
        Print plan and return without writing.
    """
    if not os.path.exists(in_path):
        sys.exit("Input not found: " + in_path)

    nc_in = Dataset(in_path, 'r')

    # Sanity: check existing source attribute
    src = getattr(nc_in, 'nonrad_source', None)
    if src == 'tendency_direct':
        nc_in.close()
        sys.exit("Refusing to overwrite tendency_direct file. "
                 "If you really want to, edit nonrad_source attribute first.")

    plev = np.array(nc_in.variables['lev'][:], dtype=np.float64)
    lats = np.array(nc_in.variables['lat'][:])
    lons = np.array(nc_in.variables['lon'][:])
    nlev, nlat, nlon = len(plev), len(lats), len(lons)

    print("=== Expanding nonrad forcing ===")
    print("File: %s -> %s" % (in_path, out_path))
    print("Grid: %d lev x %d lat x %d lon" % (nlev, nlat, nlon))
    print("Pressure: %.1f -> %.1f hPa" % (plev.min(), plev.max()))

    summary = {}
    new_data = {}

    for term in ('lhflx', 'shflx'):
        arr = read_var(nc_in, term)
        if arr is None:
            print("  %s: not present, skipping" % term)
            continue
        mode = detect_mode(arr)
        summary[term] = {'mode_in': mode}
        print("  %s: input mode=%s" % (term, mode))

        if mode == 'column':
            print("    -> already column, copying as-is")
            new_data[term] = arr
            summary[term]['mode_out'] = 'column'
            summary[term]['action'] = 'passthrough'
            continue

        if mode == 'empty':
            print("    -> empty, leaving zeros")
            new_data[term] = np.zeros_like(arr)
            summary[term]['mode_out'] = 'empty'
            summary[term]['action'] = 'zeros'
            continue

        # mode == 'surface' — find the level with valid data
        valid_per_level = (np.abs(arr) < TOL).any(axis=(1, 2))
        sfc_k = int(np.where(valid_per_level)[0][0])
        sfc_2d = arr[sfc_k, :, :]
        sfc_2d = np.where(np.abs(sfc_2d) < TOL, sfc_2d, 0.0)
        print("    -> sfc data at lev[%d]=%.1f hPa, mean=%.3f W/m^2" %
              (sfc_k, plev[sfc_k], float(np.nanmean(sfc_2d))))

        prof_name = lh_profile if term == 'lhflx' else sh_profile
        kwargs = {}
        ps_used = ps_hpa_2d
        if prof_name == 'pbl_exp' and ps_used is None:
            ps_used = np.full((nlat, nlon), 1000.0)
            print("    -> %s requires ps; using uniform 1000 hPa" % prof_name)

        # Use energy-conserving expansion (Lu & Cai 2009 Q^lh structure):
        # surface row keeps original sfc_flux (energy lost from surface),
        # atm rows hold -sfc_flux × profile (energy regained by column),
        # column total = 0 at every grid point.
        col = expand_surface_to_column_energy_conserved(
            sfc_2d, plev, ps_hpa=ps_used, profile=prof_name,
            surface_level_index=sfc_k)

        col_sum = col.sum(axis=0)
        max_err = float(np.max(np.abs(col_sum)))
        print("    -> energy-conserving expansion (column total ~0)")
        print("    -> column-sum |max| = %.3e W/m^2 (should be roundoff)" % max_err)

        new_data[term] = col
        summary[term]['mode_out'] = 'column'
        summary[term]['action'] = 'expanded_with_' + prof_name
        summary[term]['profile'] = prof_name
        summary[term]['sfc_input_level'] = float(plev[sfc_k])
        summary[term]['sfc_input_mean'] = float(np.nanmean(sfc_2d))
        summary[term]['column_max_abs_err'] = max_err

    nc_in.close()

    if dry_run:
        print("\n[DRY RUN] No file written.")
        return summary

    # Write output: copy structure from input, replace lhflx/shflx
    if os.path.abspath(in_path) != os.path.abspath(out_path):
        shutil.copy(in_path, out_path)

    nc_out = Dataset(out_path, 'r+')
    for term, arr in new_data.items():
        nc_out.variables[term][0, :, :, :] = arr

    # Tag provenance
    nc_out.nonrad_source = 'surface_expanded'
    nc_out.lhflx_profile = lh_profile
    nc_out.shflx_profile = sh_profile
    nc_out.expansion_tool = 'scripts/expand_nonrad_to_column.py'
    nc_out.expansion_summary = repr(summary)
    nc_out.close()

    print("\nWrote column-mode forcing -> %s" % out_path)
    print("Provenance attributes set: nonrad_source='surface_expanded'")
    return summary


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--case', help='Case name (uses cases/<name>/input/nonrad_forcing.nc)')
    g.add_argument('--in', dest='infile', help='Explicit input path')
    ap.add_argument('--out', dest='outfile',
                    help='Output path (default: in-place overwrite)')
    ap.add_argument('--lh-profile', default='lu_cai',
                    choices=['lu_cai', 'pbl_exp'],
                    help='Profile for latent heat (default: lu_cai)')
    ap.add_argument('--sh-profile', default='pbl_exp',
                    choices=['lu_cai', 'pbl_exp'],
                    help='Profile for sensible heat (default: pbl_exp)')
    ap.add_argument('--ps-from-surf', action='store_true',
                    help='Read ps from cases/<case>/input/base_surf.nc for terrain awareness')
    ap.add_argument('--dry-run', action='store_true')
    args = ap.parse_args()

    if args.case:
        cfg = load_case(args.case)
        in_path = cfg['input']['nonrad_forcing']
        out_path = args.outfile or in_path
        ps2d = None
        if args.ps_from_surf:
            surf_path = cfg['input']['base_surf']
            ncs = Dataset(surf_path)
            if 'ps' in ncs.variables:
                ps2d = np.array(ncs.variables['ps'][0, :, :], dtype=np.float64) / 100.0
            ncs.close()
    else:
        in_path = args.infile
        out_path = args.outfile or args.infile
        ps2d = None

    expand_file(in_path, out_path,
                lh_profile=args.lh_profile, sh_profile=args.sh_profile,
                ps_hpa_2d=ps2d, dry_run=args.dry_run)


if __name__ == '__main__':
    main()
