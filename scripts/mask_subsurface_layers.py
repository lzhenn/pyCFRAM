#!/usr/bin/env python3
"""Mask subsurface layers in a case's pres files based on ps.

Why
---
Standard pressure-level input files (e.g., 37 levels 1000→1 hPa) have layers
at lev > ps for high-elevation columns (Tibet, Andes, polar plateaus). Some
upstream interpolators fill these "subsurface" layers with smooth extrapolated
values (T, q, clouds), which RRTMG then ingests as if they were real
atmosphere — pathologizing the Planck matrix drdt_inv.

This tool sets, per column, all layers with lev > ps to:
  ta_lay = ts   (zero-thickness layer at ground temperature)
  q, cliq, cice, camt, o3, bc, ocphi, ocpho, sulf, ss, dust = 0
  co2 = unchanged (uniform field, keeps a small but nonzero Planck response)

Also clips globally negative cliq/cice (numerical advection ghosts) to 0.

Files modified in place: base_pres.nc, perturbed_pres.nc.
Backup written with .bak_premask suffix on first run only.
"""
import os
import sys
import argparse
import shutil
import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case

# Subsurface strategy:
#   HOLD = copy lowest-real-layer (just above ps) value down through subsurface.
#         RRTMG cannot handle exactly-zero H2O / O3 in any layer; holding the
#         last real value avoids the zero-discontinuity while keeping the
#         column physically below-ground (radiatively shadowed by surface).
#   ZERO = set to 0 (only for fields RRTMG tolerates as zero: clouds/aerosols).
MASK_VARS_HOLD = ['q', 'o3']
MASK_VARS_ZERO = ['cliq', 'cice', 'camt',
                  'bc', 'ocphi', 'ocpho', 'sulf', 'ss', 'dust']
CLIP_NEG_VARS = ['cliq', 'cice']


def mask_one(pres_path, surf_path, dry_run=False):
    """Apply ps mask to one (pres_file, surf_file) pair, in place."""
    nc_p = Dataset(pres_path, 'r' if dry_run else 'r+')
    nc_s = Dataset(surf_path, 'r')

    lev = np.array(nc_p.variables['lev'][:], dtype=np.float64)  # hPa
    ps = np.array(nc_s.variables['ps'][0, :, :], dtype=np.float64) / 100.0  # Pa→hPa
    ts = np.array(nc_s.variables['ts'][0, :, :], dtype=np.float64)
    nlev, nlat, nlon = len(lev), ps.shape[0], ps.shape[1]

    # Per-column mask: True where layer is below surface (lev > ps).
    submask_ps = lev[:, None, None] > ps[None, :, :]  # (nlev, nlat, nlon)
    # Also catch any cells with NaN/fillvalue in ta_lay regardless of ps
    # (CMIP6 plev19 fills extra cells based on a topography mask, beyond
    # what monthly ps comparison would catch).
    ta_raw = np.array(nc_p.variables['ta_lay'][0, :, :, :], dtype=np.float64)
    submask_nan = np.isnan(ta_raw) | (np.abs(ta_raw) > 1e15)
    submask_3d = submask_ps | submask_nan
    n_sub = int(submask_3d.sum())
    n_extra_nan = int((submask_nan & ~submask_ps).sum())
    print('  %s: %d/%d subsurface cells (%.1f%%; +%d NaN-only beyond ps mask)' %
          (os.path.basename(pres_path), n_sub, submask_3d.size,
           100.0 * n_sub / submask_3d.size, n_extra_nan))

    summary = {}

    # ta_lay → ts in subsurface (reuse already-loaded ta_raw)
    ta = ta_raw
    ta_old_min = float(np.nanmin(ta[submask_3d])) if n_sub > 0 else float('nan')
    ta_old_max = float(np.nanmax(ta[submask_3d])) if n_sub > 0 else float('nan')
    ts_b = np.broadcast_to(ts[None, :, :], ta.shape)
    ta_new = np.where(submask_3d, ts_b, ta)
    summary['ta_lay'] = (ta_old_min, ta_old_max, ta_new[submask_3d].min(), ta_new[submask_3d].max())

    masked_arrays = {'ta_lay': ta_new}

    # Helper: index of "lowest real layer" = layer just above ground (highest
    # pressure among non-subsurface layers). Convention here: lev is in hPa
    # and may run either ascending or descending. The lowest-real layer is
    # the one whose lev is the LARGEST satisfying lev <= ps.
    real_mask_3d = ~submask_3d  # (nlev, nlat, nlon)
    # For each (lat, lon) column, find argmax of lev among real layers.
    lev_3d = lev[:, None, None] * np.ones_like(submask_3d, dtype=np.float64)
    lev_masked = np.where(real_mask_3d, lev_3d, -np.inf)
    lowest_real_k = lev_masked.argmax(axis=0)  # (nlat, nlon)

    # HOLD strategy: copy value at lowest_real_k down into subsurface layers
    for v in MASK_VARS_HOLD:
        if v in nc_p.variables:
            arr = np.array(nc_p.variables[v][0, :, :, :], dtype=np.float64)
            # Build held value (nlat, nlon)
            held = np.take_along_axis(arr, lowest_real_k[None, :, :], axis=0)[0]  # (nlat, nlon)
            held_b = np.broadcast_to(held[None, :, :], arr.shape)
            arr_new = np.where(submask_3d, held_b, arr)
            masked_arrays[v] = arr_new
            summary[v + '_hold_sample'] = float(held.mean())

    # ZERO strategy: zero out subsurface (clouds/aerosols)
    for v in MASK_VARS_ZERO:
        if v in nc_p.variables:
            arr = np.array(nc_p.variables[v][0, :, :, :], dtype=np.float64)
            old_max = arr[submask_3d].max() if n_sub > 0 else 0.0
            arr_new = np.where(submask_3d, 0.0, arr)
            # global clip negatives for cliq/cice
            if v in CLIP_NEG_VARS:
                neg_count = int((arr_new < 0).sum())
                if neg_count > 0:
                    arr_new = np.where(arr_new < 0, 0.0, arr_new)
                    summary[v + '_neg_clipped'] = neg_count
            masked_arrays[v] = arr_new
            summary[v + '_sub_max'] = old_max

    if dry_run:
        print('  [DRY RUN] no write')
        nc_p.close()
        nc_s.close()
        return summary

    # Write back
    for v, arr in masked_arrays.items():
        nc_p.variables[v][0, :, :, :] = arr

    # Tag provenance
    nc_p.subsurface_masked = 'true'
    nc_p.subsurface_mask_tool = 'scripts/mask_subsurface_layers.py'

    nc_p.close()
    nc_s.close()
    return summary


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', required=True)
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--no-backup', action='store_true',
                    help='skip backup creation (use if already backed up)')
    args = ap.parse_args()

    cfg = load_case(args.case)
    pairs = [
        (cfg['input']['base_pres'],      cfg['input']['base_surf']),
        (cfg['input']['perturbed_pres'], cfg['input']['perturbed_surf']),
    ]

    print('=== Subsurface masking: %s ===' % args.case)
    for pres, surf in pairs:
        if not os.path.exists(pres):
            sys.exit('Missing: ' + pres)
        if not os.path.exists(surf):
            sys.exit('Missing: ' + surf)

        # Check for prior masking
        nc = Dataset(pres)
        already = getattr(nc, 'subsurface_masked', None)
        nc.close()
        if already == 'true':
            print('  %s already tagged subsurface_masked, skipping' %
                  os.path.basename(pres))
            continue

        # Backup
        bak = pres + '.bak_premask'
        if not args.dry_run and not args.no_backup and not os.path.exists(bak):
            shutil.copy(pres, bak)
            print('  backup: ' + bak)

        summary = mask_one(pres, surf, dry_run=args.dry_run)
        print('    ta_lay subsurface T: was [%.1f, %.1f] -> now [%.1f, %.1f] K' %
              summary['ta_lay'])
        for k, v in summary.items():
            if k.endswith('_neg_clipped'):
                print('    %s: %d negative cells clipped to 0' % (k, v))
            elif k.endswith('_sub_max'):
                print('    %s subsurface max (pre-mask): %.3e' %
                      (k.replace('_sub_max', ''), v))

    print('Done.')


if __name__ == '__main__':
    main()
