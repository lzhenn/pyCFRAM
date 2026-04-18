#!/usr/bin/env python3
"""Phase 5 verification: domain-mean per-species surface dT, self vs paper.

Regrids self output onto paper grid (nearest-neighbor) before comparing.
Key test: spatial correlation > 0.6 on overlapping region for each species.

Usage:
    python3 scripts/verify_phase5_vs_paper.py --case eh22 eh13
"""
import os, sys, argparse
import numpy as np
from netCDF4 import Dataset

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import load_case, PROJECT_ROOT

PAPER_CASES = {
    'eh13': 'case_eh13_c20250102',
    'eh22': 'case_eh22_c20250118',
}


def load_paper(case_name):
    paper_dir = os.path.join(PROJECT_ROOT, 'paper_data', 'cfram_out', PAPER_CASES[case_name])
    pt_file = [f for f in os.listdir(paper_dir) if 'partial_t' in f][0]
    nc = Dataset(os.path.join(paper_dir, pt_file))
    lats = np.array(nc.variables['lat'][:])
    lons = np.array(nc.variables['lon'][:])
    data = {}
    for sp in ['bc', 'oc', 'sulf', 'seas', 'dust']:
        if sp in nc.variables:
            a = np.array(nc.variables[sp][0, 0, :, :], dtype=np.float64)
            data[sp] = np.where(np.abs(a) > 900, np.nan, a)
    nc.close()
    return lats, lons, data


def load_self(case_name):
    cfg = load_case(case_name)
    result = os.path.join(cfg['_output_dir'], 'cfram_result.nc')
    nc = Dataset(result)
    lats = np.array(nc.variables['lat'][:])
    lons = np.array(nc.variables['lon'][:])

    def ld(v):
        if v not in nc.variables: return None
        a = np.array(nc.variables[v][-1, :, :], dtype=np.float64)
        return np.where(np.abs(a) > 900, np.nan, a)

    data = {}
    data['bc']   = ld('dT_bc')
    data['sulf'] = ld('dT_sulf')
    data['dust'] = ld('dT_dust')
    ocphi = ld('dT_ocphi'); ocpho = ld('dT_ocpho')
    if ocphi is not None and ocpho is not None:
        data['oc'] = np.nan_to_num(ocphi) + np.nan_to_num(ocpho)
    data['seas'] = ld('dT_ss')
    nc.close()
    return lats, lons, data


def regrid_nearest(src_lat, src_lon, src_val, tgt_lat, tgt_lon):
    out = np.full((len(tgt_lat), len(tgt_lon)), np.nan)
    for i, la in enumerate(tgt_lat):
        ilat = int(np.argmin(np.abs(src_lat - la)))
        for j, lo in enumerate(tgt_lon):
            ilon = int(np.argmin(np.abs(src_lon - lo)))
            out[i, j] = src_val[ilat, ilon]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--case', nargs='+', required=True)
    args = ap.parse_args()

    for cn in args.case:
        print('\n=== %s ===' % cn)
        p_lat, p_lon, pdata = load_paper(cn)
        s_lat, s_lon, sdata = load_self(cn)
        print('paper grid: %d x %d; self grid: %d x %d' % (
            len(p_lat), len(p_lon), len(s_lat), len(s_lon)))
        print('%-8s %12s %12s %10s %10s %10s' % ('species', 'paper_mean', 'self_mean', 'corr', 'rmse', 'bias'))
        for sp in ['bc', 'oc', 'sulf', 'seas', 'dust']:
            if sp not in pdata or sp not in sdata or sdata[sp] is None:
                print('  %s MISSING' % sp)
                continue
            p = pdata[sp]
            s_regrid = regrid_nearest(s_lat, s_lon, sdata[sp], p_lat, p_lon)
            mask = np.isfinite(p) & np.isfinite(s_regrid)
            if mask.sum() == 0:
                print('  %s no overlap' % sp)
                continue
            pm = float(np.nanmean(p))
            sm = float(np.nanmean(s_regrid))
            corr = float(np.corrcoef(p[mask], s_regrid[mask])[0, 1])
            rmse = float(np.sqrt(np.mean((p[mask] - s_regrid[mask])**2)))
            bias = float(np.mean(s_regrid[mask] - p[mask]))
            print('  %-8s %+12.4f %+12.4f %10.4f %10.4f %+10.4f' % (sp, pm, sm, corr, rmse, bias))


if __name__ == '__main__':
    main()
