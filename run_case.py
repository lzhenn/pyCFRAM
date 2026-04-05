#!/usr/bin/env python3
"""pyCFRAM entry point: run CFRAM decomposition for a case.

Usage:
    python3 run_case.py eh13                    # full workflow
    python3 run_case.py eh13 --step extract     # only extract input
    python3 run_case.py eh13 --step run         # only run CFRAM
    python3 run_case.py eh13 --step plot        # only plot results
    python3 run_case.py eh13 --nproc 20         # override parallel workers
"""
import os
import sys
import argparse
import subprocess

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from core.config import load_case, get_nproc, PROJECT_ROOT


def run_step(script, args_list):
    """Run a script with arguments."""
    cmd = [sys.executable, os.path.join(PROJECT_ROOT, 'scripts', script)] + args_list
    print("\n>>> %s" % ' '.join(cmd))
    ret = subprocess.call(cmd)
    if ret != 0:
        print("ERROR: %s failed with exit code %d" % (script, ret))
        sys.exit(ret)


def main():
    parser = argparse.ArgumentParser(description='pyCFRAM: run CFRAM decomposition')
    parser.add_argument('case', help='Case name (directory under cases/)')
    parser.add_argument('--step', choices=['extract', 'run', 'plot', 'all'],
                        default='all', help='Which step to run (default: all)')
    parser.add_argument('--nproc', type=int, default=None,
                        help='Number of parallel workers (default: from config or all CPUs)')
    args = parser.parse_args()

    # Validate case exists
    cfg = load_case(args.case)
    nproc = args.nproc or get_nproc(cfg)
    print("=" * 60)
    print("pyCFRAM: %s (%s)" % (cfg.get('case_name', args.case), cfg.get('description', '')))
    print("Workers: %d" % nproc)
    print("=" * 60)

    # Check input files exist
    for key in ['base_pres', 'base_surf', 'perturbed_pres', 'perturbed_surf']:
        path = cfg['input'].get(key)
        if not path or not os.path.exists(path):
            print("ERROR: Input file missing: %s" % path)
            print("Run: python3 scripts/prepare_from_paper_data.py --help")
            sys.exit(1)

    steps = ['extract', 'run', 'plot'] if args.step == 'all' else [args.step]

    if 'extract' in steps:
        run_step('extract_full_field.py', ['--case', args.case])

    if 'run' in steps:
        run_step('run_parallel_python.py', ['--case', args.case, '--nproc', str(nproc)])

    if 'plot' in steps:
        run_step('plot_fig3_self.py', ['--case', args.case])

    print("\n" + "=" * 60)
    print("Done! Results in: cases/%s/output/ and cases/%s/figures/" % (args.case, args.case))
    print("=" * 60)


if __name__ == '__main__':
    main()
