#!/usr/bin/env python3
"""Convert paper_data format to pyCFRAM standard input format.

Creates symlinks (or copies) from paper_data input_check.nc files
to the standard base_pres/base_surf/perturbed_pres/perturbed_surf names.
Also extracts non-radiative forcing from partial_forcing.nc.

Usage:
    python3 scripts/prepare_from_paper_data.py --paper_dir paper_data/cfram_out/case_eh13_c20250102 --case eh13
"""
import os
import sys
import argparse
import glob

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from core.config import PROJECT_ROOT


def find_file(directory, pattern):
    """Find a file matching pattern in directory."""
    matches = glob.glob(os.path.join(directory, pattern))
    if not matches:
        return None
    return matches[0]


def main():
    parser = argparse.ArgumentParser(description='Convert paper_data to pyCFRAM input format')
    parser.add_argument('--paper_dir', required=True,
                        help='Path to paper_data case directory (e.g., paper_data/cfram_out/case_eh13_c20250102)')
    parser.add_argument('--case', required=True,
                        help='Case name (e.g., eh13)')
    parser.add_argument('--copy', action='store_true',
                        help='Copy files instead of creating symlinks')
    args = parser.parse_args()

    paper_dir = os.path.abspath(args.paper_dir)
    if not os.path.isabs(args.paper_dir):
        paper_dir = os.path.join(PROJECT_ROOT, args.paper_dir)

    case_input = os.path.join(PROJECT_ROOT, 'cases', args.case, 'input')
    os.makedirs(case_input, exist_ok=True)

    # Mapping: standard name -> paper_data file pattern
    mapping = {
        'base_pres.nc':      '*baseline_pres*input_check*',
        'base_surf.nc':      '*baseline_surf*input_check*',
        'perturbed_pres.nc': '*all_pres*input_check*',
        'perturbed_surf.nc': '*all_surf*input_check*',
        'nonrad_forcing.nc': '*partial_forcing*',
    }

    print("Paper data: %s" % paper_dir)
    print("Case input: %s" % case_input)
    print()

    for std_name, pattern in mapping.items():
        src = find_file(paper_dir, pattern)
        dst = os.path.join(case_input, std_name)

        if src is None:
            if std_name == 'nonrad_forcing.nc':
                print("  [skip] %s (optional, not found)" % std_name)
            else:
                print("  [WARN] %s not found matching '%s'" % (std_name, pattern))
            continue

        # Remove existing link/file
        if os.path.exists(dst) or os.path.islink(dst):
            os.remove(dst)

        if args.copy:
            import shutil
            shutil.copy2(src, dst)
            print("  [copy] %s -> %s" % (os.path.basename(src), std_name))
        else:
            os.symlink(src, dst)
            print("  [link] %s -> %s" % (os.path.basename(src), std_name))

    print("\nDone. Run: python3 run_case.py %s" % args.case)


if __name__ == '__main__':
    main()
