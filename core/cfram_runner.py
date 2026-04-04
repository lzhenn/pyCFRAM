"""Subprocess-based wrapper for the Fortran cfram_rrtmg executable.

Provides a Python interface to run the compiled Fortran CFRAM-A program
and read its outputs. This is the "Strategy B" approach: generate binary
input files -> run executable -> parse binary output files.

Usage:
    runner = CFRAMRunner('/path/to/CFRAM_RRTMG')
    results = runner.run()
    print(results['dT_co2'])  # partial temperature change from CO2
"""

import os
import subprocess
import numpy as np
from scipy.io import FortranFile


class CFRAMRunner:
    """Wrapper for the Fortran cfram_rrtmg executable."""

    # 9 CFRAM decomposition terms (matching Fortran output)
    TERMS = ['co2', 'q', 'ts', 'o3', 'solar', 'albedo', 'cloud', 'aerosol', 'warm']

    def __init__(self, work_dir, exe_name='cfram_rrtmg'):
        """Initialize runner.

        Args:
            work_dir: directory containing the executable and data_prep/
            exe_name: name of the Fortran executable
        """
        self.work_dir = os.path.abspath(work_dir)
        self.exe_path = os.path.join(self.work_dir, exe_name)
        self.data_prep = os.path.join(self.work_dir, 'data_prep')
        self.data_output = os.path.join(self.work_dir, 'data_output')

        if not os.path.isfile(self.exe_path):
            raise FileNotFoundError(f"Executable not found: {self.exe_path}")

    def run(self, timeout=300):
        """Execute cfram_rrtmg and return results.

        Args:
            timeout: max runtime in seconds

        Returns:
            dict with keys 'frc_<term>' and 'dT_<term>' for each term,
            each value is a 1D numpy array of length nlev+1.
        """
        os.makedirs(self.data_output, exist_ok=True)

        result = subprocess.run(
            [self.exe_path],
            cwd=self.work_dir,
            capture_output=True,
            text=True,
            timeout=timeout,
        )

        if result.returncode != 0:
            raise RuntimeError(
                f"cfram_rrtmg failed (rc={result.returncode}):\n"
                f"STDOUT: {result.stdout}\n"
                f"STDERR: {result.stderr}"
            )

        return self._read_outputs()

    def _read_outputs(self):
        """Read all output .dat files from data_output/."""
        results = {}

        for term in self.TERMS:
            for prefix in ('frc', 'dT'):
                fname = f'{prefix}_{term}.dat'
                fpath = os.path.join(self.data_output, fname)
                if os.path.exists(fpath):
                    data = self._read_sequential_unformatted(fpath)
                    results[f'{prefix}_{term}'] = data

        return results

    @staticmethod
    def _read_sequential_unformatted(filepath):
        """Read gfortran sequential unformatted binary file.

        gfortran writes: [8-byte record length][data][8-byte record length]
        """
        try:
            f = FortranFile(filepath, 'r')
            data = f.read_reals(dtype=np.float64)
            f.close()
            return data
        except Exception:
            # Fallback: try reading raw
            return np.fromfile(filepath, dtype=np.float64)

    def get_valid_layers(self, data, fill_value=-999.0):
        """Extract valid (non-fill) atmospheric layers from output.

        Args:
            data: 1D array of length nlev+1 (38)
            fill_value: fill value for inactive layers above surface

        Returns:
            valid_data: values for active layers (excluding fill values)
            nlayer: number of active layers
        """
        mask = data != fill_value
        return data[mask], np.sum(mask)

    def print_summary(self, results):
        """Print a summary of CFRAM results."""
        print(f"{'Term':<12} {'frc range':>24}  {'dT range':>24}")
        print("-" * 64)
        for term in self.TERMS:
            frc_key = f'frc_{term}'
            dt_key = f'dT_{term}'
            if frc_key in results and dt_key in results:
                frc = results[frc_key]
                dt = results[dt_key]
                frc_valid, _ = self.get_valid_layers(frc)
                dt_valid, _ = self.get_valid_layers(dt)
                if len(frc_valid) > 0 and len(dt_valid) > 0:
                    print(f"{term:<12} [{frc_valid.min():>10.4f}, {frc_valid.max():>10.4f}]  "
                          f"[{dt_valid.min():>10.4f}, {dt_valid.max():>10.4f}]")


if __name__ == '__main__':
    import sys
    work_dir = sys.argv[1] if len(sys.argv) > 1 else 'CFRAM_RRTMG'
    runner = CFRAMRunner(work_dir)
    print(f"Running CFRAM in {work_dir}...")
    results = runner.run()
    runner.print_summary(results)
