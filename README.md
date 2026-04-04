# pyCFRAM

Python-driven CFRAM (Climate Feedback-Response Analysis Method) framework using Fortran RRTMG as the radiation engine.

Reproduces the temperature decomposition analysis from Wu et al. (2025, *J. Climate*) for extreme heat events EH13 (Aug 2013) and EH22 (Aug 2022) over East China.

## Architecture

```
Fortran RRTMG (radiation engine)
  └── 11 RRTMG calls per grid point → 9 radiative forcings + Planck matrix (∂R/∂T)

Python (decomposition + analysis)
  ├── dT_i = -(∂R/∂T)⁻¹ × frc_i          # radiative terms (WV, cloud, aerosol, ...)
  ├── dT_lhflx/shflx/sfcdyn               # non-radiative via same Planck matrix
  ├── dT_atmdyn = dT_obs - Σ(all others)  # atmospheric dynamics (residual)
  └── multiprocessing parallel execution   # embarrassingly parallel over grid points
```

## Prerequisites

- Python 3.8+ with: `numpy`, `netCDF4`, `matplotlib`, `cartopy`, `scipy`
- `gfortran` (for compiling RRTMG)
- LAPACK/BLAS libraries

## Setup

### 1. Clone

```bash
git clone git@github.com:lzhenn/pyCFRAM.git
cd pyCFRAM
```

### 2. Compile Fortran RRTMG

```bash
cd fortran
make
cd ..
```

This builds the `cfram_rrtmg` executable. The RRTMG lookup tables (`rrtmg_lw.nc`, `rrtmg_sw.nc`) are tracked via git-lfs and must be present after clone.

### 3. Obtain paper_data

Download the author's validation dataset and place it under `pyCFRAM/paper_data/`:

```
paper_data/
└── cfram_out/
    ├── case_eh13_c20250102/
    │   ├── merra2_eh13_baseline_pres_input_check.nc
    │   ├── merra2_eh13_baseline_surf_input_check.nc
    │   ├── merra2_eh13_all_pres_input_check.nc
    │   ├── merra2_eh13_all_surf_input_check.nc
    │   ├── merra2_eh13_partial_forcing.nc
    │   ├── merra2_eh13_partial_t.nc
    │   └── ...
    └── case_eh22_c20250118/
        └── (same structure)
```

## Quickstart: Reproduce Fig.3

### Step 1: Extract full-field input

```bash
python3 scripts/extract_full_field.py --case eh13
python3 scripts/extract_full_field.py --case eh22
```

Reads `paper_data/` input files, computes aerosol optical properties, writes Fortran binary to `fortran/data_prep/`.

### Step 2: Run CFRAM decomposition (parallel)

```bash
python3 scripts/run_parallel_python.py --case eh13 --nproc 40
python3 scripts/run_parallel_python.py --case eh22 --nproc 40
```

9801 grid points (81 lat x 121 lon), each running an independent single-column RRTMG calculation. On 40 cores, each case takes ~20 minutes.

Output: `cfram_output/cfram_eh13_python.nc`, `cfram_eh22_python.nc`

### Step 3: Plot

```bash
python3 scripts/plot_fig3_self.py
```

Output: `figures/fig3_self_computed.png`

## Project Structure

```
pyCFRAM/
├── core/           # Python modules: radiation prep, aerosol optics,
│                   #   Planck matrix, decomposition, Fortran I/O
├── fortran/        # Modified CFRAM-RRTMG + upstream RRTMG LW/SW
├── scripts/        # Workflow scripts (extract → run → plot)
├── plotting/       # Matplotlib/Cartopy visualization
├── data/           # ERA5 data loader
├── tests/          # Unit tests
├── configs/        # Experiment parameters (defaults.yaml)
└── docs/           # Algorithm specification
```

## Key Scripts

| Script | Purpose |
|--------|---------|
| `extract_full_field.py` | Paper_data NetCDF → Fortran binary input |
| `run_parallel_python.py` | Parallel CFRAM decomposition (multiprocessing) |
| `plot_fig3_self.py` | Fig.3: 6-row spatial decomposition maps |
| `plot_fig3.py` | Fig.3 directly from paper_data (validation) |
| `plot_fig4.py` | Fig.4: PAP bar charts |
| `extract_paper_data_column.py` | Single-column extraction for debugging |

## Reference

Wu et al. (2025), *Journal of Climate*.
