# pyCFRAM

A Python-encapsulated implementation of the **Climate Feedback–Response Analysis Method (CFRAM)**. pyCFRAM aims to build a user-friendly, high-efficient, and extendable interface for CFRAM users. It uses `multiprocessing` (and future `f2py`/`mpi4py` support) to exploit the maximum potential of code parallelism without dramatically modifying the original Fortran RRTMG radiation code.

Currently reproduces the temperature decomposition analysis from Wu et al. (2025, *J. Climate*) for extreme heat events EH13 (Aug 2013) and EH22 (Aug 2022) over East China.

For the CFRAM methodology, please refer to:
- Lu, J., and M. Cai, 2009: A new framework for isolating individual feedback processes in coupled general circulation climate models. *Climate Dynamics*.
- Cai, M., and J. Lu, 2009: A new framework for isolating individual feedback processes in coupled general circulation climate models. Part II. *Climate Dynamics*.

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

### Step 1: Prepare case input from paper_data

```bash
python3 scripts/prepare_from_paper_data.py \
    --paper_dir paper_data/cfram_out/case_eh13_c20250102 --case eh13
python3 scripts/prepare_from_paper_data.py \
    --paper_dir paper_data/cfram_out/case_eh22_c20250118 --case eh22
```

### Step 2: Run CFRAM (one command per case)

```bash
python3 run_case.py eh13
python3 run_case.py eh22
```

This extracts input, runs parallel CFRAM decomposition (all CPUs by default, ~20 min/case), and plots results. Output in `cases/eh13/output/` and `cases/eh13/figures/`.

Or run individual steps:

```bash
python3 run_case.py eh13 --step extract   # only extract input
python3 run_case.py eh13 --step run       # only run CFRAM
python3 run_case.py eh13 --step plot      # only plot
```

### Adding a new case

1. Create `cases/my_case/case.yaml` (see `cases/eh13/case.yaml` for template)
2. Place input NetCDF files in `cases/my_case/input/` (see `docs/input_spec.md`)
3. Run: `python3 run_case.py my_case`

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
