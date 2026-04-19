# pyCFRAM

A Python-encapsulated implementation of the **Climate Feedback–Response Analysis Method (CFRAM)**. pyCFRAM aims to build a user-friendly, high-efficient, and extendable interface for CFRAM users. It uses `multiprocessing` (and future `f2py`/`mpi4py` support) to exploit the maximum potential of code parallelism without dramatically modifying the original Fortran RRTMG radiation code.

Currently reproduces the temperature decomposition analysis from Wu et al. (2025, *J. Climate*) for extreme heat events EH13 (Aug 2013) and EH22 (Aug 2022) over East China.

For the CFRAM methodology, please refer to:
- Lu, J., and M. Cai, 2009: A new framework for isolating individual feedback processes in coupled general circulation climate models. *Climate Dynamics*.
- Cai, M., and J. Lu, 2009: A new framework for isolating individual feedback processes in coupled general circulation climate models. Part II. *Climate Dynamics*.

## Architecture

```
Fortran RRTMG (radiation engine, per grid point)
  ├── base + warm + 8 partial-perturbation rad_driver calls
  ├── Planck matrix (∂R/∂T) via nlayer+1 LW-only perturbations
  ├── 6 per-species aerosol perturbation calls (bc, ocphi, ocpho, sulf, ss, dust)
  └── cloud LW/SW component snapshot (zero extra RRTMG calls)
  → writes 9 bulk + 6 per-species + 2 cloud_lw/sw forcings + Planck inverse

Python (decomposition + analysis)
  ├── dT_i = -(∂R/∂T)⁻¹ × frc_i           # radiative terms (all 17, incl. splits)
  ├── dT_lhflx/shflx                       # non-radiative via same Planck matrix
  ├── dT_sfcdyn/ocndyn                     # energy-conservation residuals
  ├── dT_atmdyn = dT_obs − Σ(all others)   # atmospheric dynamics (residual)
  └── multiprocessing parallel execution   # embarrassingly parallel over grid points
```

## Prerequisites

- Python 3.8+ with: `numpy`, `netCDF4`, `matplotlib`, `cartopy`, `scipy`
- `gfortran` (for compiling RRTMG)
- LAPACK/BLAS libraries (on HPC clusters, may require `module load lapack`; see `docs/technical_notes_zh.md` §13.9)

## Setup

### 1. Clone and install dependencies

```bash
git clone git@github.com:lzhenn/pyCFRAM.git
cd pyCFRAM
pip install -r requirements.txt
```

### 2. Compile Fortran RRTMG

```bash
cd fortran
make          # builds cfram_rrtmg (full-field) and cfram_rrtmg_1col (single-column)
cd ..
```

The single-column executable `cfram_rrtmg_1col` is used by the parallel Python runner. The RRTMG lookup tables (`rrtmg_lw.nc`, `rrtmg_sw.nc`) are included in the repository.

### 3. Obtain source data (ERA5 + MERRA-2)

The standard pipeline drives CFRAM directly from ERA5 (state) and MERRA-2 (aerosol). Place these under `pyCFRAM/era5_data/`:

```
era5_data/
├── daily/                                 # ERA5 daily, 2003-2022 Aug
│   ├── era5_pl_{var}_{YYYY}08.nc          # PL: t, q, o3, cc, ciwc, clwc (6-hourly, 37 lev)
│   └── era5_sl_{YYYY}08/                  # SL: subdirectory per month
│       ├── *stepType-instant.nc           # skt, sp
│       ├── *stepType-accum.nc             # ssrd, ssr, tisr, slhf, sshf
│       └── *stepType-max.nc               # mx2t
└── merra2/
    └── M2I3NVAER_{YYYY}08/*.nc4           # 3-hourly, 72 model levels, 13 aerosol species
```

Helper downloaders (CDS / NASA Earthdata credentials required):

```bash
python3 scripts/download_era5_flux.py        # ERA5 PL+SL via CDS API
python3 scripts/download_merra2_aerosol.py   # MERRA-2 M2I3NVAER via GES DISC
```

Optional: to validate against Wu et al. results, also download `paper_data/` (not required for standard runs).

## Quickstart: Reproduce Fig.3

### Step 1: Build case input from ERA5 + MERRA-2

```bash
python3 scripts/build_case_input.py --case eh13
python3 scripts/build_case_input.py --case eh22
```

Reads `cases/<case>/case.yaml` → generates `cases/<case>/input/{base,perturbed}_{pres,surf}.nc` and `nonrad_forcing.nc`.

### Step 2: Run CFRAM (one command per case)

```bash
python3 run_case.py eh13
python3 run_case.py eh22
```

This extracts input, runs parallel CFRAM decomposition (all CPUs by default, ~20 min/case on 80 cores), and plots results. Output in `cases/eh13/output/` and `cases/eh13/figures/`.

Or run individual steps:

```bash
python3 run_case.py eh13 --step build     # only build ERA5+MERRA-2 input
python3 run_case.py eh13 --step extract   # only extract to Fortran binary
python3 run_case.py eh13 --step run       # only run CFRAM
python3 run_case.py eh13 --step plot      # only plot
```

### Validation against Wu et al. paper_data

```bash
python3 scripts/validate_vs_paper.py eh22           # surface dT comparison
python3 scripts/plot_fig3_independent.py eh22       # 2-column Fig.3 (paper | indep)
```

### Adding a new case

1. Create `cases/my_case/case.yaml` (see `cases/eh13/case.yaml` for template — includes `source:` block for ERA5/MERRA-2 driving)
2. Run: `python3 scripts/build_case_input.py --case my_case` to generate inputs
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
| `run_case.py` | Unified entry point: build → extract → run → plot |
| `build_case_input.py` | ERA5 + MERRA-2 → standard NetCDF input (driven by case.yaml) |
| `download_era5_flux.py` | Download ERA5 PL+SL data from Copernicus CDS |
| `download_merra2_aerosol.py` | Download MERRA-2 M2I3NVAER from NASA GES DISC |
| `extract_full_field.py` | Input NetCDF → Fortran binary (bulk + per-species aerosol optics) |
| `run_parallel_python.py` | Parallel CFRAM decomposition (multiprocessing) |
| `plot_fig3_independent.py` | Fig.3: 2-column (paper vs independent) spatial maps |
| `plot_fig3_self.py` | Fig.3: 6-row spatial decomposition maps |
| `plot_fig3.py` | Fig.3 directly from paper_data (validation) |
| `plot_fig4.py` | Fig.4: PAP bar charts |
| `plot_fig5.py` | Fig.5: per-species aerosol decomposition maps |
| `validate_vs_paper.py` | Compare surface dT against Wu et al. results |
| `verify_phase1_spc.py` | Sanity-check per-species optical-property additivity |
| `verify_phase4_spc.py` | Per-species aerosol forcing/dT additivity + regression |
| `verify_phase5_vs_paper.py` | Per-species spatial correlation vs Wu et al. paper_data |
| `verify_cloud_split.py` | Cloud LW/SW additivity (`bulk == lw + sw`) + regression |

## Decomposition Output

`cases/<case>/output/cfram_result.nc` contains partial temperature changes `dT_*` and the underlying radiative forcings `frc_*` (W/m²), shape `(lev, lat, lon)` with surface at `lev[-1]`:

| Group | Terms |
|-------|-------|
| Radiative (bulk) | `co2`, `q`, `ts`, `o3`, `solar`, `albedo`, `cloud`, `aerosol`, `warm` |
| Cloud LW/SW split | `cloud_lw`, `cloud_sw` — exact additive (`cloud == cloud_lw + cloud_sw`) |
| Aerosol species | `bc`, `ocphi`, `ocpho`, `sulf`, `ss`, `dust` — sum ≈ bulk (small non-linear residual) |
| Non-radiative | `lhflx`, `shflx` |
| Derived | `atmdyn`, `sfcdyn`, `ocndyn`, `observed` |

## Input Data Format

See [docs/input_spec.md](docs/input_spec.md) for the standard NetCDF input format. pyCFRAM takes any two atmospheric states (base, perturbed) — it does not prescribe how they are defined.

## References

**CFRAM methodology**
- Lu, J. and Cai, M., 2009. A new framework for isolating individual feedback processes in coupled general circulation climate models. Part I: Formulation. *Climate Dynamics*, 32, 873–885.
- Cai, M. and Lu, J., 2009. A new framework for isolating individual feedback processes in coupled general circulation climate models. Part II: Method demonstrations and comparisons. *Climate Dynamics*, 32, 887–900.

**CFRAM with aerosol (CFRAM-A framework origin)**
- Zhang, T., Deng, Y., Chen, J., Yang, S., Gao, P. and Zhang, H., 2022. Disentangling physical and dynamical drivers of the 2016/17 record-breaking warm winter in China. *Environmental Research Letters*, 17(7), 074024. https://doi.org/10.1088/1748-9326/ac79c1
  First CFRAM study to include an explicit aerosol radiative perturbation term: "the effect of aerosols ... has not been included in the previous CFRAM analysis" (Zhang et al. 2022). The work introduces the CFRAM-A variant with five-species aerosol decomposition (BC, OC, sulfate, sea salt, dust) driven by MERRA-2 reanalysis and off-line RRTMG v5 radiative transfer. pyCFRAM adopts the same conceptual framework (CFRAM-A, off-line RRTMG, MERRA-2 driving) and reimplements the pipeline with a Python orchestration layer: automated ERA5 + MERRA-2 preprocessing, Python-side full Planck-matrix solve, and `multiprocessing`-based per-grid-point parallelism. Note: the aerosol optical properties in Zhang et al. 2022 come from MAM4/CAM6; pyCFRAM uses GOCART lookup tables (matching the Wu et al. 2025 reference code used for reproduction).

**CFRAM applications (extreme-event attribution)**
- Wu, Q., Li, Q., Zhang, T., Sun, X., Yang, S. and Hu, X., 2025. Quantitative attribution of 2013 and 2022 extremely hot events in China: insights from a climate feedback–response perspective. *Journal of Climate*, 38(17), pp.4331–4349.