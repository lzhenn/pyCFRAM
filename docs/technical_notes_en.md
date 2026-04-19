# pyCFRAM Technical Notes (concise English version)

This is a condensed counterpart of [`technical_notes_zh.md`](technical_notes_zh.md).
The Chinese version is authoritative for deep background and exhaustive
troubleshooting; this file covers the essentials a new user needs to set up,
run, and interpret pyCFRAM.

---

## 1. What pyCFRAM does

pyCFRAM is a Python-orchestrated implementation of the **Climate Feedback–
Response Analysis Method (CFRAM)** that decomposes an observed temperature
change into partial contributions from individual physical processes:

```
ΔT_i = −(∂R/∂T)⁻¹ · F_i          (per process i)
ΔT_total ≈ Σ_i ΔT_i              (linearised)
```

`F_i` is the radiative forcing (W m⁻²) for process *i* obtained by perturbing
one variable at a time through RRTMG. `∂R/∂T` is the Planck response matrix,
also computed by RRTMG via +1 K layer-by-layer perturbations. The matrix
inverse is applied in Python.

**References**: Lu & Cai 2009 / Cai & Lu 2009 (CFRAM formulation);
Zhang et al. 2022 (first CFRAM study to include aerosol, "CFRAM-A" framework);
Wu et al. 2025 (the extreme-heat-event application reproduced here).

---

## 2. Architecture

```
Fortran RRTMG engine (per grid point, fortran/cfram_rrtmg_1col)
  ├── base + perturbed + 8 partial-perturbation rad_driver calls
  ├── Planck matrix (∂R/∂T) via nlayer+1 LW-only perturbations
  ├── 6 per-species aerosol perturbation calls
  └── cloud LW/SW snapshot (zero extra RRTMG invocations)
  → writes 17 forcings (9 bulk + 2 cloud_lw/sw + 6 species) + drdt_inv

Python (scripts/run_parallel_python.py)
  ├── dT = −(∂R/∂T)⁻¹ · F for every radiative forcing
  ├── dT for non-radiative lhflx/shflx via the same Planck matrix
  ├── sfcdyn/ocndyn derived from energy-balance residuals
  ├── atmdyn residual: dT_obs − Σ(all other terms)
  └── multiprocessing.Pool, one worker per grid point
```

**Why no OpenMP inside Fortran?** RRTMG has `SAVE` state in its absorption-
coefficient modules that is not thread-safe. Process-level isolation via
`multiprocessing` avoids the race cleanly. Each worker writes its own
`tmpdir/data_prep/`, invokes the single-column executable, and reads back
the forcings.

---

## 3. Install

Prerequisites:
- Python 3.8+ with `numpy`, `netCDF4`, `matplotlib`, `cartopy`, `scipy`, `pyyaml`
- `gfortran` ≥ 9 (or `ifort`, see `fortran/makefile.ifort`)
- LAPACK / BLAS — **on HPC clusters you must load the module** (e.g.
  `module load lapack`). A silent-exit symptom (absurdly high `pts/s`
  throughput and all-NaN results) almost always means `ldd cfram_rrtmg_1col`
  shows `not found`. See `technical_notes_zh.md` §13.9 for the full diagnosis.
- Accounts: Copernicus CDS (ERA5) + NASA Earthdata (MERRA-2) if you want to
  run the end-to-end pipeline; not needed to read pre-built cases.

Build:

```bash
cd fortran && make && cd ..   # produces cfram_rrtmg and cfram_rrtmg_1col
```

---

## 4. End-to-end run (EH22 as example)

```bash
python3 scripts/download_era5_flux.py          # first time only, ~5 GB
python3 scripts/download_merra2_aerosol.py     # first time only, ~35 GB
python3 scripts/build_case_input.py --case eh22    # ERA5+MERRA-2 → NetCDF
python3 run_case.py eh22                           # extract → run → plot
```

`run_case.py` chains `build → extract → run → plot`; use `--step <name>`
to run a single stage. A full `--step run` on a 256-core host with 80
workers finishes one EH case in about 10 minutes (~20 pts/s).

---

## 5. Input format

Per case, under `cases/<case>/input/`:

| File | Contents |
|---|---|
| `base_pres.nc` / `perturbed_pres.nc` | (time, lev, lat, lon) 3D fields: `ta_lay, q, o3, camt, cliq, cice, co2, bc, ocphi, ocpho, sulf, ss, dust` |
| `base_surf.nc` / `perturbed_surf.nc` | (time, lat, lon) 2D fields: `ts, ps, solar, albedo` |
| `nonrad_forcing.nc` (optional) | 3D ΔF: `lhflx, shflx` — `sfcdyn/ocndyn/atmdyn` are computed internally |

`lev` is surface → TOA (e.g. 1000, 975, …, 1 hPa). pyCFRAM flips internally
to RRTMG's bottom-to-top order. Full reference: [`input_spec.md`](input_spec.md).

---

## 6. Output variables in `cfram_result.nc`

Shape `(lev, lat, lon)` with surface at `lev[-1]`. Fill value `-999`.

| Group | Variables | Notes |
|---|---|---|
| Radiative (bulk) | `frc_/dT_{co2, q, ts, o3, solar, albedo, cloud, aerosol, warm}` | 9 terms |
| Cloud LW/SW split | `frc_/dT_{cloud_lw, cloud_sw}` | Exactly additive: `cloud = lw + sw` to float64 precision |
| Aerosol species | `frc_/dT_{bc, ocphi, ocpho, sulf, ss, dust}` | Sum ≈ bulk; small residual = non-linear coupling |
| Non-radiative | `dT_{lhflx, shflx}` | From the input `nonrad_forcing.nc` |
| Derived | `dT_{atmdyn, sfcdyn, ocndyn, observed}` | Residual + energy-balance splits |

---

## 7. Key algorithmic choices

**Per-species aerosol via RRTMG (path A), not external pre-computed forcings (path B).** For each species the Fortran engine builds a mixed state
where only that species is swapped base → warm, using AOD-weighted
combination for SSA and *g*:

```
tau_mix      = tau_bulk_base − tau_base_spc(i) + tau_warm_spc(i)
tau·ssa_mix  = (Σ tau·ssa)_base − (tau·ssa)_base_spc(i) + (tau·ssa)_warm_spc(i)
tau·ssa·g_mix = same with g
ssa_mix = tau·ssa_mix / tau_mix
g_mix   = tau·ssa·g_mix / tau·ssa_mix
```

This retains aerosol–aerosol and aerosol–cloud non-linear coupling that a
linear (path B) Planck-inverse split would discard. The price is a 2–5 %
additivity residual between `Σspecies` and the bulk aerosol forcing.

**Cloud LW/SW split — no extra RRTMG calls.** The LW and SW flux
components of the cloud perturbation are snapshotted from the shared
workspace immediately after the cloud `rad_driver` call. The resulting
`dT_cloud_lw + dT_cloud_sw` matches `dT_cloud` to float64 precision
(residual ~10⁻¹⁵).

**Surface `dT` in Python.** Unlike the reference Fortran code, which zeroed
surface `dT` by solving only the atmospheric sub-block, pyCFRAM inverts the
full `(nlayer+1) × (nlayer+1)` Planck matrix in Python and produces a
physically meaningful surface response.

**ERA5 6-hourly accumulation ×6.** CDS returns 1-hour increments per
6-hour step for `ssrd/ssr/tisr/slhf/sshf`. `data/era5_source.py` multiplies
by 6 before dividing by 86400. Missing this factor yields a 6× undercount
of `lhflx/shflx` and flips the sign of `dT_cloud`.

---

## 8. Known limitations

- **Dust species correlation vs paper** is poor (|corr| < 0.25 on EH22 /
  EH13). Two factors: dust surface `dT` magnitudes are tiny (~0.01–0.1 K),
  and the GOCART single effective-size approximation loses skill in that
  small-signal regime. Future work could use the five DU001–DU005 size
  bins independently.
- **Sulfate bias is systematically high** (~+0.2 K). Single effective
  radius of 0.16 μm for SO₄ optics may not match Wu et al.'s treatment.
- **Cloud / aerosol complementary bias**: path A preserves non-linear
  optical coupling, paper's path B does not. They differ by 10–30 % on
  individual terms but the `cloud + aerosol` sum correlation stays
  above 0.95.
- **`core/radiation.py`, `core/decomposition.py`, and
  `tests/test_radiation_prep.py`** are stale (import constants that no
  longer exist). The actual pipeline doesn't use them — they are
  pending cleanup.

---

## 9. Three things that will bite you

Distilled from `technical_notes_zh.md` §13:

1. **HPC LAPACK/BLAS not loaded → silent Fortran exit → NaN results + 1000+
   pts/s throughput.** Check `ldd fortran/cfram_rrtmg_1col`.
2. **ERA5 accumulation ×6 factor.** Diagnose via `build_case_input.py`
   log line `solar domain mean`. Correct value ≈ 420 W m⁻²; broken
   value ≈ 70 W m⁻².
3. **`ndarray.tofile()` is always C-order**, regardless of the array's
   Fortran flags. `np.asfortranarray(arr).tofile(...)` is a no-op for
   byte layout. The Fortran reader's implicit-DO loop order must match
   C-order disk bytes (innermost axis varies fastest on disk).

---

## 10. Remote operations (HKUST-specific)

Production runs happen on the `hqlx*` cluster:

- Edit locally → `rsync` to `mini` (macOS jump box) → `rsync` to
  `hqlx204` or similar. **Always start the sync from the Mac**; a direct
  `ssh mini 'rsync ...'` will use mini's possibly-stale copy.
- Compile on `hqlx127` (gfortran 4.8.5 is what the pre-built object files
  in `linux_gfortran_gfortran_dbl.obj/` were generated with); run on
  `hqlx204` (256 cores, `netCDF4` works).
- On `hqlx74` / `hqlx127` the system GLIBC is too old for modern `netCDF4`
  wheels; do only compilation there.

This is mentioned here because the checked-in compiled binaries
(`cfram_rrtmg`, `cfram_rrtmg_1col`) are ephemeral artefacts of that
workflow; a fresh clone on any system should run `cd fortran && make`.
