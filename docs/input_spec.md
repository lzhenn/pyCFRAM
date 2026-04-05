# pyCFRAM Input Data Specification

pyCFRAM requires two atmospheric states (**base** and **perturbed**) as input. Each state is represented by two NetCDF files: one for pressure-level (3D) variables and one for surface (2D) variables.

Optionally, a third file provides non-radiative forcing for surface process decomposition.

## File Structure

```
cases/<case_name>/input/
├── base_pres.nc           # Base state: pressure-level variables
├── base_surf.nc           # Base state: surface variables
├── perturbed_pres.nc      # Perturbed state: pressure-level variables
├── perturbed_surf.nc      # Perturbed state: surface variables
└── nonrad_forcing.nc      # Non-radiative forcing (optional)
```

## Pressure-level file (`*_pres.nc`)

### Dimensions

| Dimension | Description |
|-----------|-------------|
| `time` | Time (typically 1 for mean state) |
| `lev` | Pressure levels (hPa), ordered surface → TOA |
| `lat` | Latitude (degrees north) |
| `lon` | Longitude (degrees east) |

### Required Variables

| Variable | Dimensions | Units | Description |
|----------|------------|-------|-------------|
| `ta_lay` | (time, lev, lat, lon) | K | Layer-mean temperature |
| `q` | (time, lev, lat, lon) | kg/kg | Specific humidity |
| `o3` | (time, lev, lat, lon) | kg/kg | Ozone mass mixing ratio |
| `camt` | (time, lev, lat, lon) | 0–1 | Cloud fraction |
| `cliq` | (time, lev, lat, lon) | kg/kg | Cloud liquid water content |
| `cice` | (time, lev, lat, lon) | kg/kg | Cloud ice water content |
| `co2` | (time, lev, lat, lon) | mol/mol | CO2 volume mixing ratio |
| `bc` | (time, lev, lat, lon) | kg/kg | Black carbon mixing ratio |
| `ocphi` | (time, lev, lat, lon) | kg/kg | Hydrophilic organic carbon |
| `ocpho` | (time, lev, lat, lon) | kg/kg | Hydrophobic organic carbon |
| `sulf` | (time, lev, lat, lon) | kg/kg | Sulfate mixing ratio |
| `ss` | (time, lev, lat, lon) | kg/kg | Sea salt mixing ratio |
| `dust` | (time, lev, lat, lon) | kg/kg | Dust mixing ratio |

### Notes

- `lev` must be ordered **surface → TOA** (e.g., 1000, 975, ..., 3, 2, 1 hPa)
- pyCFRAM internally flips to TOA → surface for Fortran RRTMG
- The 37 standard ERA5 pressure levels are: 1, 2, 3, 5, 7, 10, 20, 30, 50, 70, 100, 125, 150, 175, 200, 225, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000 hPa

## Surface file (`*_surf.nc`)

### Required Variables

| Variable | Dimensions | Units | Description |
|----------|------------|-------|-------------|
| `ts` | (time, lat, lon) | K | Skin/surface temperature |
| `ps` | (time, lat, lon) | Pa | Surface pressure |
| `solar` | (time, lat, lon) | W/m² | TOA incident solar radiation |
| `albedo` | (time, lat, lon) | 0–1 | Surface albedo |

## Non-radiative forcing file (`nonrad_forcing.nc`, optional)

If provided, pyCFRAM uses the Planck matrix from RRTMG to convert these forcings into temperature contributions.

### Variables

| Variable | Dimensions | Units | Description |
|----------|------------|-------|-------------|
| `lhflx` | (lev, lat, lon) | W/m² | Latent heat flux forcing (perturbed − base) |
| `shflx` | (lev, lat, lon) | W/m² | Sensible heat flux forcing |
| `sfcdyn` | (lev, lat, lon) | W/m² | Surface dynamics forcing |

### Notes

- Forcing is defined as: `F = flux_perturbed - flux_base`
- For surface-only fluxes (LH, SH), typically only the surface level is nonzero
- If this file is absent, pyCFRAM outputs only radiative decomposition terms

## Data Sources

These input files can be constructed from:

- **ERA5 / MERRA-2 reanalysis**: Average over desired periods to create base and perturbed states
- **CMIP6 model output**: Use historical vs. future scenarios
- **Sensitivity experiments**: Control vs. perturbed runs
- **Author-provided paper_data**: Use `scripts/prepare_from_paper_data.py` to convert

pyCFRAM does **not** prescribe how the two states are defined — any two atmospheric states can be compared.
