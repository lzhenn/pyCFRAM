# CFRAM-A Algorithm Specification

**Source:** `OA/202604020900/CFRAM_RRTMG/` Fortran code analysis  
**Reference:** Wu et al. (2025), J. Climate, Eq. 1

## 1. Overview

CFRAM decomposes total temperature change ΔT into partial contributions from individual physical processes:

```
ΔT = (∂R/∂T)⁻¹ × [Σ Δ(S-R)_i + ΔQ_ATD + ΔQ_SRF]
```

The Fortran code implements the radiative part (9 terms). Non-radiative terms (ATD, SRF) are NOT in this code.

## 2. Call Flow

```
cfram_rrtmg.f90 (main)
├── Read binary inputs: t, q, o3, ts, ps, solin, ssrd, ssru, cloud, aerosol, co2
├── Compute derived quantities: albedo_sw = ssru/ssrd, zenith = solin/scon
├── Determine active layers: nlayer (plev < ps)
│
├── rad_driver(base)    → rad_1d_base, fd_1d_base, fu_1d_base
├── rad_driver(warm)    → rad_1d_warm, fd_1d_warm, fu_1d_warm  [ALL vars changed]
│
├── calc_drdt(base)     → drdt(nlayer+1, nlayer+1)  [Planck matrix, LW only]
├── drdt_atm_inv = inv(drdt(1:nlayer, 1:nlayer))
│
├── rad_driver(co2)     → rad_1d_co2      [only co2 changed]
├── rad_driver(q)       → rad_1d_q        [only humidity changed]
├── rad_driver(ts)      → rad_1d_ts       [only surface T changed]
├── rad_driver(o3)      → rad_1d_o3       [only ozone changed]
├── rad_driver(solar)   → rad_1d_solar    [only zenith angle changed]
├── rad_driver(albedo)  → rad_1d_albedo   [only SW albedo changed]
├── rad_driver(cloud)   → rad_1d_cloud    [only cloud changed]
├── rad_driver(aerosol) → rad_1d_aerosol  [only aerosol changed]
│
├── Compute forcing: frc_i = rad_perturbed - rad_base  (per layer)
│                    frc_i(sfc) = net_flux_perturbed(sfc) - net_flux_base(sfc)
│
├── Solve: dT_i = drdt_atm_inv × (-frc_i)   [for each of 9 terms]
│
└── Write outputs: frc_*.dat, dT_*.dat
```

## 3. Coordinate Systems

| System | Index 1 | Index nlayer | Notes |
|--------|---------|-------------|-------|
| ERA5 input | TOA (1 hPa) | Surface (1000 hPa) | Top-to-bottom |
| RRTMG internal | Surface | TOA | Bottom-to-top |
| Output | TOA | Surface | Back to ERA order |

**Transformation:** `rrtmg_idx = nlayer - era_idx + 1`

**Interface pressures:**
- `pint(0) = 0.5 hPa` (TOA boundary)
- `pint(i) = (plev(i) + plev(i+1)) / 2` for i=1..nlayer-1
- `pint(nlayer) = ps` (surface)

## 4. Key Physical Computations

### 4.1 Column Density (molecules/cm²)

```
Mmoist = (1 - qv) × 28.966 + qv × 18.016    [g/mol, moist air molecular weight]
coldry = (Δp × 1e3 × 6.02e23) / (1e2 × 9.8 × Mmoist × (1 + qv))
```

where qv = specific humidity converted to volume mixing ratio = q × 28.966/18.016

### 4.2 Gas Mixing Ratios (wkl array, 7 species)

| Index | Species | Source | Conversion |
|-------|---------|--------|------------|
| 1 | H₂O | ERA5 specific humidity | q × 28.966/18.016 → volume mixing ratio |
| 2 | CO₂ | co2.txt (ppmv) | × 1e-6 |
| 3 | O₃ | ERA5 ozone mass mixing ratio | × 28.966/48.0 |
| 4 | N₂O | hardcoded 0.28 ppmv | × 1e-6 |
| 5 | CO | not used (0) | — |
| 6 | CH₄ | hardcoded 1.6 ppmv | × 1e-6 |
| 7 | O₂ | constant 0.209 | volume fraction |

**Final column amount:** `wkl_col = coldry × mixing_ratio` [molecules/cm²]

### 4.3 Cloud Water Path (g/m²)

```
coef = (1/9.8) × 100 × 1000 = 10204.08
ciwp = coef × cldiwc × Δp    [where Δp in hPa, cldiwc in kg/kg]
clwp = coef × cldlwc × Δp
```

Cloud effective radii: ice = 5 μm, liquid = 20 μm (hardcoded)

### 4.4 Surface Albedo

```
albedo_sw = swup_surf / swdn_surf    (0 if swdn=0)
albedo_lw = 0.0                       (hardcoded)
```

### 4.5 Solar Zenith Angle

```
zenith = solin / scon    (cosine of zenith angle)
```

where solin = TOA incident solar radiation (W/m²), scon = 1360.98 W/m²

## 5. Radiative Flux Definitions

For each atmospheric layer i (1 to nlayer):

```
sw_base(i) = fds(i) + fus(i+1) - fds(i+1) - fus(i)    [SW convergence]
lw_base(i) = fdl(i) + ful(i+1) - fdl(i+1) - ful(i)    [LW convergence]
rad_base(i) = lw_base(i) + sw_base(i)                    [total convergence]
```

Flux indices: 1 = TOA boundary, nlayer+1 = surface boundary (in RRTMG bottom-to-top order, but code handles the reversal)

## 6. Forcing Computation

For atmospheric layers (1:nlayer):
```
frc_X(1:nlayer) = rad_perturbed_X(1:nlayer) - rad_base(1:nlayer)
```

For surface (nlayer+1):
```
frc_X(nlayer+1) = [fd_perturbed(sfc) - fu_perturbed(sfc)] - [fd_base(sfc) - fu_base(sfc)]
```

## 7. Planck Matrix (drdt)

**Dimension:** (nlayer+1) × (nlayer+1)  
**Method:** LW-only, +1K perturbation per layer

```
For each layer k = 1..nlayer:
    T_perturbed = T_base
    T_perturbed(k) += 1.0 K
    Call rad_driver_lw → lw_1k, fdl_1k, ful_1k
    drdt(k, 1:nlayer) = lw_1k - lw_base
    drdt(k, nlayer+1) = [fdl_1k(sfc) - ful_1k(sfc)] - [fdl_base(sfc) - ful_base(sfc)]

For surface:
    Ts_perturbed = Ts_base + 1.0 K
    Call rad_driver_lw → lw_1k, fdl_1k, ful_1k
    drdt(nlayer+1, 1:nlayer) = lw_1k - lw_base
    drdt(nlayer+1, nlayer+1) = [fdl_1k(sfc) - ful_1k(sfc)] - [fdl_base(sfc) - ful_base(sfc)]
```

**Indexing convention:** drdt(row, col) where row = radiative response location, col = T perturbation location

## 8. Temperature Response

```
drdt_atm_inv = inv(drdt(1:nlayer, 1:nlayer))   [atmospheric only]

For each layer i = 1..nlayer:
    dT_X(i) = Σ_j [-frc_X(j) × drdt_atm_inv(j, i)]
            = dot_product(-frc_X(1:nlayer), drdt_atm_inv(:, i))
```

Note: only frc(1:nlayer) is used in the dot product, NOT frc(nlayer+1). The surface forcing is excluded from the atmospheric temperature solve.

Surface dT is set to 0.0 in the output.

## 9. Aerosol Optical Properties (from Matlab)

**Input:** MERRA-2 inst3_3d_aer_Nv (kg/kg mixing ratios on 72 model levels)  
**Lookup tables:** opticsBands_*.nc (GOCART, bext/bsca/g vs radius and RH, per RRTMG band)

```
For each species (BC, OC, SU, SS, DU) and each size bin:
    1. Read bext(band, rh_index, radius_index) from opticsBands_*.nc
    2. Find nearest RH index for each layer
    3. AOD(band, layer) = mixing_ratio × air_density × 1e3 × layer_thickness × bext × 1e-3
    4. SSA(band, layer) = bsca / bext
    5. g(band, layer) from lookup table
```

**Bands:** SW 1-14 (jpband), LW 15-30 (nbndlw)

## 10. Differences from Wu et al. (2025) Paper

| Paper (10 terms) | Code (9 terms) | Note |
|-------------------|----------------|------|
| ΔS^solar | frc_solar | ✓ |
| Δ(S-R)^CO₂ | frc_co2 | ✓ |
| Δ(S-R)^CH₄ | — | **Missing**: CH₄ hardcoded, no separate perturbation |
| Δ(S-R)^O₃ | frc_o3 | ✓ |
| Δ(S-R)^WV | frc_q | ✓ |
| Δ(S-R)^CLD | frc_cloud | ✓ |
| Δ(S-R)^AER | frc_aerosol | ✓ (total aerosol; paper also decomposes by species) |
| ΔS^α | frc_albedo | ✓ |
| ΔQ^ATD | — | **Missing**: atmospheric dynamics, not in this code |
| ΔQ^SRF | — | **Missing**: surface processes, not in this code |
| — | frc_ts | Extra: surface temperature perturbation (not in paper Eq.1) |

## 11. NCL Preprocessing Notes

- **Base period:** 1986-1995 (code default; paper uses 2003-2022 climatology)
- **Warm period:** 2004-2013 (code default; paper uses specific event year)
- **Grid:** 1°×1° global (code default; paper uses 0.25°×0.25° regional)
- **Solar radiation unit conversion:** ERA5 J/m² → W/m² (÷ 86400)
- **Bug in cloud.ncl/cldliq.ncl:** warm case output uses base case index (0 instead of 1)
