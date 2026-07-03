# kwline — Fully Relativistic Iron Line Model for Kerr-like Wormholes

## Overview

This project implements a custom XSPEC spectral model (`kwline` / `kwconv`) that computes relativistic iron line emission profiles from accretion disks around **Kerr-like wormholes**. Building upon the established Kerr black hole ray-tracing framework (YNOGK; Yang & Wang 2012), the code introduces a wormhole deformation parameter $\bar{\lambda}$ (`lambdaBar`) that smoothly interpolates between a standard Kerr black hole ($\bar{\lambda}=0$) and a traversable wormhole geometry. The model is designed for X-ray spectral fitting with NASA's XSPEC software, enabling observational constraints on wormhole spacetime parameters from iron K$\alpha$ line profiles.


## Code Architecture

### Source Files

| File | Description |
|------|-------------|
| `amodules.f90` | Core modules: `internal_grids`, `constants`, `rootsfinding`, `ellfunction`, `blcoordinate` (geodesic solver), `mod6_raytrace_camera`. Contains all ray-tracing machinery with wormhole corrections (~8300 lines). |
| `kerrwhline.f90` | Main subroutines `kwline` (line model) and `kwconv` (convolution model), plus a standalone test wrapper program. |
| `lmodel.dat` | XSPEC model definition file specifying `kwline` (9 params, additive) and `kwconv` (8 params, convolution). |
| `kerrwhlineFunctionMap.cxx` | C++ function map registering `kwline` and `kwconv` as XSPEC model functions. |
| `lpack_kerrwhline.cxx` | TCL package initializer for loading the model into XSPEC. |
| `Makefile` | Build system using HEASoft's `HD_*` variables for shared library compilation. |

### Key Modules in `amodules.f90`

- **`blcoordinate`**: Full geodesic integration with wormhole corrections. Key routines:
  - `lambdaq()` — Computes motion constants $\lambda$, $q$ using separated $\Delta_K$ and $\hat{\Delta}$
  - `YNOGK()` — Main geodesic solver (Weierstrass $\wp$-function based)
  - `INTRPART()` / `INTRPART_WH_CORR()` — Radial integrals with $\sqrt{\Delta_K/\hat{\Delta}}$ wormhole correction
  - `INTTPART()` — Angular integrals (unchanged from Kerr)
  - `radiustp()` — Radial turning points from Kerr quartic with throat boundary
  - `Pemdisk()` — Finds disk-photon intersections
  - `metricg()` — Evaluates metric components, separating $\Delta_K$ and $\hat{\Delta}$

### Data Flow

```
Parameters (a, inc, E_line, indices, radii, lambdaBar)
  └─→ impactgrid()        : Generate (α, β) impact parameter grid
      └─→ dGRtrace()      : Ray-trace each photon geodesic
          ├─→ lambdaq()    : Compute (λ, q) constants of motion
          ├─→ Pemdisk()    : Find disk crossing parameter p
          └─→ YNOGK()      : Integrate geodesic to find r_emission
              └─→ INTRPART() + WH_PATH_COMPENSATION()
                             : Apply wormhole radial correction
      └─→ dlgfac()        : Compute g-factor at each emission radius
      └─→ Accumulate line profile on energy grid
  └─→ FFT convolution (kwconv) or direct output (kwline)
```

## XSPEC Models

### `kwline` — Additive Line Model (9 parameters)

| # | Parameter | Unit | Default | Description |
|---|-----------|------|---------|-------------|
| 1 | `a` | — | 0.5 | Black hole / wormhole spin |
| 2 | `inc` | deg | 49 | Observer inclination angle |
| 3 | `Eline` | keV | 6.4 | Rest-frame line energy |
| 4 | `gamma_i` | — | 3.0 | Inner emissivity index |
| 5 | `gamma_o` | — | 3.0 | Outer emissivity index |
| 6 | `r_br` | $r_g$ | 15 | Break radius (negative → units of $\max(r_{\rm th}, r_{\rm ISCO})$) |
| 7 | `r_i` | $r_g$ | −1 | Inner disk radius (negative → units of $\max(r_{\rm th}, r_{\rm ISCO})$) |
| 8 | `r_o` | $r_g$ | 50 | Outer disk radius |
| 9 | `lambdaBar` | — | 0.0 | Wormhole deformation parameter |

### `kwconv` — Convolution Model (8 parameters)

Same as `kwline` but without `Eline`; designed to convolve with reflection spectra (e.g., `xillverCp`).

| # | Parameter | Unit | Default | Description |
|---|-----------|------|---------|-------------|
| 1 | `a` | — | 0.5 | Spin |
| 2 | `inc` | deg | 49 | Inclination |
| 3 | `index_i` | — | 3.0 | Inner emissivity index |
| 4 | `index_o` | — | 3.0 | Outer emissivity index |
| 5 | `r_br` | $r_g$ | 15 | Break radius |
| 6 | `r_i` | $r_g$ | −1 | Inner disk radius |
| 7 | `r_o` | $r_g$ | 50 | Outer disk radius |
| 8 | `lambdaBar` | — | 0.0 | Wormhole deformation parameter |

## Build & Run

### Prerequisites

- **HEASoft/XSPEC** (v6.30+) with development headers
- Fortran compiler (gfortran) and C++ compiler

### Compile the XSPEC Model

```bash
# Source HEASoft environment
source $HEADAS/headas-init.sh

# Build shared library
make
```

### Load in XSPEC

```tcl
# Method 1: Using load script
@load.xcm

# Method 2: Manual
initpackage kerrwhline lmodel.dat kerrwhline.f90
lmod kerrwhline .
```

### Fit Example

```tcl
# Load data and model
data 1 fakedata_lambda09_50ks_nobkg.pi
lmod kerrwhline .

# Iron line only
model kwline

# Full reflection spectrum convolution
model kwconv*atable{xillverCp_v3.6.fits}

# Fit
statistic chi
ignore **:**-3., 79.0-**
fit 1000
plot ldata delchi
```

## Output Format

Spectral output files (`.dat`) contain two columns:
- Column 1: Energy $E$ (keV), bin center
- Column 2: $E^2 \times \text{photar} / \Delta E$ (normalized flux density)


## Availability & Citation

The complete source code for **kwline** and **kwconv** is available in this repository.

If you use this software or any part of its implementation in your research, please cite:

> C. Liu, H. Siew, H.-X. Jiang, Y. Mizuno, and T. Zhu, *Signature of iron line profile from a Kerr-like wormhole*, **Astronomy & Astrophysics** (2026), DOI: 10.1051/0004-6361/202557027.

This work builds upon the YNOGK framework (Yang & Wang 2012) and the iron-line fitting approach of Mummery & Ingram (2024). We also gratefully acknowledge Andrew Mummery for his invaluable guidance and correspondence during the development of this work.

