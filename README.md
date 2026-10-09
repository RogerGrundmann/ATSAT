# ATSAT

Saturn atmospheric general circulation model based on the numerical framework of the ATOM
climate model. Solves the 3-D Navier-Stokes equations in a spherical shell extending from
the deep water-cloud level (~10 bar) up through the hydrocarbon haze layers, including full
cloud microphysics and chemistry for the H₂/He atmosphere of Saturn.

Saturn presents two phenomena that motivated this work. The first is the persistent
hexagonal polar vortex observed in Saturn's northern hemisphere by Voyager and Cassini —
a six-sided standing wave pattern surrounding the pole that has remained stable for
decades. The second is Saturn's super-rotating equatorial jet, which reaches wind speeds
near 470 m/s at the cloud tops, far exceeding anything observed on Jupiter. Both features
suggest deep-seated dynamical structures that a 3-D circulation model can illuminate.

For all relevant data concerning Saturn and its atmosphere the book *Planetary Sciences*
by Imke de Pater and Jack J. Lissauer was indispensable.

**Where this model sits in the literature.** Roughly half a dozen 3-D GCMs have been run on
Saturn, of which DYNAMICO-Saturn (Spiga, Guerlet, Cabanes *et al.*, Icarus Parts I–V) is the
only sustained Saturn-dedicated programme. Those papers state that the model **does not
parametrise moist convection**, and note that it "occurs in giant planets' atmospheres and
carries substantial energy into the jet system". That is the gap ATSAT is built for: a
three-category microphysics scheme for H₂O, NH₃ and CH₄ plus NH₄SH sedimentation, coupled
into the temperature and moisture equations. The moist physics is the reason for the code,
not the dynamical core.

---

## Physics & Numerics

- **Domain:** spherical shell, 41 × 181 × 361 grid points (r × θ × φ), ~2.7 million cells
- **Vertical extent:** 500 km atmospheric shell (12.5 km per radial level)
- **Dynamics:** finite-difference discretisation of the 3-D Navier-Stokes equations in spherical coordinates
- **Time integration:** 4th-order Runge-Kutta (inner loop) with a Poisson pressure solver (outer loop)
- **Parallelism:** OpenMP shared-memory threading across the chemistry, diffusion, and saturation-adjustment hot paths
- **Thermodynamics:**
  - Temperature initialised as a parabolic pole-to-pole profile (zonally uniform); reference temperature 134 K
  - Boussinesq buoyancy, carried as an **anomaly** (horizontal mean removed at each level)
  - Clausius-Clapeyron / Sanchez-Lavega SVP formulation for saturation vapour pressures
  - Mixed-phase (liquid + ice) saturation adjustment with iterative convergence
- **Microphysics:** COSMO-derived **three-category** scheme (rain / snow / graupel) run per
  condensable species, plus Stokes sedimentation for NH₄SH crystals
- **Boundary conditions:** measured zonal wind profiles from QuikSCAT and OSCAR datasets;
  temperature/pressure profiles from Voyager (1980/81), Cassini (2004–2017), and Hubble observations
- **Planetary constants:** g = 10.0 m/s², Ω = 1.63 × 10⁻⁴ rad/s, lapse rate 0.7 K/km

### Body forces

The Coriolis, centrifugal and buoyancy terms each need a different factor to reach the
non-dimensional unit system the RHS works in (`ATSAT_NONDIM` switches all three on together;
the individual `ATSAT_ND_*` switches exist for attribution, not for production runs). All three
read the **true** sine of colatitude, reconstructed from the cosine, never the metric's polar
floor — a floor is a property of the difference scheme and has no business in a body force.

---

## Chemical Species

| Species | Phases modelled | Precipitates as |
|---------|-----------------|-----------------|
| H₂O | vapour · cloud · ice | rain · snow · graupel |
| NH₃ | vapour · cloud · ice | rain · snow · graupel |
| CH₄ | vapour · cloud · ice | rain · snow · graupel |
| H₂S | vapour | — (consumed by the NH₄SH reaction) |
| NH₄SH | solid | Stokes sedimentation |

Methane condensation, absent from the warmer Jovian troposphere, becomes a first-order
process at Saturn's lower temperatures: the model reaches 61 K at the top of the shell,
well below methane's 90.69 K triple point, and produces a genuine CH₄ ice deck and methane
precipitation reaching the ground.

---

## Modules

Ported from ATJUP, the sister Jupiter model, which shares this numerical framework. Each is
gated by an environment knob, so switching one on is a measurable change rather than a new
baseline.

| module | file | knob | default |
|---|---|---|---|
| Precipitation microphysics | `PrecipitationSat.h` | `ATSAT_PRECIP` | **on** |
| └ coupling into rhs_t and the moisture RHS | `RHS_Sat_Turb.cpp` | `ATSAT_PRECIP_COUPLING` | **1.0** |
| Mirrored saturation adjustment | `SaturationAdjustmentSat.*` | `ATSAT_SATADJ` | off |
| Grey multi-layer radiation | `RadiationSat.h` | `ATSAT_RADIATION` | off |
| └ coupling into rhs_t | | `ATSAT_RAD_COUPLING` | 0.0 |
| Turbulence closure (k-ε / k-ω / k-ω SST) | `TurbulenceSat.h` | `ATSAT_TURB` | off |
| └ coupling into momentum and scalars | | `ATSAT_TURB_COUPLING` | 0.0 |
| Dry convective adjustment | `ConvectiveAdjustmentSat.h` | `ATSAT_CONV_ADJ` | off |
| Thermal-wind residual (diagnostic only) | `ThermalWindDiagSat.h` | `ATSAT_TW_DIAG` | off |

### Numerical safety nets

| | knob | default |
|---|---|---|
| Non-finite census, per field, with index extent | `ATSAT_NANCHECK` | off |
| SIGFPE on the first invalid operation | `ATSAT_FPE` | off |
| Temperature limiter (clip to the neighbour range) | `ATSAT_T_LIMITER` | off |
| Shapiro filter on u, v, w | `ATSAT_VEL_SHAPIRO_INLOOP` | 0 passes |
| Zero floor on the condensable species | `ATSAT_NO_CLAMP` to disable | **on** |
| Polar metric floor on sin θ | `ATSAT_SINTHE_MIN` | 0.0 (no floor) |
| Metric radius (see the note in `cSaturnModel.h`) | `ATSAT_METRIC_RADIUS` | 0.0 (off) |
| Timestep override | `ATSAT_DT` | formula (0.032 s / iteration) |

---

## Diagnostics

Printed at the `checkpoint` cadence:

- **Equatorial column profile** — `i`, `p[bar]`, `T[K]`, `eps`, `netRad`, `Q_rad` from the model top
  down to the deep boundary. ATSAT was the last model to gain it. It is the diagnostic that
  localises profile faults where a column *mean* cannot, and it bears directly on Saturn's open
  under-emission below: `eps` against `p[bar]` down the column is the direct view of the opacity that
  question is about.
- **Photosphere line** (with `ATSAT_RADIATION=1`) — mean OLR against the input budget, the τ=1 level
  in bar, and the temperature there beside the blackbody flux it implies. `T(τ=1) − T_eff(OLR)` is
  the *scheme's* excess; `T(τ=1) − T_eff(in)` is how far the *column* sits from the energy budget.
- **ParaView** — `Radiation`, `Q_rad_mW_m3` and `Emissivity` alongside the turbulence six and the
  precipitation eleven, in all four views.

- **printMinMax** — max/min with location for every prognostic and diagnostic field.
- **steadyQuery** (`ATSAT_STEADY`, default on) — for each prognostic field, max|f − f_n| over
  the grid with the cell it occurs in, i.e. the largest change that field underwent in one
  iteration, plus the largest residual of the continuity equation. The "how far from steady,
  and where" measurement. It runs *before* `restoreVar`, which is what makes the differences
  non-zero.
- **negative-value clamp budget** — gross mass clipped per species as a percentage of that
  field's own mass, **split by whether it occurred on the two radial boundary planes**. That
  split matters: `bcRadius` sets those planes by f = (4/3)f[a] − (1/3)f[b], which returns a
  negative value wherever a field decays steeply towards the boundary. That is extrapolation
  overshoot, not transport undershoot, and reading one for the other is expensive — `ch4_ice`
  once appeared to clip 1848 % of its own mass when 100 % of it sat on one plane the equations
  never touch.

---

## Repository Layout

```
ATSAT/
├── planet/          # core model (RHS, RK4, thermodynamics, chemistry, microphysics, I/O)
├── lib/             # array types, config parser, FFT, utilities
├── cli/             # command-line driver (sat)
├── python/          # Cython bindings (pyatsat)
├── saturn/          # run directory (XML config, observational data, output)
│   ├── oscar/       # OSCAR ocean-current dataset (used as zonal-wind template)
│   └── windspeed/   # QuikSCAT surface wind data
├── tinyxml2/        # vendored XML library
├── param.py         # code-generation script (auto-generates parameter files)
└── Makefile
```

`param.py` is the single source of truth for every configuration parameter: it generates
`SaturnParams.h.inc`, `SaturnLoadConfig.cpp.inc`, `cSaturnDefaults.cpp.inc`, the three
`config_atsat.xml` copies and the Cython accessors. **A parameter added to the C++ but not to
`param.py` is silently dropped the next time anything regenerates.**

---

## Build

**Dependencies:** C++11 compiler with OpenMP support, Python 3, Cython, NumPy.

```bash
make            # parameter files + CLI + Python extension
make sat        # CLI binary only
make python     # Python extension only
make clean
```

### Shared physics headers

Eleven headers in `planet/` are **byte-identical across ATJUP, ATSAT, ATNEPT and ATURAN**:

```
ATPhys.h  BoundaryConditions.h  ConvectiveAdjustment.h  FluxLimiter.h  ParaViewWriter.h
Precipitation.h  PressureSolver.h  Radiation.h  Reporting.h  SaturationAdjustment.h  Turbulence.h
```

There is no submodule and no symlink holding them together — Synology Drive has silently reverted
a working tree once, and a submodule costs friction on every clone. They are plain copies, and
`planet/SHARED.md5` is what makes a divergence loud:

```bash
make check-shared
```

Editing one means: edit it in one repo, copy it to the other three, regenerate its line in **all
four** manifests, and rebuild each.

**`make check-shared` cannot catch everything, and this is the part to read before trusting it.**
It verifies a repo against *its own* manifest, so two repos holding different copies of the same
header both report OK — a state that has already occurred once. The check that does catch it is a
diff between repos, which is why the checksum lines are kept sorted by filename:

```bash
diff ../ATJUP/planet/SHARED.md5 planet/SHARED.md5
```

---

## Usage

### Command-line — note the TWO arguments

```bash
./cli/sat <path> <config file name>
./cli/sat . config_atsat.xml            # config in the current directory
```

A single combined path fails with `couldn't load config file inside cSaturnModel`. This is the
reverse of ATOM's `cli/atm`.

```bash
OMP_NUM_THREADS=12 ./cli/sat . config_atsat.xml > run.log 2>&1
```

### Python

```python
import sys
sys.path.insert(0, "saturn")
import pyatsat

model = pyatsat.SaturnModel()
model.load_config("saturn/config_atsat.xml")
model.run()
```

---

## Output

With `checkpoint = 100` a run writes all three output kinds on the same iterations (100, 200, …)
into `<output_path>-Saturn/`:

| file | what |
|---|---|
| `Saturn_radial_20_<n>.vtk` | one radial level, all fields, plus the 2-D surface maps |
| `Saturn_zonal_180_<n>.vtk` | meridional slice at k = 180 |
| `Saturn_longal_112_<n>.vtk` | longitudinal slice at j = 112 |
| `Saturn_panorama_<n>.vts` | the full 3-D volume (~900 MB), needs `paraview_flag` |
| `sat_restart_<n>.bin` | binary restart, 37 prognostic arrays (~756 MB) |

Restart is controlled by `checkpoint_save_iter` (explicit dump) and `restart_from_iter`
(resume), both default −1. The periodic dump is written **only when the state is finite**, so a
diverged run can never overwrite a good restart point. The file magic is `"SAT1"`: ATJUP
restarts have identical grid dimensions and would otherwise be read as nonsense.

---

## Known limitations

None of these stops a run; all of them affect what a result means.

0. **`ATSAT_THERMAL_MASSFLUX` exists as a measurement instrument** (default 1.0, inert). It scales the
   diffusive-enthalpy sink in `rhs_t`. On this model that term is well-behaved — it *balances* the
   transport term (+0.163 against +0.158) and switching it off moves the photosphere by 0.24 K and
   the emitted flux by 0.6 % (1.234 → 1.226 W/m²). That is not true of the ice giants, where the same term dominates
   `rhs_t` by three to six orders of magnitude; ATSAT escapes because it builds its species
   diffusivities from per-species `mue_nh3/rg_nh3` and never consumes `mue_mix`. Worth knowing before
   the formula is "harmonised" across models.

1. **The microphysics rate coefficients are not calibrated to Saturn.** They are ATOM's
   terrestrial values rescaled *once* to **Jupiter's** energy budget and carried here unchanged,
   against a condensate loading ~4× higher. Where `P_snow` reads 8.6400 mm/day it is
   `P_max_flux`, the safety clamp, not the microphysics — at 250 km that is 65 337 of 65 341
   cells. Read those cells as an upper bound and the rest as physics. **This is the
   highest-value open item**, since the moist physics is the reason for the code; the observable
   is Lv·P against Saturn's emitted flux, the same one ATJUP used.
2. **The mirrored saturation adjustment is off** (`ATSAT_SATADJ`). Switching it on moves max
   `h2o_cloud` by −16 % and lifts the deck a layer. It carries a genuine repair — the legacy
   routine's saturation target was frozen. Until 2026-10-09 its ice-phase coefficients were the
   **liquid** ones; the model now has ice pairs for H₂O and NH₃ (the substance's, as on ATJUP;
   CH₄ has one curve), so that objection is gone. What is still owed before it becomes the
   default is a run of useful length comparing the two routines.
3. **`t_00_ch4` = 67.36 K is a borrowed proportion, not a measurement** — `t_0_ch4` scaled by
   the mean of the H₂O and NH₃ ratios (0.7427). It was 190.56 K, methane's *critical*
   temperature, which made the ice-autoconversion band `(T < 90.69 && T ≥ 190.56)` empty by
   construction and left methane ice with no sink at all.
4. **The timestep is tiny.** The default `dt` = 3.01e-5 is 0.032 s of Saturn time per iteration,
   so a 200-iteration run spans **6.4 seconds** and nothing horizontal can develop. Raising it
   is not simply available: `ATSAT_DT=0.001` goes non-finite at iteration 1 in the bottom 11
   levels, and that is *below* the advective CFL limit (~0.044), so whatever limits the deep
   atmosphere is not advection. Unresolved.
5. **`nh3_ice` clips 11.1 % of its own mass per 100 iterations with 0.0 % of it on the
   boundary** — genuine interior transport undershoot, the largest real one in either model.
6. **Not mirrored from ATJUP:** the per-level momentum census (`ATJUP_WPROFILE`), the
   single-cell term probes (`ATJUP_PROBE*`), and roughly 25 further knobs. Some are
   inapplicable — ATSAT has no obstacle, so the whole SeaMount family is meaningless here — but
   the momentum census in particular has repeatedly been the thing that settled an argument in
   ATJUP.

---

## Author

Roger Grundmann — roger.grundmann@web.de
