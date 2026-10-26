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
| Shared saturation adjustment (0 = the inherited routine) | `SaturationAdjustmentSat.*` | `ATSAT_SATADJ` | **on** |
| Grey multi-layer radiation | `RadiationSat.h` | `ATSAT_RADIATION` | **on** since 2026-10-11 (0 = off) |
| └ coupling into rhs_t | | `ATSAT_RAD_COUPLING` | 0.0 |
| Turbulence closure (k-ε / k-ω / k-ω SST) | `TurbulenceSat.h` | `turb_model` in the config (`ATSAT_TURB=0` switches it off) | **on** since 2026-10-11 |
| └ coupling into momentum and scalars | | `ATSAT_TURB_COUPLING` | 0.0 |
| Dry convective adjustment | `ConvectiveAdjustmentSat.h` | `ATSAT_CONV_ADJ` | off |
| Moist convective adjustment (Ledoux; H2O, NH3, CH4; H2S / NH4SH as tracers), meaningful with `ATSAT_LOCAL_RHO=1` | shared `ConvectiveAdjustment.h` | `ATSAT_MCA=<stride>` | off |
| Fitted NH3 liquid / NH3 ice / H2O ice vapour pressures (NIST Antoine, Fray & Schmitt 2009, Murphy & Koop 2005) | shared `ATPhys.h` | `ATSAT_SVP_FIT` | off |
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
| Metric radius (see the note in `cSaturnModel.h`) | `ATSAT_METRIC_RADIUS` | 58232 km since 2026-10-09 (0 = the unit sphere of before) |
| Timestep | `ATSAT_DT` | 0.001 = 1.06 s / iteration since 2026-10-09 (−1 = the old formula, 0.032 s / iteration) |
| Pole boundary form | `ATSAT_BC_POLE_COPY` | 1 since 2026-10-09 (0 = extrapolate, as before) |
| Rows next to the pole | `ATSAT_BC_POLE_FILL` | 1 since 2026-10-11 (0 restores; sets `j = 1` and `j = jm-2`, which the Runge-Kutta step does not integrate and which otherwise keep the initial condition for the whole run, from `j = 2` / `jm-3`: copy for the scalars and `u`, half the value for `v`, `w`) |
| Lid: species | `ATSAT_BC_SPECIES_LID_COPY` | 1 since 2026-10-11 (0 = the cubic extrapolation, which overshot on the lid: CH4 ice 93 g/m3 there against 42 one cell below) |
| Bottom wall: species | `ATSAT_BC_SPECIES_WALL_COPY` | 1 since 2026-10-11 (0 = the cubic extrapolation; the copy closes the wall for the diffusion stencil, as on ATJUP. The species are flat at the bottom today, so no effect of it has been measured) |
| Lid: v and w | `ATSAT_BC_TOP_TAPER` | 2 since 2026-10-11 (0 = extrapolated, as before; 2 = tapered to zero over the top three layers) |
| Adjustment kept by the time step | `ATSAT_RK_SCALAR_SYNC` | 2 since 2026-10-11 (0 = before: the Runge-Kutta step started from copies made before the saturation adjustment and discarded its condensate and its latent heat every step; 1 = keep the condensate only) |
| Temperature budget by level (print-only) | `ATSAT_TB_DIAG=<n>` | off (prints every n iterations which stage of the loop and which term of `rhs_t` changed the level-mean temperature; `planet/TbSat.h`) |
| Continuity residuum in a floored geometry (print-only) | `ATSAT_RES_FLOOR=<s>` | off (print-only (2026-10-17, shared `planet/Reporting.h`): at every residuum printout a second line with the continuity residuum in a geometry whose sin θ is floored at s, the maximum and its cell, then the five parts (du/dr, 2u/r, dv/dθ / r, cot v / r, dw/dφ / (r sin)) at both reported cells in both geometries. The standard line uses the TRUE sine; a model that clamps the polar metric integrates cos θ / s instead of cot θ poleward of asin(s), 23× less at 89° for s = 0.4. ATSAT integrates the true sine (`sinthe_min()` = 0), so its standard line is already in the model's geometry. Output identical unset and set, satchk/rescheck_2026-10-17) |
| Residuum by part around the reported cell (print-only) | `ATSAT_RES_PARTS` | 0 (shared `planet/Reporting.h`, from ATURAN 2026-10-16; output identical with and without) |
| Column budget of H2O, NH3, CH4 (print-only) | `ATSAT_CWB_DIAG=<n>` | off (prints every n iterations which stage of the loop and which term of the species equations changed each column of vapour + cloud + ice, in g/m², with `q·div(u)` as a reference; `planet/CwbSat.h`, ported from `ATJUP_CWB_DIAG`) |
| Radial branch of the initial cells | `ATSAT_U_INIT_SCALE` | 0 since 2026-10-12 (1 = before: a prescribed radial wind of up to 80 m/s without a matching horizontal branch, which was the velocity divergence of the initial state; with it the H2O and NH3 columns lost 0.45 % in 224 iterations through the advective-form transport, without it 0.21 %) |
| Flux form of the species transport | `ATSAT_SPECIES_DIVU=<s>` | 0 = off (adds `−s·q·div(u)` to all eleven species; at 1, 224 iterations: H2O column +0.003 %, NH3 +0.007 %, but CH4 −0.019 % instead of +0.020 % — not closed at the lid) |
| Initial dynamic pressure | `ATSAT_PDYN_INIT` | 2 since 2026-10-12 (later the same day): `p_dyn` starts in balance with the centrifugal term `(v²+w²)/r` of the radial momentum equation, as on the other three models; 224 iterations against the zero start: radial wind −0.004..+0.142 → −0.014..+0.082 m/s, columns, continuity residuum (0.026) and precipitation unchanged. 0 = zero start (the default of that morning), with which the columns change by +0.0007 / +0.0049 / −0.0002 % against the kinetic-energy field (1 = before: the kinetic energy density of the prescribed wind, whose radial gradient accelerated u by up to 1.5 m/s per iteration under the lid; with 0 the H2O / NH3 / CH4 columns change by +0.0007 / +0.0049 / −0.0002 % in 224 iterations instead of −0.21 / −0.21 / +0.02 %, and the continuity residuum stays at 0.026) |
| Rigid lid for the radial wind | `ATSAT_BC_RIGID_LID` | 1 since 2026-10-12, as on ATJUP (0 = before: u on both radial planes extrapolated from the interior; the same knob sets the pressure solver's wall condition). 224 iterations, 0 against 1: columns and continuity residuum unchanged (0.026), radial wind −0.007..+0.140 → −0.004..+0.142 m/s |
| Periodic phi seam | `ATSAT_BC_PHI_PERIODIC` | 1 since 2026-10-13 (0 restores the reconstructed seam, byte-identical to the head before). The longitude grid carries the same meridian twice, at k = 0 and k = km-1, and no loop that differences in phi integrated it: `bcPhi` set it to the mean of two one-sided extrapolations, the pressure solver did the same to `p_dyn`, the Shapiro filter took the two copies for each other's neighbour, and the output routine reconstructed T at the seam in the state at every write. With 1 the meridian is an ordinary cell — Runge-Kutta, RHS, chemistry, latent heat, forces, pressure solver, turbulence closure, flux limiter and filter take their phi neighbours round the ring of the km-1 distinct meridians — and k = km-1 is its copy. Measured here at 16 iterations (ATJUP/satchk/giants, SPH0 against SPH16): the seam deviations from the zonal median (nue_t 65.8, massflux_h2s 0.083, NH3 8e-4) are gone, 63 of 100 printed extrema identical, the continuity residuum 0.0265 unchanged (it sits at 89° S) |
| Metric split | `ATSAT_METRIC_SPLIT` | 1 since 2026-10-13 (0 restores, byte-identical to the head before): the radial metric term (v²+w²)/r of `rhs_u` is carried by the radial integral of itself, rebuilt every step, instead of by the pressure solver — ATOM_Precipitation's hydrostatic split applied to the other radial body force. `rhs_u` drops the term (it and the radial gradient of that pressure cancel by construction), `rhs_v` and `rhs_w` take its horizontal gradient, and `p_dyn` starts at zero. Why: the solver cannot hold the gradient that balances a radial body force (measured on ATNEPT with `ATNEPT_UB_DIAG`: after 224 iterations the balanced initial `p_dyn` had lost nine tenths of its gradient next to the radial walls; residuum there 0.008808 → 0.000341 with the split). The pressure is the new array `p_metric`. Measured here at 16 iterations (ATJUP/satchk/giants, SPDd against SMDd): 81 of 100 printed extrema identical, max u 0.00395 → 0.00118, max `p_dyn` 0.00086 → 0.000006, continuity residuum 0.026454 → 0.026450 (it sits at 89° S and is not this). 224 iterations, 0 against 1 (SM0_224, SM1_224): continuity residuum 0.026152 → 0.026098, printed max u 0.0819 → 0.0112, max `p_dyn` 0.000665 → 0.000031, columns unchanged |
| Radial Coriolis term in the split | `ATSAT_CORIOLIS_SPLIT` | 1 since 2026-10-13 (0 restores, byte-identical to the head before; acts only with `ATSAT_METRIC_SPLIT`): the radial Coriolis term of `rhs_u` joins the metric term in the same column integral and leaves `rhs_u` — it is a radial body force too, and with the metric term gone it was what still drove u (ATNEPT, equatorial band: −0.014 m/s in 224 iterations at every level). It enters the integral with the true sine, not a metric floor. 224 iterations, 0 against 1 (ATJUP/satchk/giants/SM1_224, SCS224): continuity residuum 0.026098 in both (it sits at 89° S), printed max u 0.0112 → 0.0043, columns unchanged in the third digit |
| Geostrophically balanced start | `ATSAT_GEO_BALANCE` | **1 since 2026-10-25** (0 restores the earlier default: Coriolis 1064× too weak, i.e. Saturn nearly non-rotating; every ATSAT number before it was made that way). With it `ATSAT_CELL_HOLD` (follows it; the zonal Coriolis term on the departure from the initial zonal-mean cells): ATJUP/satchk/cor_2026-10-25 SG0_1000 vs SGH_1000, 1000 it: bottom polar-cell w drift +10.8 → 0.00 m/s, residuum flat 0.0250; the 6.2 m/s w left at i30 86–88° is vertical diffusion of the kinked prescribed polar jet (corr +0.97 with the 2nd radial difference), not the hold. Built 2026-10-17 with default 0 (lead (b)): ATJUP's port of the ice giants' balanced start. Switches on `ATSAT_ND_COR` 1 (Coriolis × L/u₀ = 1064), `ATSAT_PDYN_GEO` 1 (`initGeoBalance()` in InitVariables_Sat.cpp: p_dyn the barotropic balance of each latitude's jet-level zonal-mean wind in `rhs_v`, `p_baro` the frozen per-level rest in `p_metric`; p_dyn −3.39 … +1.90, p_baro −5.36 … 0.005, model unit), `ATSAT_PRESS_INCR` 1 (shared solver, defect-correction source) and `ATSAT_DAMP_PDYN` 0 (no Shapiro pass on p_dyn). WHY: with `ATSAT_ND_COR=1` alone v grows at exactly f·w (−4.10 m/s per 100 it at 87 S; ATJUP/satchk/cor_2026-10-17). 224 iterations, 6 threads (ATJUP/satchk/geo_2026-10-17/SGEO224b against cor_2026-10-17/SCF1000, SNC1000): largest v change against the default at 200, zonal slice 1.00 m/s (ND_COR alone 8.4), level 20 0.46 (5.83); residuum at 208 0.02628 (0.02609). The first try ran away under the lid (SGEO224: u −70 m/s at 33 N by 208, ×3.7 per 32 it): with PRESS_INCR the source read `rhs_u`, which `prepareProjectionBoundaries` extrapolates with the 3-point cubic at the radial walls while it zeroes `aux_u`; the shared solver now takes the rigid-lid zero there (attribution SG_NOINCR / SG_DAMP / SG_INCR / SG_NOBARO). u under the lid still grows linearly, 0.045 m/s at 208 |
| Radial-wind budget | `ATSAT_UB_DIAG=<n>` | off; print-only (2026-10-13, `planet/UbSat.h`, ported from ATNEPT): every n iterations, which stage of the loop and which term of `rhs_u` changed the level-mean radial wind, for the whole level and for the equatorial band (`ATSAT_UB_BAND=<degrees>`, default 5), with the radial pressure force and what the pressure solver and the filter on `p_dyn` did to it. This model has no filter on u in the loop and no sponge, so those two rows stay zero. Unset and set both leave every output file identical (16 iterations) |
| Precipitation rate scale | `ATSAT_PRECIP_SCALE` | 0.0045 since 2026-10-12 (0.006 before, calibrated on 10-09 on a model that discarded the condensate; on today's defaults it read 6.00 W/m² of latent flux against 4.45 emitted at 224 iterations, 0.0045 reads 4.45) |

---

## Diagnostics

Printed at the `checkpoint` cadence:

- **Equatorial column profile** — `i`, `p[bar]`, `T[K]`, `eps`, `netRad`, `Q_rad` from the model top
  down to the deep boundary. ATSAT was the last model to gain it. It is the diagnostic that
  localises profile faults where a column *mean* cannot, and it bears directly on Saturn's open
  under-emission below: `eps` against `p[bar]` down the column is the direct view of the opacity that
  question is about.
- **Photosphere line** (printed unless `ATSAT_RADIATION=0`) — mean OLR against the input budget, the τ=1 level
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

0. **`ATSAT_THERMAL_MASSFLUX` defaults to 0.0 since 2026-10-11: the diffusive-enthalpy term is out of the temperature equation**, as on ATJUP (1.0 restores it). With it the column drifted — over 224 iterations the level-mean temperature changed by +1.5 K at the bottom, −3.4 K at level 20 and +8.2 K under the lid, and the term itself supplied +1.5 / −2.5 / +7.9 of that; without it the changes are −0.06 / −1.05 / +0.17 K and max T stays at its initial value. `ATSAT_TB_DIAG` is the instrument that showed it. What follows describes the knob as it was introduced (default 1.0, inert). It scales the
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
2. **The shared saturation adjustment is the default since 2026-10-10** (`ATSAT_SATADJ`, 0 =
   the inherited routine). It carries a genuine repair — the legacy routine's saturation target
   was frozen — and real ice pairs for H₂O and NH₃. Over 224 iterations against the inherited
   routine: temperature, velocity and pressure extrema identical, max `h2o_cloud` −6.7 %, max
   `nh3_ice` +26 %, max `ch4_ice` −23 %; it deletes no ice and leaves no cell unconverged.
   Which of the two cloud amounts is right is not settled by a measurement.
3. **`t_00_ch4` = 67.36 K is a borrowed proportion, not a measurement** — `t_0_ch4` scaled by
   the mean of the H₂O and NH₃ ratios (0.7427). It was 190.56 K, methane's *critical*
   temperature, which made the ice-autoconversion band `(T < 90.69 && T ≥ 190.56)` empty by
   construction and left methane ice with no sink at all.
4. **The timestep was tiny until 2026-10-09** (3.01e-5, i.e. 0.032 s of Saturn time per
   iteration). What limited it was the unit sphere the model ran on: the zonal spacing next to
   the poles was ~116 times too small, so `ATSAT_DT=0.001` went non-finite at iteration 3, first
   at 88–90° in the lowest levels. With Saturn's radius in the metric (`ATSAT_METRIC_RADIUS`) and
   the pole copied instead of extrapolated (`ATSAT_BC_POLE_COPY`) — both defaults now — a step of
   0.001 and of 0.005 runs 224 iterations clean; 0.025 does not. The default is 0.001: 1.06 s
   per iteration, 224 iterations = 4 minutes. That is still far from anything climatic.
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
