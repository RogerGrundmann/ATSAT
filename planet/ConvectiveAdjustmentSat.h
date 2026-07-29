/*
 * Saturn Atmosphere Circulation Model (ATSAT)
 * Dry convective adjustment.
 *
 * Reference algorithm:
 *   Manabe, S. and Strickler, R. F.: "Thermal Equilibrium of the Atmosphere with a Convective
 *   Adjustment", J. Atmos. Sci. 21, 361-385, 1964.
 *
 * PROVENANCE. This is a port of ConvectiveAdjustmentJup from ATJUP, structurally unchanged. It is
 * the first module carried across between the two models, and it was chosen for that precisely
 * because it is self-contained: it reads t, p_stat and SeaMount, writes t, and touches nothing
 * else. Whatever it does here is therefore attributable to the scheme and not to a coupling.
 *
 * WHAT IS AND IS NOT MEASURED. In ATJUP the need for it was established from that model's own
 * restart files: a superadiabatic layer appeared and deepened with the run, and nothing in the
 * model restored a column to its adiabat. **No equivalent measurement exists for ATSAT yet.** Do
 * not read the Jupiter result as a Saturn result — the two models differ in the numbers that
 * decide the criterion. A layer here is L_atm/(im-1) = 500/40 = 12.5 km against Jupiter's 3.5 km,
 * and Saturn's g is about 10.4 m/s2 against Jupiter's 25.9, so the per-layer adiabatic drop works
 * out near 12 K rather than 8 K, from a lapse rate roughly two and a half times gentler. Whether
 * ATSAT develops a superadiabatic layer at all is the first thing this class is here to answer,
 * and until a run says so it is switched off (ATSAT_CONV_ADJ, default 0, bit-identical).
 *
 * WHAT IT DOES. Each column is swept from the bottom up. Wherever a layer pair is steeper than the
 * dry adiabat, the whole unstable segment is mixed to exactly the adiabatic lapse rate, conserving
 * the mass-weighted enthalpy of the segment. Sweeps repeat until the column is stable. Its two
 * properties are the ones that matter: it removes the instability completely rather than damping
 * it, and it does not create or destroy energy — the run reports the residual drift so that claim
 * stays checkable rather than asserted.
 *
 * WHAT IT DELIBERATELY DOES NOT DO. It mixes temperature only, not composition. Real convection
 * carries the species with it, and a moist adjustment would use the saturated adiabat where cloud
 * is present rather than the dry one. Both are reasonable extensions; neither is done here, because
 * each is a modelling decision with consequences for the microphysics in SaturationAdjustmentSat.
 */

#pragma once

#include <cmath>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

class cSaturnModel;

class ConvectiveAdjustmentSat {
public:
    explicit ConvectiveAdjustmentSat(cSaturnModel& model) : m(model) {}

    void run();

private:
    cSaturnModel& m;
};
