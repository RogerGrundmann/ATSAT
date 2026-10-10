/*
 * Saturn Atmosphere General Circulation Modell (ATSAT)
 * Saturation adjustment and mixed-phase partitioning, mirrored from ATJUP's
 * SaturationAdjustmentJup.h / .cpp.
 *
 * Reference algorithm:
 *   Tao, W.-K., Simpson, J., and McCumber, M.:
 *   "An Ice-Water Saturation Adjustment", AMS Notes and Correspondence, 1988.
 *
 * ===== WHY THIS FILE EXISTS =====
 *
 * This header was a declaration-only shim whose run() forwarded to
 * cSaturnModel::Saturation_Adjustment() in Weather_Sat.cpp. That function is still there;
 * until 2026-10-10 it was the default. THE DEFAULT IS NOW THE SHARED ALGORITHM (the user's
 * decision); ATSAT_SATADJ=0 selects the inherited routine, so the two can still be measured
 * against each other.
 *
 * ATSAT's version and ATJUP's had drifted apart in ten places. In descending order of
 * consequence:
 *
 * 1. THE SATURATION VAPOUR PRESSURE WAS FROZEN AT THE CELL'S ENTRY TEMPERATURE.
 *    Inside the mixed-phase iteration, Weather_Sat.cpp:177-180 recomputes E_Rain and E_Ice
 *    from t_u — the temperature the cell had on entry — while the iteration's whole purpose
 *    is that condensation releases latent heat, which raises T, which raises the saturation
 *    vapour pressure, which limits how much more can condense. ATJUP passes T. With t_u the
 *    saturation targets q_Rain and q_Ice are CONSTANT across all thirty passes, so the
 *    feedback the Tao scheme is built around is absent and the loop merely relaxes q_v_hyp
 *    toward a fixed number. The two lines immediately above the loop (Weather_Sat.cpp:122-125)
 *    legitimately use t_u, because there T == t_u; this looks like those lines copied
 *    downward. The commented-out Clausius-Clapeyron lines beside them use T.
 *
 * 2. NO GUARD ON UNPHYSICAL CELLS. ATJUP skips any cell with t_u <= 0 or p_u <= 0, written
 *    as !(x > 0.0) so an already-arrived NaN is skipped rather than propagated, and its
 *    comment records this as one of two seeds of a domain-wide NaN. ATSAT has no such test,
 *    and both formulas are singular there: saturation_vapour_pressure forms -L0/T (-inf at
 *    T=0) and del_alf*log(T) (0*-inf = NaN), and q_Rain divides by p_u. The existing
 *    `c >= q_Rain` entry test cannot catch it, because every comparison against NaN is false.
 *
 * 3. THE DENSITY IS INCONSISTENT INSIDE ATSAT — but the fix is NOT simply "use the local one".
 *    ATJUP's corresponding line reads m.rho_at(i,j,k), and ATJUP's rho_at() begins
 *    `if(!local_rho()) return r_mix;` — its default is the REFERENCE density, with the local
 *    one opt-in behind ATJUP_LOCAL_RHO, and its saturation adjustment and its precipitation
 *    both go through that one gate so they always agree. ATSAT's rho_at() has no gate and
 *    always returns the local value. So PrecipitationSat.h:212 saturates on the local density
 *    while Weather_Sat.cpp saturates on r_mix, and PrecipitationSat's own comment — "Same
 *    density choice as SaturationAdjustmentSat — the two must agree or the microphysics and
 *    the adjustment would saturate at different values" — is not true of the code as it
 *    stands. This routine reproduces ATJUP's gate locally (ATSAT_LOCAL_RHO, default off =
 *    r_mix); giving rho_at() the gate so that all its callers move together is the real fix
 *    and is left as its own change, because PrecipitationSat and the radiative coupling read
 *    it too. Copying ATJUP's line without its gate diverges the run in eight iterations —
 *    measured, see the note at the call site.
 *
 * 4. NO ICE-PHASE COEFFICIENTS. E_Ice is computed from the LIQUID pair, so deposition and
 *    sublimation see the saturation vapour pressure over liquid where they should see the one
 *    over ice — the vapour-pressure difference that drives the Bergeron process is
 *    understated and ice grows at the expense of droplets too slowly. Recorded when the
 *    precipitation scheme was ported (35077b5) and still true. run_mirrored() below takes the
 *    ice quadruple as ATJUP's does, so the structure is complete. SINCE 2026-10-09 the call
 *    site passes real ice values: the coefficients are properties of the substance, so they are
 *    ATJUP's pairs for H2O and NH3 (cSaturnModel.h has the continuity check); CH4 has one curve.
 *    Before that it passed the liquid values and this point described the model.
 *
 * 5. q_v_b IS NEVER CLAMPED. ATSAT clamps q_c_b and q_i_b to >= 0 inside the loop but not the
 *    vapour itself; ATJUP clamps all three.
 *
 * 6. THE CONVERGENCE TEST DIVIDES BY A POSSIBLY-ZERO q_v_hyp: fabs(q_v_b/q_v_hyp - 1.0),
 *    where ATJUP writes fabs(q_v_b - q_v_hyp)/(q_v_hyp + 1e-20).
 *
 * 7. THE DIAGNOSTIC CAPTURE CONDITION IS MEANINGLESS. Weather_Sat.cpp:204 reads
 *    `if((j == 90) && (k == 180) && (i == iter_prec_found))` — a grid index compared against
 *    an iteration counter, so which cell gets reported depends on how many passes it took.
 *    ATJUP reports the last converged cell.
 *
 * 8. THE SHARED DIAGNOSTICS ARE WRITTEN FROM INSIDE THE PARALLEL REGION WITHOUT SYNCHRONISATION
 *    (sat_found, i_sat, saturation, ...). ATJUP writes them under omp critical. This is one
 *    concrete, locatable contributor to ATSAT not reproducing itself run to run.
 *
 * 9. ITERATION BUDGET AND TOLERANCE: ATSAT 30 passes to 1e-4, ATJUP 15 to 1e-3.
 *
 * 10. THE ENTRY TEST does not reject a saturation vapour pressure that has underflowed to
 *     zero in a very cold cell; ATJUP's `q_Rain_0 <= 0.0` does.
 *
 * ===== ONE DELIBERATE NON-MIRROR =====
 *
 * ATSAT's routine ends by rewriting p_stat from the adjusted temperature
 * (Weather_Sat.cpp:304). ATJUP does not touch p_stat here at all. That is kept, because
 * removing it is a separate decision about where the hydrostatic pressure is allowed to
 * respond to latent heating, and bundling it into this change would make the comparison
 * unreadable. It is flagged rather than settled.
 */

#pragma once

#include "cSaturnModel.h"
#include "SaturationAdjustment.h"   // the SHARED algorithm this class dispatches to

#include <cmath>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

class Array;

class SaturationAdjustmentSat {
public:
    explicit SaturationAdjustmentSat(cSaturnModel& model) : m(model) {}

    // Selects between the inherited routine and the mirrored one. DEFAULT 1 = the shared
    // SaturationAdjustment.h SINCE 2026-10-10 (the user's decision); ATSAT_SATADJ=0 restores the
    // legacy cSaturnModel::Saturation_Adjustment(). Every ATSAT number before that commit was
    // made with the legacy routine unless its run set the knob.
    //
    // THE RUN THAT WAS OWED, 224 iterations, 8 threads (ATJUP/satchk/giants/SAT224i, SAT224s),
    // legacy -> shared: 45 of 100 printed extrema identical, temperature, velocity and pressure
    // among them; max h2o_cloud 79.45 -> 74.10 g/m3 (-6.7 %), max nh3_cloud 1.665 -> 1.506,
    // max nh3_ice 2.80 -> 3.54 (+26 %), max ch4_ice 48.18 -> 37.29 (-23 %), max latent heat
    // 1.472 -> 1.549 W/m3. The shared routine's own budget over its 113 calls: no ice deleted,
    // no cell unconverged, column change of each gas below 1e-6 g/m2. The legacy routine has
    // no such instrument, so what it does to the column is not known. Which of the two cloud
    // amounts is right is still not settled by a measurement; the older note follows.
    //
    // WHAT IT DOES WHEN SET, measured 2026-07-31 over 50 iterations (config_m50, 12 threads),
    // against the legacy routine with everything else identical:
    //
    //     max h2o_cloud   113.254 -> 94.812 g/m3   (-16.3 %), and the peak moves 175 -> 187 km
    //     max latent heat   1.1115 -> 1.0179 W/m2
    //     max/min t, max nh4sh: unchanged to the last digit printed
    //
    // So it condenses less water and puts the deck a layer higher. Which of the two is right is
    // NOT settled by this measurement, and the choice is a real trade rather than a fix:
    //   FOR turning it on  — the mirror carries the repair from the port commit, where the
    //     legacy routine's saturation target was found to be frozen, and it writes its shared
    //     diagnostic state under omp critical instead of from every thread at once (point 8).
    //   AGAINST, until 2026-10-09 — its ice-phase quadruple was the LIQUID one (point 4 below).
    //     The ice coefficients exist now, so that objection is gone, and the run of useful
    //     length is the one at the top of this note.
    static int mirrored_enabled(){
        static const int v = [](){
            const char* e = getenv("ATSAT_SATADJ"); return e ? atoi(e) : 1; }();
        return v;
    }

    // The call-site signature is unchanged, so cSaturnModel.cpp needs no edit. The ice-phase
    // coefficients ATJUP carries are supplied by run() to run_mirrored() below; see point 4 in
    // the note above for why they are currently the liquid ones.
    void run(std::string gas,
             double &coeff_A,   double &coeff_B,
             double &coeff_A_i, double &coeff_B_i,
             double &t_0,       double &t_00,
             double &ep,        double &lv,  double &ls,
             double &cp,        double &r,
             double &C,         double &L0,  double &R,
             double &del_alf,   double &del_bet,  double &m_mol,
             Array &c, Array &cloud, Array &ice);

    // The three static helpers that used to sit here are gone. saturation_vapour_pressure and
    // clausius_clapeyron are ATPhys:: functions now — one implementation for every planet instead
    // of a copy per class — and humility_critical had no caller in either model, which is why it
    // was deleted rather than moved: a shared header is the wrong place to preserve dead code.
    // (cSaturnModel::Humility_critical, capital H, is a DIFFERENT function and is still live.)

private:
    cSaturnModel& m;


};
