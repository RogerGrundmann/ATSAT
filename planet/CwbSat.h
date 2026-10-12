/*
 * Saturn Atmosphere Circulation Model (ATSAT)
 * ATSAT_CWB_DIAG=<n> -- the column budget of H2O, NH3 and CH4 (vapour + cloud + ice each),
 * PRINT-ONLY, default off; printed every n iterations. ATJUP's ATJUP_CWB_DIAG (planet/CwbJup.h),
 * ported 2026-10-12 and widened from H2O to the three condensables.
 *
 * WHY. At the defaults of 2a67d19 the columns fall between the saturation adjustment's calls:
 * H2O -0.454 % and NH3 -0.460 % in 224 iterations, accelerating; CH4 -0.151 % in the first step.
 * Arms with the precipitation off, the species wall and lid copy off, and the pole fill off all
 * show the same loss (satchk/giants/PS006, CD0, CDp, CDw, CDf), so it is in the time step, and
 * nothing in the model said in which term.
 *
 * This instrument marks each column after every stage of the time loop that can write the nine
 * fields and charges the difference to that stage. The Runge-Kutta step is split further:
 *   rk_reset   it integrates from the n-level copies, which restoreVar refreshed at the END of
 *              the previous iteration -- every direct write made since then and not copied into
 *              them (ATSAT_RK_SCALAR_SYNC) is discarded here;
 *   rk terms   the RK4-weighted column integrals of each right-hand-side term of the three
 *              equations of a gas, accumulated per thread inside RHSSat: radial / meridional /
 *              zonal transport (advective form, -u.grad q), diffusion, the chemistry source
 *              (NH3 only), the precipitation source;
 *   rk_rest    what the step changed beyond reset + terms: the floor of the final update, the
 *              limiters and round-off. An identity check when none of them binds.
 * q*div(u) is accumulated beside them and NOT applied: the transport is written in advective
 * form, and -u.grad q = -div(q u) + q div(u), so on a closed column the advective form creates or
 * destroys exactly the integral of q*div(u). Printed so the two can be compared.
 *
 * Unit: g/m2, the mean over the sphere weighted by sin(colatitude), layer thickness in metres on
 * densities in kg/m3 -- the same sum the saturation adjustment's instrument prints. The top level
 * has no thickness and is not counted. The plane i = 0 is not integrated (bcRadius sets it) and
 * its part of the change is printed apart from the integrated cells.
 *
 * The state lives in a function-local static so that the three translation units share one
 * instance without a model member (C++11 has no inline variables).
 */

#pragma once

namespace SatCwb {

    enum { G_H2O = 0, G_NH3, G_CH4, NGAS };
    enum { T_RAD = 0, T_THE, T_PHI, T_DIFF, T_CHEM, T_PRECIP, T_QDIV, NTERM };
    const int MAXTHR = 256;

    struct State {
        bool   on;                               // ATSAT_CWB_DIAG
        double stage_wgt;                        // RK4 weight of the stage being evaluated (1, 2, 2, 1)
        double acc[MAXTHR][NGAS][8];             // [thread][gas][term], weighted sums of the tendencies
    };

    inline State& st(){
        static State s = State();
        return s;
    }
}
