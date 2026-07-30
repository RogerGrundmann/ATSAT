#pragma once

#include "cSaturnModel.h"

#include <vector>
#include <cstdlib>
#include <cmath>
#include <chrono>
#include <cstdio>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

/*
 * Poisson solver for the dynamic pressure, mirrored from ATJUP's PressureSolverJup.h.
 *
 * Until now this header was a fifteen-line shim whose run() forwarded to
 * cSaturnModel::computePressure() in Pressure_Sat.cpp. That function is still there and is
 * still the default; this class is the ATJUP-shaped solver beside it, selected with
 * ATSAT_PRESS_SOLVER=1, so the two can be run against each other.
 *
 * ===== WHAT THE MIRROR CHANGES, AND WHY IT IS NOT COSMETIC =====
 *
 * THE SOURCE TERM. This is the substantive one. A projection method solves
 *
 *     div(grad p) = div(u*)
 *
 * where u* is the velocity tendency WITHOUT the pressure gradient — which is exactly what
 * ATSAT's aux_u/aux_v/aux_w hold: RHS_Sat_Turb.cpp ends with aux_u = rhs_u + dpdr, i.e. the
 * right-hand side with the pressure term added back out. ATJUP uses div(aux) and so do
 * ATOM_turb and ATOM_Precipitation.
 *
 * ATSAT's computePressure() instead forms
 *
 *     (aux_u[i+1] - aux_u[i-1]) - (rhs_u[i+1] - rhs_u[i-1])
 *
 * and aux - rhs IS dpdr, identically, by the line that defined it. So the source is not the
 * divergence of the intermediate velocity at all — it is the divergence of the pressure
 * gradient the previous stage already computed. The equation actually being relaxed is
 *
 *     Laplacian_compact(p) = Laplacian_wide(p)
 *
 * the left side being the 7-point stencil below, the right side the same operator built from
 * centred first derivatives of centred first derivatives, which spans i-2..i+2. Two
 * discretisations of the same operator, set equal to each other. Any smooth field satisfies
 * that to truncation error, so the sweep is a checkerboard filter on p_dyn and nothing more.
 * THE VELOCITY DIVERGENCE NEVER ENTERS. Whatever ATSAT's pressure field is doing, it is not
 * enforcing mass conservation, and that is true no matter how many sweeps it is given.
 *
 * This is inherited, not invented here: old ATOM/atmosphere/Pressure_Atm.cpp has the same
 * subtraction, and ATOM_turb dropped it. ATSAT was branched from the older one.
 * ATSAT_PRESS_SRC=0 restores the legacy source inside this solver, so the source term can be
 * attributed separately from everything else the mirror brings.
 *
 * THE POISSON METRIC IS ALREADY RIGHT HERE, and no knob is provided for it. ATJUP carries
 * ATJUP_POISSON_METRIC because its stencil weights had the metric factors of the DIVERGENCE
 * (1/r and 1/(r sin)) where the Laplacian needs 1/r^2 and 1/(r^2 sin^2). ATSAT's
 * computePressure() writes num2 = 1/(rm^2 dthe^2) and num3 = 1/((rm sin)^2 dphi^2) — the
 * corrected form, already. The bug ATJUP had to find does not exist in this file and the
 * mirror must not import it.
 *
 * THE RADIAL AND POLAR BOUNDARIES of p_dyn move from the 3-point cubic
 * (p[3] - 3p[2] + 3p[1]) to the 2-point Neumann form (4/3)p[1] - (1/3)p[2]. The cubic
 * amplifies an alternating error 7x per call, and p_dyn is not touched by BC_Sat, so a value
 * it puts at i=0 survives untouched through the whole following RK step and drives dpdr at
 * i=1. Same reasoning as the boundaries given to k* and dis*.
 *
 * THE BOUNDARY PASSES MOVE INSIDE THE SWEEP LOOP. They are the outer boundary condition of
 * the relaxation; with more than one sweep, leaving them outside lets the interior run away
 * from its own edges for all but the last pass.
 *
 * ===== A DEFECT CARRIED OVER DELIBERATELY =====
 *
 * ATJUP's relaxation is `#pragma omp parallel for collapse(2)` over (i,j) while writing
 * p_dyn in place and reading p_dyn[i±1][j±1] — threads read cells their neighbours are
 * concurrently writing. That makes the sweep neither Jacobi nor Gauss-Seidel but a hybrid
 * whose result depends on the thread count and the scheduling, i.e. a data race. ATSAT's
 * computePressure() has a SERIAL Poisson loop and does not have this problem.
 *
 * It is mirrored here rather than quietly fixed, because it is ATJUP's and the point of this
 * class is to be comparable to it. But it must not be switched on without knowing: ATSAT is
 * already not reproducible run to run at 24 threads (measured 2026-07-30), and this would add
 * a second source. Run this solver at OMP_NUM_THREADS=1 for anything that has to be compared.
 *
 * ===== WHAT IS NOT MIRRORED =====
 *
 * Everything to do with the obstacle: the land mask, ATJUP_PRESS_WALL, and the one-sided
 * divergence stencils at SeaMount faces. ATSAT has no solid body — BC_seamount is never
 * called and TurbulenceSat::is_land returns false outright — so all of it would be inert. The
 * rigid-lid condition on aux_u is likewise not imported; ATSAT's own polar treatment
 * (aux_v = aux_w = 0 at the poles) is kept instead, since that is this model's condition and
 * the mirror is not the place to change it.
 */
class PressureSolverSat {
public:
    explicit PressureSolverSat(cSaturnModel& model) : m(model) {}

    // Selects between the inherited solver and the mirrored one. Default 0 = the legacy
    // cSaturnModel::computePressure(), so the model is bit-identical until this is set.
    static int mirrored_enabled(){
        static const int v = [](){
            const char* e = getenv("ATSAT_PRESS_SOLVER"); return e ? atoi(e) : 0; }();
        return v;
    }

    void run(){
        if(mirrored_enabled() == 0){ m.computePressure(); return; }
        runMirrored();
    }

private:
    cSaturnModel& m;

    void runMirrored(){
        using namespace std;
        cout << endl << "      ATSAT: PressureSolverSat (mirrored from ATJUP)" << endl;

        auto begin = std::chrono::high_resolution_clock::now();

        // sin(theta) table — depends on j alone. The floor is ATSAT's own rule from
        // computePressure(): exactly zero is replaced by 1e-5. ATJUP clamps with a configurable
        // sinthe_min instead; ATSAT has no such setting and this is not the place to add one.
        std::vector<double> sinthe_table(m.jm);
        for(int j = 0; j < m.jm; j++){
            sinthe_table[j] = sin(m.the.z[j]);
            // Same floor the momentum equations use (ATSAT_SINTHE_MIN, 0 by default), so the
            // two halves of the projection cannot disagree about the polar metric — which is a
            // split ATJUP had to find and fix in its own pair of files.
            if(sinthe_table[j] < cSaturnModel::sinthe_min())
                sinthe_table[j] = cSaturnModel::sinthe_min();
            if(sinthe_table[j] == 0.0) sinthe_table[j] = 1.0e-5;
        }

        // Divergence-source clamp for the p_dyn update (discrete max principle, from ATOM via
        // ATJUP). Bounds each cell's source contribution so a velocity spike cannot drive p_dyn
        // unbounded at the converging-meridian pole. ATSAT's p_dyn runs to ~0.14 bar, so the
        // 0.2 default is a backstop and not a working limit.
        static const double p_dyn_cap = [](){
            const char* e = getenv("ATSAT_PDYN_CAP"); return e ? atof(e) : 0.2; }();

        // Source term: 1 = div(aux), the projection source ATJUP and ATOM_turb use;
        // 0 = div(aux - rhs), the inherited form, which is div(grad p). See the header note.
        static const int src_projection = [](){
            const char* e = getenv("ATSAT_PRESS_SRC"); return e ? atoi(e) : 1; }();

        // One Gauss-Seidel sweep per call is what this model has always done, and a single
        // sweep moves information one cell — the elliptic problem is never solved, p_dyn is a
        // local smoothing of the source. Raising this is the honest way to ask whether the
        // pressure response is merely unconverged. Cost is linear in the count.
        static const int n_sweeps = [](){
            const char* e = getenv("ATSAT_PRESS_SWEEPS");
            const int v = e ? atoi(e) : 1;
            return v > 0 ? v : 1; }();

        // ---- Boundary preparation of the intermediate velocity ----
        // Kept as computePressure() had it, including the 3-point extrapolation, because these
        // are ATSAT's conditions on aux/rhs and the divergence stencil below reads them. rhs_*
        // is prepared too: with ATSAT_PRESS_SRC=0 the legacy source needs it, and nothing else
        // in the model writes rhs_* at the boundaries.
        // Rigid radial walls (ATSAT_BC_RIGID_LID, default 0 = off, as ATSAT has always run).
        // aux_u is the wall-NORMAL intermediate velocity and it feeds du_dr in the divergence
        // source below, so it carries whatever condition the radial walls are meant to impose.
        // Extrapolating it re-injects a wall-normal flux into the projection and leaves the
        // column mass budget open; setting it to zero closes it, which is what a rigid lid and
        // a rigid deep boundary mean. ATJUP defaults this ON; ATSAT does not, because ATSAT's
        // i=0 is the deep interior of a gas giant rather than a floor and whether a lid belongs
        // there at all is a modelling question this port does not get to settle.
        static const int rigid_lid = [](){
            const char* e = getenv("ATSAT_BC_RIGID_LID"); return e ? atoi(e) : 0; }();

        #pragma omp parallel for
        for(int j = 1; j < m.jm-1; j++){          // r-direction
            for(int k = 1; k < m.km-1; k++){
                if(rigid_lid){
                    m.aux_u.x[0][j][k]      = 0.0;
                    m.aux_u.x[m.im-1][j][k] = 0.0;
                } else {
                m.aux_u.x[0][j][k] = m.aux_u.x[3][j][k]
                    - 3.0 * m.aux_u.x[2][j][k] + 3.0 * m.aux_u.x[1][j][k];
                m.aux_u.x[m.im-1][j][k] = m.aux_u.x[m.im-4][j][k]
                    - 3.0 * m.aux_u.x[m.im-3][j][k] + 3.0 * m.aux_u.x[m.im-2][j][k];
                }

                m.aux_v.x[0][j][k] = m.aux_v.x[3][j][k]
                    - 3.0 * m.aux_v.x[2][j][k] + 3.0 * m.aux_v.x[1][j][k];
                m.aux_v.x[m.im-1][j][k] = m.aux_v.x[m.im-4][j][k]
                    - 3.0 * m.aux_v.x[m.im-3][j][k] + 3.0 * m.aux_v.x[m.im-2][j][k];

                m.aux_w.x[0][j][k] = m.aux_w.x[3][j][k]
                    - 3.0 * m.aux_w.x[2][j][k] + 3.0 * m.aux_w.x[1][j][k];
                m.aux_w.x[m.im-1][j][k] = m.aux_w.x[m.im-4][j][k]
                    - 3.0 * m.aux_w.x[m.im-3][j][k] + 3.0 * m.aux_w.x[m.im-2][j][k];

                m.rhs_u.x[0][j][k] = m.rhs_u.x[3][j][k]
                    - 3.0 * m.rhs_u.x[2][j][k] + 3.0 * m.rhs_u.x[1][j][k];
                m.rhs_u.x[m.im-1][j][k] = m.rhs_u.x[m.im-4][j][k]
                    - 3.0 * m.rhs_u.x[m.im-3][j][k] + 3.0 * m.rhs_u.x[m.im-2][j][k];

                m.rhs_v.x[0][j][k] = m.rhs_v.x[3][j][k]
                    - 3.0 * m.rhs_v.x[2][j][k] + 3.0 * m.rhs_v.x[1][j][k];
                m.rhs_v.x[m.im-1][j][k] = m.rhs_v.x[m.im-4][j][k]
                    - 3.0 * m.rhs_v.x[m.im-3][j][k] + 3.0 * m.rhs_v.x[m.im-2][j][k];

                m.rhs_w.x[0][j][k] = m.rhs_w.x[3][j][k]
                    - 3.0 * m.rhs_w.x[2][j][k] + 3.0 * m.rhs_w.x[1][j][k];
                m.rhs_w.x[m.im-1][j][k] = m.rhs_w.x[m.im-4][j][k]
                    - 3.0 * m.rhs_w.x[m.im-3][j][k] + 3.0 * m.rhs_w.x[m.im-2][j][k];
            }
        }

        #pragma omp parallel for
        for(int k = 1; k < m.km-1; k++){          // theta-direction
            for(int i = 1; i < m.im-1; i++){
                m.aux_u.x[i][0][k] = m.aux_u.x[i][3][k]
                    - 3.0 * m.aux_u.x[i][2][k] + 3.0 * m.aux_u.x[i][1][k];
                m.aux_v.x[i][0][k] = 0.0;
                m.aux_w.x[i][0][k] = 0.0;

                m.aux_u.x[i][m.jm-1][k] = m.aux_u.x[i][m.jm-4][k]
                    - 3.0 * m.aux_u.x[i][m.jm-3][k] + 3.0 * m.aux_u.x[i][m.jm-2][k];
                m.aux_v.x[i][m.jm-1][k] = 0.0;
                m.aux_w.x[i][m.jm-1][k] = 0.0;

                m.rhs_u.x[i][0][k] = m.rhs_u.x[i][3][k]
                    - 3.0 * m.rhs_u.x[i][2][k] + 3.0 * m.rhs_u.x[i][1][k];
                m.rhs_v.x[i][0][k] = 0.0;
                m.rhs_w.x[i][0][k] = 0.0;

                m.rhs_u.x[i][m.jm-1][k] = m.rhs_u.x[i][m.jm-4][k]
                    - 3.0 * m.rhs_u.x[i][m.jm-3][k] + 3.0 * m.rhs_u.x[i][m.jm-2][k];
                m.rhs_v.x[i][m.jm-1][k] = 0.0;
                m.rhs_w.x[i][m.jm-1][k] = 0.0;
            }
        }

        #pragma omp parallel for
        for(int i = 0; i < m.im; i++){            // phi-direction
            for(int j = 0; j < m.jm; j++){
                m.aux_u.x[i][j][0] = m.c43 * m.aux_u.x[i][j][1] - m.c13 * m.aux_u.x[i][j][2];
                m.aux_u.x[i][j][m.km-1] = m.c43 * m.aux_u.x[i][j][m.km-2] - m.c13 * m.aux_u.x[i][j][m.km-3];
                m.aux_u.x[i][j][0] = m.aux_u.x[i][j][m.km-1] =
                    (m.aux_u.x[i][j][0] + m.aux_u.x[i][j][m.km-1])/2.0;

                m.aux_v.x[i][j][0] = m.c43 * m.aux_v.x[i][j][1] - m.c13 * m.aux_v.x[i][j][2];
                m.aux_v.x[i][j][m.km-1] = m.c43 * m.aux_v.x[i][j][m.km-2] - m.c13 * m.aux_v.x[i][j][m.km-3];
                m.aux_v.x[i][j][0] = m.aux_v.x[i][j][m.km-1] =
                    (m.aux_v.x[i][j][0] + m.aux_v.x[i][j][m.km-1])/2.0;

                m.aux_w.x[i][j][0] = m.c43 * m.aux_w.x[i][j][1] - m.c13 * m.aux_w.x[i][j][2];
                m.aux_w.x[i][j][m.km-1] = m.c43 * m.aux_w.x[i][j][m.km-2] - m.c13 * m.aux_w.x[i][j][m.km-3];
                m.aux_w.x[i][j][0] = m.aux_w.x[i][j][m.km-1] =
                    (m.aux_w.x[i][j][0] + m.aux_w.x[i][j][m.km-1])/2.0;

                m.rhs_u.x[i][j][0] = m.c43 * m.rhs_u.x[i][j][1] - m.c13 * m.rhs_u.x[i][j][2];
                m.rhs_u.x[i][j][m.km-1] = m.c43 * m.rhs_u.x[i][j][m.km-2] - m.c13 * m.rhs_u.x[i][j][m.km-3];
                m.rhs_u.x[i][j][0] = m.rhs_u.x[i][j][m.km-1] =
                    (m.rhs_u.x[i][j][0] + m.rhs_u.x[i][j][m.km-1])/2.0;

                m.rhs_v.x[i][j][0] = m.c43 * m.rhs_v.x[i][j][1] - m.c13 * m.rhs_v.x[i][j][2];
                m.rhs_v.x[i][j][m.km-1] = m.c43 * m.rhs_v.x[i][j][m.km-2] - m.c13 * m.rhs_v.x[i][j][m.km-3];
                m.rhs_v.x[i][j][0] = m.rhs_v.x[i][j][m.km-1] =
                    (m.rhs_v.x[i][j][0] + m.rhs_v.x[i][j][m.km-1])/2.0;

                m.rhs_w.x[i][j][0] = m.c43 * m.rhs_w.x[i][j][1] - m.c13 * m.rhs_w.x[i][j][2];
                m.rhs_w.x[i][j][m.km-1] = m.c43 * m.rhs_w.x[i][j][m.km-2] - m.c13 * m.rhs_w.x[i][j][m.km-3];
                m.rhs_w.x[i][j][0] = m.rhs_w.x[i][j][m.km-1] =
                    (m.rhs_w.x[i][j][0] + m.rhs_w.x[i][j][m.km-1])/2.0;
            }
        }

        // Grid-spacing reciprocals — constant for the entire grid.
        const double inv_2dr   = 1.0 / (2.0 * m.dr);
        const double inv_2dthe = 1.0 / (2.0 * m.dthe);
        const double inv_2dphi = 1.0 / (2.0 * m.dphi);
        const double inv_dr2   = 1.0 / (m.dr   * m.dr);
        const double inv_dthe2 = 1.0 / (m.dthe * m.dthe);
        const double inv_dphi2 = 1.0 / (m.dphi * m.dphi);

        for(int sweep = 0; sweep < n_sweeps; sweep++){

            // See the header note on the race this pragma carries over from ATJUP.
            #pragma omp parallel for collapse(2) schedule(dynamic, 4)
            for(int i = 1; i < m.im-1; i++){
                for(int j = 1; j < m.jm-1; j++){

                    // Geometry once per (i,j), in the same CellGeometry RHSSat is given.
                    cSaturnModel::CellGeometry geo;
                    geo.rm       = m.metricRadius(m.rad.z[i]);
                    geo.rm2      = geo.rm * geo.rm;
                    geo.exp_rm   = m.coord_stretching ? 1.0 / (geo.rm + 1.0) : 1.0;
                    geo.exp_2_rm = geo.exp_rm * geo.exp_rm;
                    geo.sinthe   = sinthe_table[j];
                    geo.sinthe2  = geo.sinthe * geo.sinthe;
                    geo.costhe   = cos(m.the.z[j]);
                    geo.cotanthe             = geo.costhe / geo.sinthe;
                    geo.inv_rm               = 1.0 / geo.rm;
                    geo.inv_rm2              = 1.0 / geo.rm2;
                    geo.inv_rmsinthe         = 1.0 / (geo.rm * geo.sinthe);
                    geo.inv_rm2sinthe        = geo.inv_rm2 / geo.sinthe;
                    geo.inv_rm2sinthe2       = geo.inv_rm2 / geo.sinthe2;
                    geo.costhe_inv_rm2sinthe = geo.costhe * geo.inv_rm2sinthe;
                    geo.inv_2dr   = inv_2dr;
                    geo.inv_2dthe = inv_2dthe;
                    geo.inv_2dphi = inv_2dphi;
                    geo.inv_dr2   = inv_dr2;
                    geo.inv_dthe2 = inv_dthe2;
                    geo.inv_dphi2 = inv_dphi2;

                    // Stencil weights. These are the Laplacian's own metric factors — 1/r^2 on
                    // theta and 1/(r^2 sin^2) on phi — which is what computePressure() already
                    // used. ATJUP's ATJUP_POISSON_METRIC exists because its weights were the
                    // divergence's; there is nothing to correct here.
                    const double num1 = geo.exp_2_rm     * inv_dr2;
                    const double num2 = geo.inv_rm2      * inv_dthe2;
                    const double num3 = geo.inv_rm2sinthe2 * inv_dphi2;
                    const double denom = 2.0 * num1 + 2.0 * num2 + 2.0 * num3;

                    for(int k = 1; k < m.km-1; k++){

                        // ---- Divergence source ----
                        double du_dr, dv_dthe, dw_dphi;
                        if(src_projection != 0){
                            du_dr   = (m.aux_u.x[i+1][j][k] - m.aux_u.x[i-1][j][k]) * inv_2dr;
                            dv_dthe = (m.aux_v.x[i][j+1][k] - m.aux_v.x[i][j-1][k]) * inv_2dthe;
                            dw_dphi = (m.aux_w.x[i][j][k+1] - m.aux_w.x[i][j][k-1]) * inv_2dphi;
                        } else {
                            du_dr   = ((m.aux_u.x[i+1][j][k] - m.aux_u.x[i-1][j][k])
                                     - (m.rhs_u.x[i+1][j][k] - m.rhs_u.x[i-1][j][k])) * inv_2dr;
                            dv_dthe = ((m.aux_v.x[i][j+1][k] - m.aux_v.x[i][j-1][k])
                                     - (m.rhs_v.x[i][j+1][k] - m.rhs_v.x[i][j-1][k])) * inv_2dthe;
                            dw_dphi = ((m.aux_w.x[i][j][k+1] - m.aux_w.x[i][j][k-1])
                                     - (m.rhs_w.x[i][j][k+1] - m.rhs_w.x[i][j][k-1])) * inv_2dphi;
                        }

                        double div_src = du_dr   * geo.exp_rm
                                       + dv_dthe * geo.inv_rm
                                       + dw_dphi * geo.inv_rmsinthe;

                        const double acc =
                              (m.p_dyn.x[i+1][j][k] + m.p_dyn.x[i-1][j][k]) * num1
                            + (m.p_dyn.x[i][j+1][k] + m.p_dyn.x[i][j-1][k]) * num2
                            + (m.p_dyn.x[i][j][k+1] + m.p_dyn.x[i][j][k-1]) * num3;

                        // Discrete max principle on the source, plus a guard against a
                        // non-finite one, so p_dyn cannot run away to NaN at the pole.
                        const double src_max = denom * p_dyn_cap;
                        if      (!std::isfinite(div_src)) div_src = 0.0;
                        else if (div_src >  src_max)      div_src =  src_max;
                        else if (div_src < -src_max)      div_src = -src_max;

                        m.p_dyn.x[i][j][k] = (acc - div_src) / denom;
                    } // k
                } // j
            } // i

            // ---- Outer boundary condition of the relaxation ----
            // 2-point Neumann, not the 3-point cubic computePressure() used: p_dyn is absent
            // from BC_Sat, so whatever lands at i=0 here survives the whole following RK step.
            #pragma omp parallel for collapse(2)
            for(int k = 0; k < m.km; k++){
                for(int j = 0; j < m.jm; j++){
                    m.p_dyn.x[0][j][k] = m.c43 * m.p_dyn.x[1][j][k] - m.c13 * m.p_dyn.x[2][j][k];
                    m.p_dyn.x[m.im-1][j][k] = m.c43 * m.p_dyn.x[m.im-2][j][k]
                                            - m.c13 * m.p_dyn.x[m.im-3][j][k];
                }
            }

            #pragma omp parallel for collapse(2)
            for(int k = 0; k < m.km; k++){
                for(int i = 0; i < m.im; i++){
                    m.p_dyn.x[i][0][k] = m.c43 * m.p_dyn.x[i][1][k] - m.c13 * m.p_dyn.x[i][2][k];
                    m.p_dyn.x[i][m.jm-1][k] = m.c43 * m.p_dyn.x[i][m.jm-2][k]
                                            - m.c13 * m.p_dyn.x[i][m.jm-3][k];
                }
            }

            #pragma omp parallel for collapse(2)
            for(int i = 0; i < m.im; i++){
                for(int j = 0; j < m.jm; j++){
                    m.p_dyn.x[i][j][0] = m.c43 * m.p_dyn.x[i][j][1] - m.c13 * m.p_dyn.x[i][j][2];
                    m.p_dyn.x[i][j][m.km-1] = m.c43 * m.p_dyn.x[i][j][m.km-2]
                                            - m.c13 * m.p_dyn.x[i][j][m.km-3];
                    m.p_dyn.x[i][j][0] = m.p_dyn.x[i][j][m.km-1] =
                        (m.p_dyn.x[i][j][0] + m.p_dyn.x[i][j][m.km-1])/2.0;
                }
            }

        } // sweep

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for PressureSolverSat\n", elapsed.count() * 1e-9);
        cout << "      ATSAT: PressureSolverSat ended" << endl;
    }
};
