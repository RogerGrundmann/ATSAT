/*
 * Atmosphere General Circulation Modell (AGCM) applied to laminar flow
 * Program for the computation of geo-atmospherical circulating flows in a spherical shell
 * Finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 2 additional transport equations to describe the water vapour and co2 concentration
 * 4. order Runge-Kutta scheme to solve 2. order differential equations
*/
#include "cSaturnModel.h"
#include "PressureSolverSat.h"

using namespace std;

// The Poisson solve itself lives in the SHARED planet/PressureSolver.h; PressureSolverSat is a
// typedef of it. What remains here of the pressure code is prepareProjectionBoundaries() below,
// this model's own boundary conditions on the aux and rhs fields, which the shared solver calls.
//
// cSaturnModel::computePressure() used to sit here as a second, inherited solver behind
// ATSAT_PRESS_SOLVER. It is deleted. Its one substantive difference is recorded in
// PressureSolver.h and survives as the <TAG>_PRESS_SRC knob: its divergence source was
// div(aux - rhs), which is div(grad p) and never sees the velocity divergence at all.

/*
 * ONE place where the mixture density is formed, mirroring ATJUP's computeMixtureDensity().
 *
 * Before this, the ideal-gas relation was written out twice — in cSaturnModel::rho_at() and
 * again as a lambda inside ThermalWindDiagSat.h — with different fallbacks for an unusable
 * cell (r_mix in one, 0.0 in the other). Two copies of one formula can drift; now there is
 * one, and rho_at() is a gated accessor over the array it fills.
 *
 * R_mix is in J/(kg K) — ChemistrySat assembles it as a mass-weighted mean of the per-species
 * R and prints it as 2909.2, which is R_universal/M_mix = 8314/2.8596 = 2907.4 to within 0.06%.
 * p_stat is in bars, hence the 1e5. There is NO further factor: this line used to carry
 * R_mix*1.0e3 on the strength of a comment claiming J/(g K), which made every density it
 * produced 1000x too small — 0.0047 kg/m3 at 55 bar and 400 K, against a reference r_mix of
 * 1.382. The buoyancy term in RHS_Sat_Turb.cpp has always divided by R_mix WITHOUT the 1e3,
 * and its diagnostic reproduces the printed 47.7 N/m3, which is what settles which of the two
 * readings was right.
 * An unusable cell stores 0.0 and every reader treats that as "no density", which is ATJUP's
 * convention: rho_at() substitutes r_mix, the direct readers skip the cell.
 *
 * DIFFERENT FROM ATJUP, DELIBERATELY AND FOR NOW: ATJUP builds this from
 * (p_stat + p_dyn*p_dyn_to_bar()), i.e. it includes the dynamic pressure. ATSAT has no
 * p_dyn_to_bar() — that is ATJUP_PDYN_UNITS, still on the list of things ATSAT lacks — and
 * ATJUP's own comment records that adding p_dyn UNCONVERTED counted 7.79x too heavily. So the
 * hydrostatic pressure alone is used here, which is exactly what ATSAT computed before, and
 * this change stays a change of structure rather than of physics. Settling p_dyn's units is
 * the prerequisite, and it is a separate piece of work.
 */
void cSaturnModel::computeMixtureDensity(){
    #pragma omp parallel for collapse(2) schedule(static)
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                const double T = t.x[i][j][k] * t_ref;
                // Written !(T > 0.0) so a NaN lands here instead of propagating.
                if(!(T > 0.0) || !(R_mix > 0.0)){ rho_mix.x[i][j][k] = 0.0; continue; }
                const double rho = (p_stat.x[i][j][k] * 1.0e5) / (R_mix * T);
                rho_mix.x[i][j][k] = std::isfinite(rho) ? rho : 0.0;
            }
        }
    }
}

/*
 * The boundary conditions this model imposes on aux_* and rhs_* before the projection takes their
 * divergence, called by the SHARED PressureSolver.h. They are computePressure()'s own conditions,
 * including the 3-point cubic extrapolation, and they stay with the model rather than moving into
 * the shared solver because they ARE this model's boundary conditions — ATJUP's are different (one
 * radial pass of 2-point Neumann on aux only) for reasons that belong to ATJUP. rhs_* is prepared
 * too: the legacy source (ATSAT_PRESS_SRC=0) needs it, and nothing else in the model writes rhs_*
 * at the boundaries.
 *
 * rigid_lid comes from the solver, which owns the knob. aux_u is the wall-NORMAL intermediate
 * velocity and feeds du_dr in the divergence source, so it carries whatever condition the radial
 * walls impose; extrapolating it re-injects a wall-normal flux and leaves the column mass budget
 * open, zeroing it closes it.
 */
void cSaturnModel::prepareProjectionBoundaries(bool rigid_lid){
        #pragma omp parallel for
        for(int j = 1; j < jm-1; j++){          // r-direction
            for(int k = 1; k < km-1; k++){
                if(rigid_lid){
                    aux_u.x[0][j][k]      = 0.0;
                    aux_u.x[im-1][j][k] = 0.0;
                } else {
                aux_u.x[0][j][k] = aux_u.x[3][j][k]
                    - 3.0 * aux_u.x[2][j][k] + 3.0 * aux_u.x[1][j][k];
                aux_u.x[im-1][j][k] = aux_u.x[im-4][j][k]
                    - 3.0 * aux_u.x[im-3][j][k] + 3.0 * aux_u.x[im-2][j][k];
                }

                aux_v.x[0][j][k] = aux_v.x[3][j][k]
                    - 3.0 * aux_v.x[2][j][k] + 3.0 * aux_v.x[1][j][k];
                aux_v.x[im-1][j][k] = aux_v.x[im-4][j][k]
                    - 3.0 * aux_v.x[im-3][j][k] + 3.0 * aux_v.x[im-2][j][k];

                aux_w.x[0][j][k] = aux_w.x[3][j][k]
                    - 3.0 * aux_w.x[2][j][k] + 3.0 * aux_w.x[1][j][k];
                aux_w.x[im-1][j][k] = aux_w.x[im-4][j][k]
                    - 3.0 * aux_w.x[im-3][j][k] + 3.0 * aux_w.x[im-2][j][k];

                rhs_u.x[0][j][k] = rhs_u.x[3][j][k]
                    - 3.0 * rhs_u.x[2][j][k] + 3.0 * rhs_u.x[1][j][k];
                rhs_u.x[im-1][j][k] = rhs_u.x[im-4][j][k]
                    - 3.0 * rhs_u.x[im-3][j][k] + 3.0 * rhs_u.x[im-2][j][k];

                rhs_v.x[0][j][k] = rhs_v.x[3][j][k]
                    - 3.0 * rhs_v.x[2][j][k] + 3.0 * rhs_v.x[1][j][k];
                rhs_v.x[im-1][j][k] = rhs_v.x[im-4][j][k]
                    - 3.0 * rhs_v.x[im-3][j][k] + 3.0 * rhs_v.x[im-2][j][k];

                rhs_w.x[0][j][k] = rhs_w.x[3][j][k]
                    - 3.0 * rhs_w.x[2][j][k] + 3.0 * rhs_w.x[1][j][k];
                rhs_w.x[im-1][j][k] = rhs_w.x[im-4][j][k]
                    - 3.0 * rhs_w.x[im-3][j][k] + 3.0 * rhs_w.x[im-2][j][k];
            }
        }

        #pragma omp parallel for
        for(int k = 1; k < km-1; k++){          // theta-direction
            for(int i = 1; i < im-1; i++){
                aux_u.x[i][0][k] = aux_u.x[i][3][k]
                    - 3.0 * aux_u.x[i][2][k] + 3.0 * aux_u.x[i][1][k];
                aux_v.x[i][0][k] = 0.0;
                aux_w.x[i][0][k] = 0.0;

                aux_u.x[i][jm-1][k] = aux_u.x[i][jm-4][k]
                    - 3.0 * aux_u.x[i][jm-3][k] + 3.0 * aux_u.x[i][jm-2][k];
                aux_v.x[i][jm-1][k] = 0.0;
                aux_w.x[i][jm-1][k] = 0.0;

                rhs_u.x[i][0][k] = rhs_u.x[i][3][k]
                    - 3.0 * rhs_u.x[i][2][k] + 3.0 * rhs_u.x[i][1][k];
                rhs_v.x[i][0][k] = 0.0;
                rhs_w.x[i][0][k] = 0.0;

                rhs_u.x[i][jm-1][k] = rhs_u.x[i][jm-4][k]
                    - 3.0 * rhs_u.x[i][jm-3][k] + 3.0 * rhs_u.x[i][jm-2][k];
                rhs_v.x[i][jm-1][k] = 0.0;
                rhs_w.x[i][jm-1][k] = 0.0;
            }
        }

        #pragma omp parallel for
        for(int i = 0; i < im; i++){            // phi-direction
            for(int j = 0; j < jm; j++){
                aux_u.x[i][j][0] = c43 * aux_u.x[i][j][1] - c13 * aux_u.x[i][j][2];
                aux_u.x[i][j][km-1] = c43 * aux_u.x[i][j][km-2] - c13 * aux_u.x[i][j][km-3];
                aux_u.x[i][j][0] = aux_u.x[i][j][km-1] =
                    (aux_u.x[i][j][0] + aux_u.x[i][j][km-1])/2.0;

                aux_v.x[i][j][0] = c43 * aux_v.x[i][j][1] - c13 * aux_v.x[i][j][2];
                aux_v.x[i][j][km-1] = c43 * aux_v.x[i][j][km-2] - c13 * aux_v.x[i][j][km-3];
                aux_v.x[i][j][0] = aux_v.x[i][j][km-1] =
                    (aux_v.x[i][j][0] + aux_v.x[i][j][km-1])/2.0;

                aux_w.x[i][j][0] = c43 * aux_w.x[i][j][1] - c13 * aux_w.x[i][j][2];
                aux_w.x[i][j][km-1] = c43 * aux_w.x[i][j][km-2] - c13 * aux_w.x[i][j][km-3];
                aux_w.x[i][j][0] = aux_w.x[i][j][km-1] =
                    (aux_w.x[i][j][0] + aux_w.x[i][j][km-1])/2.0;

                rhs_u.x[i][j][0] = c43 * rhs_u.x[i][j][1] - c13 * rhs_u.x[i][j][2];
                rhs_u.x[i][j][km-1] = c43 * rhs_u.x[i][j][km-2] - c13 * rhs_u.x[i][j][km-3];
                rhs_u.x[i][j][0] = rhs_u.x[i][j][km-1] =
                    (rhs_u.x[i][j][0] + rhs_u.x[i][j][km-1])/2.0;

                rhs_v.x[i][j][0] = c43 * rhs_v.x[i][j][1] - c13 * rhs_v.x[i][j][2];
                rhs_v.x[i][j][km-1] = c43 * rhs_v.x[i][j][km-2] - c13 * rhs_v.x[i][j][km-3];
                rhs_v.x[i][j][0] = rhs_v.x[i][j][km-1] =
                    (rhs_v.x[i][j][0] + rhs_v.x[i][j][km-1])/2.0;

                rhs_w.x[i][j][0] = c43 * rhs_w.x[i][j][1] - c13 * rhs_w.x[i][j][2];
                rhs_w.x[i][j][km-1] = c43 * rhs_w.x[i][j][km-2] - c13 * rhs_w.x[i][j][km-3];
                rhs_w.x[i][j][0] = rhs_w.x[i][j][km-1] =
                    (rhs_w.x[i][j][0] + rhs_w.x[i][j][km-1])/2.0;
            }
        }

}
