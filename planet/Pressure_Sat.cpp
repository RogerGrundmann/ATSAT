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

// PressureSolverSat::run() now lives in PressureSolverSat.h, which carries the mirrored ATJUP
// solver and dispatches back to computePressure() below unless ATSAT_PRESS_SOLVER is set.
// computePressure() is the inherited solver and stays the default; the header records what the
// two do differently, in particular that the divergence source here is div(grad p) and not
// div(u*), so this routine never sees the velocity divergence at all.

/*
 * ONE place where the mixture density is formed, mirroring ATJUP's computeMixtureDensity().
 *
 * Before this, the ideal-gas relation was written out twice — in cSaturnModel::rho_at() and
 * again as a lambda inside ThermalWindDiagSat.h — with different fallbacks for an unusable
 * cell (r_mix in one, 0.0 in the other). Two copies of one formula can drift; now there is
 * one, and rho_at() is a gated accessor over the array it fills.
 *
 * R_mix is in J/(g K) in this model family, hence the 1e3; p_stat is in bars, hence the 1e5.
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
                const double rho = (p_stat.x[i][j][k] * 1.0e5) / (R_mix * 1.0e3 * T);
                rho_mix.x[i][j][k] = std::isfinite(rho) ? rho : 0.0;
            }
        }
    }
}

void cSaturnModel::computePressure(){
    cout << endl << "      ATSAT: computePressure" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for
    for(int j = 1; j < jm-1; j++){  // r-direction
        for(int k = 1; k < km-1; k++){

            aux_u.x[0][j][k] = aux_u.x[3][j][k] 
                - 3.0 * aux_u.x[2][j][k] + 3.0 * aux_u.x[1][j][k];  // extrapolation
            aux_u.x[im-1][j][k] = aux_u.x[im-4][j][k] 
                - 3.0 * aux_u.x[im-3][j][k] + 3.0 * aux_u.x[im-2][j][k];

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
    for(int k = 1; k < km-1; k++){  // theta-direction
        for(int i = 1; i < im-1; i++){

            aux_u.x[i][0][k] = aux_u.x[i][3][k] 
                - 3.0 * aux_u.x[i][2][k] + 3.0 * aux_u.x[i][1][k];  // extrapolation
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
    for(int i = 0; i < im; i++){  // phi-direction
        for(int j = 0; j < jm; j++){

            aux_u.x[i][j][0] = c43 * aux_u.x[i][j][1] - c13 * aux_u.x[i][j][2];  // von Neumann
            aux_u.x[i][j][km-1] = c43 * aux_u.x[i][j][km-2] - c13 * aux_u.x[i][j][km-3];
            aux_u.x[i][j][0] = aux_u.x[i][j][km-1] = (aux_u.x[i][j][0] + aux_u.x[i][j][km-1])/2.0;

            aux_v.x[i][j][0] = c43 * aux_v.x[i][j][1] - c13 * aux_v.x[i][j][2];
            aux_v.x[i][j][km-1] = c43 * aux_v.x[i][j][km-2] - c13 * aux_v.x[i][j][km-3];
            aux_v.x[i][j][0] = aux_v.x[i][j][km-1] = (aux_v.x[i][j][0] + aux_v.x[i][j][km-1])/2.0;

            aux_w.x[i][j][0] = c43 * aux_w.x[i][j][1] - c13 * aux_w.x[i][j][2];
            aux_w.x[i][j][km-1] = c43 * aux_w.x[i][j][km-2] - c13 * aux_w.x[i][j][km-3];
            aux_w.x[i][j][0] = aux_w.x[i][j][km-1] = (aux_w.x[i][j][0] + aux_w.x[i][j][km-1])/2.0;


            rhs_u.x[i][j][0] = c43 * rhs_u.x[i][j][1] - c13 * rhs_u.x[i][j][2];
            rhs_u.x[i][j][km-1] = c43 * rhs_u.x[i][j][km-2] - c13 * rhs_u.x[i][j][km-3];
            rhs_u.x[i][j][0] = rhs_u.x[i][j][km-1] = (rhs_u.x[i][j][0] + rhs_u.x[i][j][km-1])/2.0;

            rhs_v.x[i][j][0] = c43 * rhs_v.x[i][j][1] - c13 * rhs_v.x[i][j][2];
            rhs_v.x[i][j][km-1] = c43 * rhs_v.x[i][j][km-2] - c13 * rhs_v.x[i][j][km-3];
            rhs_v.x[i][j][0] = rhs_v.x[i][j][km-1] = (rhs_v.x[i][j][0] + rhs_v.x[i][j][km-1])/2.0;

            rhs_w.x[i][j][0] = c43 * rhs_w.x[i][j][1] - c13 * rhs_w.x[i][j][2];
            rhs_w.x[i][j][km-1] = c43 * rhs_w.x[i][j][km-2] - c13 * rhs_w.x[i][j][km-3];
            rhs_w.x[i][j][0] = rhs_w.x[i][j][km-1] = (rhs_w.x[i][j][0] + rhs_w.x[i][j][km-1])/2.0;
        }
    }


    double rm = 0.0;
    double dr2 = dr * dr;
    double dthe2 = dthe * dthe;
    double dphi2 = dphi * dphi;
    double sinthe = 0.0;
    double rmsinthe = 0.0;
    double denom = 0.0;
    double num1 = 0.0;
    double num2 = 0.0;
    double num3 = 0.0;
    double daux_udr = 0.0;
    double daux_vdthe = 0.0;
    double daux_wdphi = 0.0;

    for(int i = 1; i < im-1; i++){
        rm = metricRadius(rad.z[i]);

        for(int j = 1; j < jm-1; j++){
            sinthe = sin(the.z[j]);
            if(sinthe == 0.0) sinthe = 1.0e-5;
            rmsinthe = rm * sinthe;
            denom = 2.0/dr2 + 2.0/(rm * rm * dthe2) 
                + 2.0/(rmsinthe * rmsinthe * dphi2);
            num1 = 1.0/dr2;
            num2 = 1.0/(rm * rm * dthe2);
            num3 = 1.0/(rmsinthe * rmsinthe * dphi2);

            for(int k = 1; k < km-1; k++){
                daux_udr = 
                     ((aux_u.x[i+1][j][k] - aux_u.x[i-1][j][k]) 
                    - (rhs_u.x[i+1][j][k] - rhs_u.x[i-1][j][k]))
                    /(2.0 * dr);
                daux_vdthe = 
                     ((aux_v.x[i][j+1][k] - aux_v.x[i][j-1][k]) 
                    - (rhs_v.x[i][j+1][k] - rhs_v.x[i][j-1][k]))
                    /(2.0 * dthe * rm);
                daux_wdphi = 
                     ((aux_w.x[i][j][k+1] - aux_w.x[i][j][k-1]) 
                    - (rhs_w.x[i][j][k+1] - rhs_w.x[i][j][k-1]))
                    /(2.0 * dphi * rmsinthe);


                p_dyn.x[i][j][k] = 
                     ((p_dyn.x[i+1][j][k] + p_dyn.x[i-1][j][k]) * num1 
                    + (p_dyn.x[i][j+1][k] + p_dyn.x[i][j-1][k]) * num2 
                    + (p_dyn.x[i][j][k+1] + p_dyn.x[i][j][k-1]) * num3 
                    - (daux_udr + daux_vdthe + daux_wdphi))/denom;

            }  // end k
        }  // end j
    }  // end i





    #pragma omp parallel for
    for(int k = 1; k < km-1; k++){  // r-direction
        for(int j = 1; j < jm-1; j++){

            p_dyn.x[0][j][k] = p_dyn.x[3][j][k] 
                - 3.0 * p_dyn.x[2][j][k] 
                + 3.0 * p_dyn.x[1][j][k];  // extrapolation
            p_dyn.x[im-1][j][k] = p_dyn.x[im-4][j][k] 
                - 3.0 * p_dyn.x[im-3][j][k] 
                + 3.0 * p_dyn.x[im-2][j][k];  // extrapolation

        }
    }

    #pragma omp parallel for
    for(int k = 1; k < km-1; k++){  // theta-direction
        for(int i = 1; i < im-1; i++){

            p_dyn.x[i][0][k] = p_dyn.x[i][3][k] 
                - 3.0 * p_dyn.x[i][2][k] + 3.0 * p_dyn.x[i][1][k];  // extrapolation
            p_dyn.x[i][jm-1][k] = p_dyn.x[i][jm-4][k] 
                - 3.0 * p_dyn.x[i][jm-3][k] + 3.0 * p_dyn.x[i][jm-2][k];  // extrapolation

        }
    }

    #pragma omp parallel for
    for(int i = 0; i < im; i++){  // phi-direction
        for(int j = 0; j < jm; j++){

            p_dyn.x[i][j][0] = c43 * p_dyn.x[i][j][1] - c13 * p_dyn.x[i][j][2];  // von Neumann
            p_dyn.x[i][j][km-1] = c43 * p_dyn.x[i][j][km-2] - c13 * p_dyn.x[i][j][km-3];
            p_dyn.x[i][j][0] = p_dyn.x[i][j][km-1] = (p_dyn.x[i][j][0] + p_dyn.x[i][j][km-1])/2.0;

        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for computePressure\n", elapsed.count() * 1e-9);

    cout << "      ATSAT: computePressure ended" << endl;
    return;
}
/*
*
*/
