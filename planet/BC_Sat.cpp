#include "cSaturnModel.h"
#include "BC_Sat.h"

using namespace std;

// ===== Boundary treatment of the turbulence fields =====
// k* and dis* became prognostic variables of the RK4 system in RungeKutta_Sat_Turb.cpp, and
// nothing here had been told about them: the three BC routines below extrapolate every other
// transported field to the domain boundaries and left tke, dis and nue frozen at whatever
// init_fields() had put there. ATJUP records the consequence at BC_Jup.h:190 — the interior
// evolves away from the frozen faces, the resulting permanent Laplacian at i=im-1, the poles and
// the phi seam drains dis* to its floor within ~20 iterations, nue* = k*/omega* then saturates
// its ceiling, and with the coupling on that eddy viscosity destroys the momentum field. So this
// is a prerequisite for ATSAT_TURB_COUPLING, not an improvement on it.
//
// Two departures from the surrounding code, both deliberate:
//
// The 2-point Neumann form f[s] = (4/3)f[a] - (1/3)f[b] is used, not the 3-point cubic
// f[d] - 3f[c] + 3f[b] that ATSAT applies to t and the species. The cubic amplifies an
// alternating error by 7x per call, which the other fields survive and these two do not: dis*
// appears in denominators throughout the closure (nue = k/dis among them). ATJUP switched its
// turbulence fields to the 2-point form for exactly this reason.
//
// The extrapolated values are clamped to the same bounds the RK4 stages enforce — k* >= 0 and
// dis* >= dis_min — because an extrapolation is free to produce a negative where the integration
// is not, and a negative dis* at one face is an infinite nue* at the next call of the closure.
//
// i=0 is also written here, though TurbulenceSat::apply_wall_bc() re-imposes its own zero-gradient
// condition there later in the same iteration. The duplication is what keeps tken/disn — which
// restoreVar copies from tke/dis between the two — from carrying a stale deep boundary.
//
// Gated by cSaturnModel::turb_active — ATSAT_TURB and the configured turb_model resolved into one
// flag. With the closure off tke, dis and nue are identically zero, so the extrapolation would be
// a no-op, and skipping it keeps the off path exactly as it was.
namespace {
    constexpr double bc_dis_min = 1.0e-10;   // matches TurbulenceSat::dis_min and the RK4 floor
}

void BC_Sat::bcRadius() { m.BC_radius(); }
void BC_Sat::bcTheta()  { m.BC_theta(); }
void BC_Sat::bcPhi()    { m.BC_phi(); }

void cSaturnModel::BC_radius(){
//    cout << endl << "      ATSAT: BC_radius" << endl;

//    auto begin = std::chrono::high_resolution_clock::now();

    const bool turb_bc = turb_active;

  #pragma omp parallel for
    for(int j = 1; j < jm-1; j++){
        for(int k = 1; k < km-1; k++){
/*
            u.x[0][j][k] = c43 * u.x[1][j][k] - c13 * u.x[2][j][k];
            v.x[0][j][k] = c43 * v.x[1][j][k] - c13 * v.x[2][j][k];
            w.x[0][j][k] = c43 * w.x[1][j][k] - c13 * w.x[2][j][k];

            u.x[0][j][k] = 0.0;
            v.x[0][j][k] = 0.0;
            w.x[0][j][k] = 0.0;

            u.x[im-1][j][k] = c43 * u.x[im-2][j][k] - c13 * u.x[im-3][j][k];
            v.x[im-1][j][k] = c43 * v.x[im-2][j][k] - c13 * v.x[im-3][j][k];
            w.x[im-1][j][k] = c43 * w.x[im-2][j][k] - c13 * w.x[im-3][j][k];
*/


            t.x[0][j][k] = t.x[3][j][k] 
                - 3.0 * t.x[2][j][k] + 3.0 * t.x[1][j][k];  // extrapolation
            u.x[0][j][k] = u.x[3][j][k] 
                - 3.0 * u.x[2][j][k] + 3.0 * u.x[1][j][k];  // extrapolation
            v.x[0][j][k] = v.x[3][j][k] 
                - 3.0 * v.x[2][j][k] + 3.0 * v.x[1][j][k];  // extrapolation
            w.x[0][j][k] = w.x[3][j][k] 
                - 3.0 * w.x[2][j][k] + 3.0 * w.x[1][j][k];  // extrapolation

            t.x[im-1][j][k] = t.x[im-4][j][k] 
                - 3.0 * t.x[im-3][j][k] + 3.0 * t.x[im-2][j][k];  // extrapolation
            u.x[im-1][j][k] = u.x[im-4][j][k] 
                - 3.0 * u.x[im-3][j][k] + 3.0 * u.x[im-2][j][k];  // extrapolation
            v.x[im-1][j][k] = v.x[im-4][j][k] 
                - 3.0 * v.x[im-3][j][k] + 3.0 * v.x[im-2][j][k];  // extrapolation
            w.x[im-1][j][k] = w.x[im-4][j][k] 
                 - 3.0 * w.x[im-3][j][k] + 3.0 * w.x[im-2][j][k];  // extrapolation


            ch4.x[0][j][k] = ch4.x[3][j][k]
                - 3.0 * ch4.x[2][j][k] + 3.0 * ch4.x[1][j][k];  // extrapolation
            ch4_cloud.x[0][j][k] = ch4_cloud.x[3][j][k]
                - 3.0 * ch4_cloud.x[2][j][k] + 3.0 * ch4_cloud.x[1][j][k];  // extrapolation
            ch4_ice.x[0][j][k] = ch4_ice.x[3][j][k]
                - 3.0 * ch4_ice.x[2][j][k] + 3.0 * ch4_ice.x[1][j][k];  // extrapolation

            ch4.x[im-1][j][k] = ch4.x[im-4][j][k]
                - 3.0 * ch4.x[im-3][j][k] + 3.0 * ch4.x[im-2][j][k];  // extrapolation
            ch4_cloud.x[im-1][j][k] = ch4_cloud.x[im-4][j][k]
                - 3.0 * ch4_cloud.x[im-3][j][k] + 3.0 * ch4_cloud.x[im-2][j][k];  // extrapolation
            ch4_ice.x[im-1][j][k] = ch4_ice.x[im-4][j][k]
                 - 3.0 * ch4_ice.x[im-3][j][k] + 3.0 * ch4_ice.x[im-2][j][k];  // extrapolation


            h2o.x[0][j][k] = h2o.x[3][j][k]
                - 3.0 * h2o.x[2][j][k] + 3.0 * h2o.x[1][j][k];  // extrapolation
            h2o_cloud.x[0][j][k] = h2o_cloud.x[3][j][k]
                - 3.0 * h2o_cloud.x[2][j][k] + 3.0 * h2o_cloud.x[1][j][k];  // extrapolation
            h2o_ice.x[0][j][k] = h2o_ice.x[3][j][k]
                - 3.0 * h2o_ice.x[2][j][k] + 3.0 * h2o_ice.x[1][j][k];  // extrapolation

            h2o.x[im-1][j][k] = h2o.x[im-4][j][k] 
                - 3.0 * h2o.x[im-3][j][k] + 3.0 * h2o.x[im-2][j][k];  // extrapolation
            h2o_cloud.x[im-1][j][k] = h2o_cloud.x[im-4][j][k] 
                - 3.0 * h2o_cloud.x[im-3][j][k] + 3.0 * h2o_cloud.x[im-2][j][k];  // extrapolation
            h2o_ice.x[im-1][j][k] = h2o_ice.x[im-4][j][k] 
                 - 3.0 * h2o_ice.x[im-3][j][k] + 3.0 * h2o_ice.x[im-2][j][k];  // extrapolation


            h2s.x[0][j][k] = h2s.x[3][j][k] 
                - 3.0 * h2s.x[2][j][k] + 3.0 * h2s.x[1][j][k];  // extrapolation
            h2s.x[im-1][j][k] = h2s.x[im-4][j][k] 
                - 3.0 * h2s.x[im-3][j][k] + 3.0 * h2s.x[im-2][j][k];  // extrapolation


            nh3.x[0][j][k] = nh3.x[3][j][k] 
                - 3.0 * nh3.x[2][j][k] + 3.0 * nh3.x[1][j][k];  // extrapolation
            nh3_cloud.x[0][j][k] = nh3_cloud.x[3][j][k] 
                - 3.0 * nh3_cloud.x[2][j][k] + 3.0 * nh3_cloud.x[1][j][k];  // extrapolation
            nh3_ice.x[0][j][k] = nh3_ice.x[3][j][k] 
                - 3.0 * nh3_ice.x[2][j][k] + 3.0 * nh3_ice.x[1][j][k];  // extrapolation

            nh3.x[im-1][j][k] = nh3.x[im-4][j][k] 
                - 3.0 * nh3.x[im-3][j][k] + 3.0 * nh3.x[im-2][j][k];  // extrapolation
            nh3_cloud.x[im-1][j][k] = nh3_cloud.x[im-4][j][k] 
                - 3.0 * nh3_cloud.x[im-3][j][k] + 3.0 * nh3_cloud.x[im-2][j][k];  // extrapolation
            nh3_ice.x[im-1][j][k] = nh3_ice.x[im-4][j][k] 
                 - 3.0 * nh3_ice.x[im-3][j][k] + 3.0 * nh3_ice.x[im-2][j][k];  // extrapolation


            nh4sh.x[0][j][k] = nh4sh.x[3][j][k] 
                - 3.0 * nh4sh.x[2][j][k] + 3.0 * nh4sh.x[1][j][k];  // extrapolation
            nh4sh.x[im-1][j][k] = nh4sh.x[im-4][j][k] 
                - 3.0 * nh4sh.x[im-3][j][k] + 3.0 * nh4sh.x[im-2][j][k];  // extrapolation



/*
            h2o.x[0][j][k] = c43 * h2o.x[1][j][k] - c13 * h2o.x[2][j][k];
            h2o_cloud.x[0][j][k] = c43 * h2o_cloud.x[1][j][k] - c13 * h2o_cloud.x[2][j][k];
            h2o_ice.x[0][j][k] = c43 * h2o_ice.x[1][j][k] - c13 * h2o_ice.x[2][j][k];

            h2o.x[im-1][j][k] = c43 * h2o.x[im-2][j][k] - c13 * h2o.x[im-3][j][k];
            h2o_cloud.x[im-1][j][k] = c43 * h2o_cloud.x[im-2][j][k] - c13 * h2o_cloud.x[im-3][j][k];
            h2o_ice.x[im-1][j][k] = c43 * h2o_ice.x[im-2][j][k] - c13 * h2o_ice.x[im-3][j][k];

            h2s.x[0][j][k] = c43 * h2s.x[1][j][k] - c13 * h2s.x[2][j][k];
            h2s_cloud.x[0][j][k] = c43 * h2s_cloud.x[1][j][k] - c13 * h2s_cloud.x[2][j][k];
            h2s_ice.x[0][j][k] = c43 * h2s_ice.x[1][j][k] - c13 * h2s_ice.x[2][j][k];

            h2s.x[im-1][j][k] = c43 * h2s.x[im-2][j][k] - c13 * h2s.x[im-3][j][k];
            h2s_cloud.x[im-1][j][k] = c43 * h2s_cloud.x[im-2][j][k] - c13 * h2s_cloud.x[im-3][j][k];
            h2s_ice.x[im-1][j][k] = c43 * h2s_ice.x[im-2][j][k] - c13 * h2s_ice.x[im-3][j][k];

            nh3.x[0][j][k] = c43 * nh3.x[1][j][k] - c13 * nh3.x[2][j][k];
            nh3_cloud.x[0][j][k] = c43 * nh3_cloud.x[1][j][k] - c13 * nh3.x[2][j][k];
            nh3_ice.x[0][j][k] = c43 * nh3_ice.x[1][j][k] - c13 * nh3_ice.x[2][j][k];

            nh3.x[im-1][j][k] = c43 * nh3.x[im-2][j][k] - c13 * nh3.x[im-3][j][k];
            nh3_cloud.x[im-1][j][k] = c43 * nh3_cloud.x[im-2][j][k] - c13 * nh3_cloud.x[im-3][j][k];
            nh3_ice.x[im-1][j][k] = c43 * nh3_ice.x[im-2][j][k] - c13 * nh3_ice.x[im-3][j][k];


            nh4sh.x[0][j][k] = c43 * nh4sh.x[1][j][k] - c13 * nh4sh.x[2][j][k];
//            nh4sh_cloud.x[0][j][k] = c43 * nh4sh_cloud.x[1][j][k] - c13 * nh4sh.x[2][j][k];
//            nh4sh_ice.x[0][j][k] = c43 * nh4sh_ice.x[1][j][k] - c13 * nh4sh_ice.x[2][j][k];

            nh4sh.x[im-1][j][k] = c43 * nh4sh.x[im-2][j][k] - c13 * nh4sh.x[im-3][j][k];
//            nh4sh_cloud.x[im-1][j][k] = c43 * nh4sh_cloud.x[im-2][j][k] - c13 * nh4sh_cloud.x[im-3][j][k];
//            nh4sh_ice.x[im-1][j][k] = c43 * nh4sh_ice.x[im-2][j][k] - c13 * nh4sh_ice.x[im-3][j][k];

//            h2.x[0][j][k] = c43 * h2.x[1][j][k] - c13 * nh4sh.x[2][j][k];
//            he.x[0][j][k] = c43 * he.x[1][j][k] - c13 * he.x[2][j][k];

//            h2.x[im-1][j][k] = c43 * h2.x[im-2][j][k] - c13 * h2.x[im-3][j][k];
//            he.x[im-1][j][k] = c43 * he.x[im-2][j][k] - c13 * he.x[im-3][j][k];
*/


            j_h2s.x[0][j][k] = j_h2s.x[3][j][k] 
                - 3.0 * j_h2s.x[2][j][k] + 3.0 * j_h2s.x[1][j][k];  // extrapolation
            j_h2s.x[im-1][j][k] = j_h2s.x[im-4][j][k] 
                - 3.0 * j_h2s.x[im-3][j][k] + 3.0 * j_h2s.x[im-2][j][k];  // extrapolation


            j_nh3.x[0][j][k] = j_nh3.x[3][j][k] 
                - 3.0 * j_nh3.x[2][j][k] + 3.0 * j_nh3.x[1][j][k];  // extrapolation
            j_nh3.x[im-1][j][k] = j_nh3.x[im-4][j][k] 
                - 3.0 * j_nh3.x[im-3][j][k] + 3.0 * j_nh3.x[im-2][j][k];  // extrapolation


            j_nh4sh.x[0][j][k] = j_nh4sh.x[3][j][k] 
                - 3.0 * j_nh4sh.x[2][j][k] + 3.0 * j_nh4sh.x[1][j][k];  // extrapolation
            j_nh4sh.x[im-1][j][k] = j_nh4sh.x[im-4][j][k] 
                - 3.0 * j_nh4sh.x[im-3][j][k] + 3.0 * j_nh4sh.x[im-2][j][k];  // extrapolation


            jT_h2s.x[0][j][k] = jT_h2s.x[3][j][k] 
                - 3.0 * jT_h2s.x[2][j][k] + 3.0 * jT_h2s.x[1][j][k];  // extrapolation
            jT_h2s.x[im-1][j][k] = jT_h2s.x[im-4][j][k] 
                - 3.0 * jT_h2s.x[im-3][j][k] + 3.0 * jT_h2s.x[im-2][j][k];  // extrapolation


            jT_nh3.x[0][j][k] = jT_nh3.x[3][j][k] 
                - 3.0 * jT_nh3.x[2][j][k] + 3.0 * jT_nh3.x[1][j][k];  // extrapolation
            jT_nh3.x[im-1][j][k] = jT_nh3.x[im-4][j][k] 
                - 3.0 * jT_nh3.x[im-3][j][k] + 3.0 * jT_nh3.x[im-2][j][k];  // extrapolation


            jT_nh4sh.x[0][j][k] = jT_nh4sh.x[3][j][k] 
                - 3.0 * jT_nh4sh.x[2][j][k] + 3.0 * jT_nh4sh.x[1][j][k];  // extrapolation
            jT_nh4sh.x[im-1][j][k] = jT_nh4sh.x[im-4][j][k] 
                - 3.0 * jT_nh4sh.x[im-3][j][k] + 3.0 * jT_nh4sh.x[im-2][j][k];  // extrapolation


/*
            j_h2s.x[0][j][k] = c43 * j_h2s.x[1][j][k] - c13 * j_h2s.x[2][j][k];
            j_h2s.x[im-1][j][k] = c43 * j_h2s.x[im-2][j][k] - c13 * j_h2s.x[im-3][j][k];

            j_nh3.x[0][j][k] = c43 * j_nh3.x[1][j][k] - c13 * j_nh3.x[2][j][k];
            j_nh3.x[im-1][j][k] = c43 * j_nh3.x[im-2][j][k] - c13 * j_nh3.x[im-3][j][k];

            j_nh4sh.x[0][j][k] = c43 * j_nh4sh.x[1][j][k] - c13 * j_nh4sh.x[2][j][k];
            j_nh4sh.x[im-1][j][k] = c43 * j_nh4sh.x[im-2][j][k] - c13 * j_nh4sh.x[im-3][j][k];

            jT_h2s.x[0][j][k] = c43 * jT_h2s.x[1][j][k] - c13 * jT_h2s.x[2][j][k];
            jT_h2s.x[im-1][j][k] = c43 * jT_h2s.x[im-2][j][k] - c13 * jT_h2s.x[im-3][j][k];

            jT_nh3.x[0][j][k] = c43 * jT_nh3.x[1][j][k] - c13 * jT_nh3.x[2][j][k];
            jT_nh3.x[im-1][j][k] = c43 * jT_nh3.x[im-2][j][k] - c13 * jT_nh3.x[im-3][j][k];

            jT_nh4sh.x[0][j][k] = c43 * jT_nh4sh.x[1][j][k] - c13 * jT_nh4sh.x[2][j][k];
            jT_nh4sh.x[im-1][j][k] = c43 * jT_nh4sh.x[im-2][j][k] - c13 * jT_nh4sh.x[im-3][j][k];
*/


            w_h2s.x[0][j][k] = w_h2s.x[3][j][k] 
                - 3.0 * w_h2s.x[2][j][k] + 3.0 * w_h2s.x[1][j][k];  // extrapolation
            w_h2s.x[im-1][j][k] = w_h2s.x[im-4][j][k] 
                - 3.0 * w_h2s.x[im-3][j][k] + 3.0 * w_h2s.x[im-2][j][k];  // extrapolation

            w_nh3.x[0][j][k] = w_nh3.x[3][j][k] 
                - 3.0 * w_nh3.x[2][j][k] + 3.0 * w_nh3.x[1][j][k];  // extrapolation
            w_nh3.x[im-1][j][k] = w_nh3.x[im-4][j][k] 
                - 3.0 * w_nh3.x[im-3][j][k] + 3.0 * w_nh3.x[im-2][j][k];  // extrapolation

            w_nh4sh.x[0][j][k] = w_nh4sh.x[3][j][k] 
                - 3.0 * w_nh4sh.x[2][j][k] + 3.0 * w_nh4sh.x[1][j][k];  // extrapolation
            w_nh4sh.x[im-1][j][k] = w_nh4sh.x[im-4][j][k] 
                - 3.0 * w_nh4sh.x[im-3][j][k] + 3.0 * w_nh4sh.x[im-2][j][k];  // extrapolation


            massflux_h2s.x[0][j][k] = massflux_h2s.x[3][j][k] 
                - 3.0 * massflux_h2s.x[2][j][k] + 3.0 * massflux_h2s.x[1][j][k];  // extrapolation
            massflux_h2s.x[im-1][j][k] = massflux_h2s.x[im-4][j][k] 
                - 3.0 * massflux_h2s.x[im-3][j][k] + 3.0 * massflux_h2s.x[im-2][j][k];  // extrapolation

            massflux_nh3.x[0][j][k] = massflux_nh3.x[3][j][k] 
                - 3.0 * massflux_nh3.x[2][j][k] + 3.0 * massflux_nh3.x[1][j][k];  // extrapolation
            massflux_nh3.x[im-1][j][k] = massflux_nh3.x[im-4][j][k] 
                - 3.0 * massflux_nh3.x[im-3][j][k] + 3.0 * massflux_nh3.x[im-2][j][k];  // extrapolation

            massflux_nh4sh.x[0][j][k] = massflux_nh4sh.x[3][j][k] 
                - 3.0 * massflux_nh4sh.x[2][j][k] + 3.0 * massflux_nh4sh.x[1][j][k];  // extrapolation
            massflux_nh4sh.x[im-1][j][k] = massflux_nh4sh.x[im-4][j][k] 
                - 3.0 * massflux_nh4sh.x[im-3][j][k] + 3.0 * massflux_nh4sh.x[im-2][j][k];  // extrapolation



            difflux_h2s.x[0][j][k] = difflux_h2s.x[3][j][k] 
                - 3.0 * difflux_h2s.x[2][j][k] + 3.0 * difflux_h2s.x[1][j][k];  // extrapolation
            difflux_h2s.x[im-1][j][k] = difflux_h2s.x[im-4][j][k] 
                - 3.0 * difflux_h2s.x[im-3][j][k] + 3.0 * difflux_h2s.x[im-2][j][k];  // extrapolation

            difflux_nh3.x[0][j][k] = difflux_nh3.x[3][j][k] 
                - 3.0 * difflux_nh3.x[2][j][k] + 3.0 * difflux_nh3.x[1][j][k];  // extrapolation
            difflux_nh3.x[im-1][j][k] = difflux_nh3.x[im-4][j][k] 
                - 3.0 * difflux_nh3.x[im-3][j][k] + 3.0 * difflux_nh3.x[im-2][j][k];  // extrapolation

            difflux_nh4sh.x[0][j][k] = difflux_nh4sh.x[3][j][k]
                - 3.0 * difflux_nh4sh.x[2][j][k] + 3.0 * difflux_nh4sh.x[1][j][k];  // extrapolation
            difflux_nh4sh.x[im-1][j][k] = difflux_nh4sh.x[im-4][j][k]
                - 3.0 * difflux_nh4sh.x[im-3][j][k] + 3.0 * difflux_nh4sh.x[im-2][j][k];  // extrapolation

            fluxlim_nh4sh.x[0][j][k] = fluxlim_nh4sh.x[3][j][k]
                - 3.0 * fluxlim_nh4sh.x[2][j][k] + 3.0 * fluxlim_nh4sh.x[1][j][k];  // extrapolation
            fluxlim_nh4sh.x[im-1][j][k] = fluxlim_nh4sh.x[im-4][j][k]
                - 3.0 * fluxlim_nh4sh.x[im-3][j][k] + 3.0 * fluxlim_nh4sh.x[im-2][j][k];  // extrapolation


            thermalmassflux.x[0][j][k] = thermalmassflux.x[3][j][k]
                - 3.0 * thermalmassflux.x[2][j][k] + 3.0 * thermalmassflux.x[1][j][k];  // extrapolation
            thermalmassflux.x[im-1][j][k] = thermalmassflux.x[im-4][j][k] 
                - 3.0 * thermalmassflux.x[im-3][j][k] + 3.0 * thermalmassflux.x[im-2][j][k];  // extrapolation

            CoriolisForce.x[0][j][k] = CoriolisForce.x[3][j][k] 
                - 3.0 * CoriolisForce.x[2][j][k] + 3.0 * CoriolisForce.x[1][j][k];  // extrapolation
            CoriolisForce.x[im-1][j][k] = CoriolisForce.x[im-4][j][k] 
                - 3.0 * CoriolisForce.x[im-3][j][k] + 3.0 * CoriolisForce.x[im-2][j][k];  // extrapolation

            CentrifugalForce.x[0][j][k] = CentrifugalForce.x[3][j][k] 
                - 3.0 * CentrifugalForce.x[2][j][k] + 3.0 * CentrifugalForce.x[1][j][k];  // extrapolation
            CentrifugalForce.x[im-1][j][k] = CentrifugalForce.x[im-4][j][k] 
                - 3.0 * CentrifugalForce.x[im-3][j][k] + 3.0 * CentrifugalForce.x[im-2][j][k];  // extrapolation

            BuoyancyForce.x[0][j][k] = BuoyancyForce.x[3][j][k] 
                - 3.0 * BuoyancyForce.x[2][j][k] + 3.0 * BuoyancyForce.x[1][j][k];  // extrapolation
            BuoyancyForce.x[im-1][j][k] = BuoyancyForce.x[im-4][j][k] 
                - 3.0 * BuoyancyForce.x[im-3][j][k] + 3.0 * BuoyancyForce.x[im-2][j][k];  // extrapolation

            PresGradForce.x[0][j][k] = PresGradForce.x[3][j][k] 
                - 3.0 * PresGradForce.x[2][j][k] + 3.0 * PresGradForce.x[1][j][k];  // extrapolation
            PresGradForce.x[im-1][j][k] = PresGradForce.x[im-4][j][k] 
                - 3.0 * PresGradForce.x[im-3][j][k] + 3.0 * PresGradForce.x[im-2][j][k];  // extrapolation

            Q_Latent.x[0][j][k] = Q_Latent.x[3][j][k] 
                - 3.0 * Q_Latent.x[2][j][k] + 3.0 * Q_Latent.x[1][j][k];  // extrapolation
            Q_Latent.x[im-1][j][k] = Q_Latent.x[im-4][j][k] 
                - 3.0 * Q_Latent.x[im-3][j][k] + 3.0 * Q_Latent.x[im-2][j][k];  // extrapolation

            Q_Sensible.x[0][j][k] = Q_Sensible.x[3][j][k] 
                - 3.0 * Q_Sensible.x[2][j][k] + 3.0 * Q_Sensible.x[1][j][k];  // extrapolation
            Q_Sensible.x[im-1][j][k] = Q_Sensible.x[im-4][j][k] 
                - 3.0 * Q_Sensible.x[im-3][j][k] + 3.0 * Q_Sensible.x[im-2][j][k];  // extrapolation


            // k*, dis*, nue* at the deep boundary and the model top (see the note at the top
            // of this file for why the form and the clamps differ from the fields above).
            if(turb_bc){
                tke.x[0][j][k] = std::max(0.0,
                    c43 * tke.x[1][j][k] - c13 * tke.x[2][j][k]);
                dis.x[0][j][k] = std::max(bc_dis_min,
                    c43 * dis.x[1][j][k] - c13 * dis.x[2][j][k]);
                nue.x[0][j][k] = std::max(0.0,
                    c43 * nue.x[1][j][k] - c13 * nue.x[2][j][k]);

                tke.x[im-1][j][k] = std::max(0.0,
                    c43 * tke.x[im-2][j][k] - c13 * tke.x[im-3][j][k]);
                dis.x[im-1][j][k] = std::max(bc_dis_min,
                    c43 * dis.x[im-2][j][k] - c13 * dis.x[im-3][j][k]);
                nue.x[im-1][j][k] = std::max(0.0,
                    c43 * nue.x[im-2][j][k] - c13 * nue.x[im-3][j][k]);
            }


/*
            w_h2s.x[0][j][k] = c43 * w_h2s.x[1][j][k] - c13 * w_h2s.x[2][j][k];
            w_h2s.x[im-1][j][k] = c43 * w_h2s.x[im-2][j][k] - c13 * w_h2s.x[im-3][j][k];

            w_nh3.x[0][j][k] = c43 * w_nh3.x[1][j][k] - c13 * w_nh3.x[2][j][k];
            w_nh3.x[im-1][j][k] = c43 * w_nh3.x[im-2][j][k] - c13 * w_nh3.x[im-3][j][k];

            w_nh4sh.x[0][j][k] = c43 * w_nh4sh.x[1][j][k] - c13 * w_nh4sh.x[2][j][k];
            w_nh4sh.x[im-1][j][k] = c43 * w_nh4sh.x[im-2][j][k] - c13 * w_nh4sh.x[im-3][j][k];


            massflux_h2s.x[0][j][k] = c43 * massflux_h2s.x[1][j][k] - c13 * massflux_h2s.x[2][j][k];
            massflux_h2s.x[im-1][j][k] = c43 * massflux_h2s.x[im-2][j][k] - c13 * massflux_h2s.x[im-3][j][k];

            massflux_nh3.x[0][j][k] = c43 * massflux_nh3.x[1][j][k] - c13 * massflux_nh3.x[2][j][k];
            massflux_nh3.x[im-1][j][k] = c43 * massflux_nh3.x[im-2][j][k] - c13 * massflux_nh3.x[im-3][j][k];

            massflux_nh4sh.x[0][j][k] = c43 * massflux_nh4sh.x[1][j][k] - c13 * massflux_nh4sh.x[2][j][k];
            massflux_nh4sh.x[im-1][j][k] = c43 * massflux_nh4sh.x[im-2][j][k] - c13 * massflux_nh4sh.x[im-3][j][k];

            difflux_h2s.x[0][j][k] = c43 * difflux_h2s.x[1][j][k] - c13 * difflux_h2s.x[2][j][k];
            difflux_h2s.x[im-1][j][k] = c43 * difflux_h2s.x[im-2][j][k] - c13 * difflux_h2s.x[im-3][j][k];

            difflux_nh3.x[0][j][k] = c43 * difflux_nh3.x[1][j][k] - c13 * difflux_nh3.x[2][j][k];
            difflux_nh3.x[im-1][j][k] = c43 * difflux_nh3.x[im-2][j][k] - c13 * difflux_nh3.x[im-3][j][k];

            difflux_nh4sh.x[0][j][k] = c43 * difflux_nh4sh.x[1][j][k] - c13 * difflux_nh4sh.x[2][j][k];
            difflux_nh4sh.x[im-1][j][k] = c43 * difflux_nh4sh.x[im-2][j][k] - c13 * difflux_nh4sh.x[im-3][j][k];

            thermalmassflux.x[0][j][k] = c43 * thermalmassflux.x[1][j][k] - c13 * thermalmassflux.x[2][j][k];
            thermalmassflux.x[im-1][j][k] = c43 * thermalmassflux.x[im-2][j][k] - c13 * thermalmassflux.x[im-3][j][k];

            CoriolisForce.x[0][j][k] = c43 * CoriolisForce.x[1][j][k] - c13 * CoriolisForce.x[2][j][k];
            CoriolisForce.x[im-1][j][k] = c43 * CoriolisForce.x[im-2][j][k] - c13 * CoriolisForce.x[im-3][j][k];

            BuoyancyForce.x[0][j][k] = c43 * BuoyancyForce.x[1][j][k] - c13 * BuoyancyForce.x[2][j][k];
            BuoyancyForce.x[im-1][j][k] = c43 * BuoyancyForce.x[im-2][j][k] - c13 * BuoyancyForce.x[im-3][j][k];

            PresGradForce.x[0][j][k] = c43 * PresGradForce.x[1][j][k] - c13 * PresGradForce.x[2][j][k];
            PresGradForce.x[im-1][j][k] = c43 * PresGradForce.x[im-2][j][k] - c13 * PresGradForce.x[im-3][j][k];
*/
        }
    }

//    auto end = std::chrono::high_resolution_clock::now();
//    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
//    printf(" time measured: %.3f seconds for BC_radius\n", elapsed.count() * 1e-9);

//    cout << "      ATSAT: BC_radius ended" << endl;
    return;
}


void cSaturnModel::BC_theta(){
//    cout << endl << "      ATSAT: BC_theta" << endl;

//    auto begin = std::chrono::high_resolution_clock::now();

    const bool turb_bc = turb_active;

    #pragma omp parallel for
    for(int k = 1; k < km-1; k++){
        for(int i = 1; i < im-1; i++){
/*
            t.x[i][0][k] = c43 * t.x[i][1][k] - c13 * t.x[i][2][k];
            t.x[i][jm-1][k] = c43 * t.x[i][jm-2][k] - c13 * t.x[i][jm-3][k];

            u.x[i][0][k] = c43 * u.x[i][1][k] - c13 * u.x[i][2][k];
            u.x[i][jm-1][k] = c43 * u.x[i][jm-2][k] - c13 * u.x[i][jm-3][k];

            v.x[i][0][k] = c43 * v.x[i][1][k] - c13 * v.x[i][2][k];
            v.x[i][jm-1][k] = c43 * v.x[i][jm-2][k] - c13 * v.x[i][jm-3][k];

            w.x[i][0][k] = c43 * w.x[i][1][k] - c13 * w.x[i][2][k];
            w.x[i][jm-1][k] = c43 * w.x[i][jm-2][k] - c13 * w.x[i][jm-3][k];
*/

//            u.x[i][0][k] = 0.0;
//            u.x[i][jm-1][k] = 0.0;
            v.x[i][0][k] = 0.0;
            v.x[i][jm-1][k] = 0.0;
            w.x[i][0][k] = 0.0;
            w.x[i][jm-1][k] = 0.0;


            t.x[i][0][k] = t.x[i][3][k] 
                - 3.0 * t.x[i][2][k] + 3.0 * t.x[i][1][k];  // extrapolation
            u.x[i][0][k] = u.x[i][3][k] 
                - 3.0 * u.x[i][2][k] + 3.0 * u.x[i][1][k];  // extrapolation
/*
            v.x[i][0][k] = v.x[3][3][k] 
                - 3.0 * v.x[i][2][k] + 3.0 * v.x[i][1][k];  // extrapolation
            w.x[i][0][k] = w.x[3][3][k] 
                - 3.0 * w.x[i][2][k] + 3.0 * w.x[i][1][k];  // extrapolation
*/
            t.x[i][jm-1][k] = t.x[i][jm-4][k] 
                - 3.0 * t.x[i][jm-3][k] + 3.0 * t.x[i][jm-2][k];  // extrapolation
            u.x[i][jm-1][k] = u.x[i][jm-4][k] 
                - 3.0 * u.x[i][jm-3][k] + 3.0 * u.x[i][jm-2][k];  // extrapolation
/*
            v.x[i][jm-1][k] = v.x[i][jm-4][k] 
                - 3.0 * v.x[i][jm-3][k] + 3.0 * v.x[i][jm-2][k];  // extrapolation
            w.x[i][jm-1][k] = w.x[i][jm-4][k] 
                 - 3.0 * w.x[i][jm-3][k] + 3.0 * w.x[i][jm-2][k];  // extrapolation
*/

            ch4.x[i][0][k] = ch4.x[i][3][k]
                - 3.0 * ch4.x[i][2][k] + 3.0 * ch4.x[i][1][k];  // extrapolation
            ch4.x[i][jm-1][k] = ch4.x[i][jm-4][k]
                - 3.0 * ch4.x[i][jm-3][k] + 3.0 * ch4.x[i][jm-2][k];  // extrapolation

            ch4_cloud.x[i][0][k] = ch4_cloud.x[i][3][k]
                - 3.0 * ch4_cloud.x[i][2][k] + 3.0 * ch4_cloud.x[i][1][k];  // extrapolation
            ch4_cloud.x[i][jm-1][k] = ch4_cloud.x[i][jm-4][k]
                - 3.0 * ch4_cloud.x[i][jm-3][k] + 3.0 * ch4_cloud.x[i][jm-2][k];  // extrapolation

            ch4_ice.x[i][0][k] = ch4_ice.x[i][3][k]
                - 3.0 * ch4_ice.x[i][2][k] + 3.0 * ch4_ice.x[i][1][k];  // extrapolation
            ch4_ice.x[i][jm-1][k] = ch4_ice.x[i][jm-4][k]
                - 3.0 * ch4_ice.x[i][jm-3][k] + 3.0 * ch4_ice.x[i][jm-2][k];  // extrapolation


            h2o.x[i][0][k] = h2o.x[i][3][k]
                - 3.0 * h2o.x[i][2][k] + 3.0 * h2o.x[i][1][k];  // extrapolation
            h2o.x[i][jm-1][k] = h2o.x[i][jm-4][k]
                - 3.0 * h2o.x[i][jm-3][k] + 3.0 * h2o.x[i][jm-2][k];  // extrapolation

            h2o_cloud.x[i][0][k] = h2o_cloud.x[i][3][k] 
                - 3.0 * h2o_cloud.x[i][2][k] + 3.0 * h2o_cloud.x[i][1][k];  // extrapolation
            h2o_cloud.x[i][jm-1][k] = h2o_cloud.x[i][jm-4][k] 
                - 3.0 * h2o_cloud.x[i][jm-3][k] + 3.0 * h2o_cloud.x[i][jm-2][k];  // extrapolation

            h2o_ice.x[i][0][k] = h2o_ice.x[i][3][k] 
                - 3.0 * h2o_ice.x[i][2][k] + 3.0 * h2o_ice.x[i][1][k];  // extrapolation
            h2o_ice.x[i][jm-1][k] = h2o_ice.x[i][jm-4][k] 
                - 3.0 * h2o_ice.x[i][jm-3][k] + 3.0 * h2o_ice.x[i][jm-2][k];  // extrapolation


            h2s.x[i][0][k] = h2s.x[i][3][k] 
                - 3.0 * h2s.x[i][2][k] + 3.0 * h2s.x[i][1][k];  // extrapolation
            h2s.x[i][jm-1][k] = h2s.x[i][jm-4][k] 
                - 3.0 * h2s.x[i][jm-3][k] + 3.0 * h2s.x[i][jm-2][k];  // extrapolation


            nh3.x[i][0][k] = nh3.x[i][3][k] 
                - 3.0 * nh3.x[i][2][k] + 3.0 * nh3.x[i][1][k];  // extrapolation
            nh3.x[i][jm-1][k] = nh3.x[i][jm-4][k] 
                - 3.0 * nh3.x[i][jm-3][k] + 3.0 * nh3.x[i][jm-2][k];  // extrapolation

            nh3_cloud.x[i][0][k] = nh3_cloud.x[i][3][k] 
                - 3.0 * nh3_cloud.x[i][2][k] + 3.0 * nh3_cloud.x[i][1][k];  // extrapolation
            nh3_cloud.x[i][jm-1][k] = nh3_cloud.x[i][jm-4][k] 
                - 3.0 * nh3_cloud.x[i][jm-3][k] + 3.0 * nh3_cloud.x[i][jm-2][k];  // extrapolation

            nh3_ice.x[i][0][k] = nh3_ice.x[i][3][k] 
                - 3.0 * nh3_ice.x[i][2][k] + 3.0 * nh3_ice.x[i][1][k];  // extrapolation
            nh3_ice.x[i][jm-1][k] = nh3_ice.x[i][jm-4][k] 
                - 3.0 * nh3_ice.x[i][jm-3][k] + 3.0 * nh3_ice.x[i][jm-2][k];  // extrapolation


            nh4sh.x[i][0][k] = nh4sh.x[i][3][k] 
                - 3.0 * nh4sh.x[i][2][k] + 3.0 * nh4sh.x[i][1][k];  // extrapolation
            nh4sh.x[i][jm-1][k] = nh4sh.x[i][jm-4][k] 
                - 3.0 * nh4sh.x[i][jm-3][k] + 3.0 * nh4sh.x[i][jm-2][k];  // extrapolation


            j_h2s.x[i][0][k] = j_h2s.x[i][3][k] 
                - 3.0 * j_h2s.x[i][2][k] + 3.0 * j_h2s.x[i][1][k];  // extrapolation
            j_h2s.x[i][jm-1][k] = j_h2s.x[i][jm-4][k] 
                - 3.0 * j_h2s.x[i][jm-3][k] + 3.0 * j_h2s.x[i][jm-2][k];  // extrapolation

            j_nh3.x[i][0][k] = j_nh3.x[i][3][k] 
                - 3.0 * j_nh3.x[i][2][k] + 3.0 * j_nh3.x[i][1][k];  // extrapolation
            j_nh3.x[i][jm-1][k] = j_nh3.x[i][jm-4][k] 
                - 3.0 * j_nh3.x[i][jm-3][k] + 3.0 * j_nh3.x[i][jm-2][k];  // extrapolation

            j_nh4sh.x[i][0][k] = j_nh4sh.x[i][3][k] 
                - 3.0 * j_nh4sh.x[i][2][k] + 3.0 * j_nh4sh.x[i][1][k];  // extrapolation
            j_nh4sh.x[i][jm-1][k] = j_nh4sh.x[i][jm-4][k] 
                - 3.0 * j_nh4sh.x[i][jm-3][k] + 3.0 * j_nh4sh.x[i][jm-2][k];  // extrapolation


            jT_h2s.x[i][0][k] = jT_h2s.x[i][3][k] 
                - 3.0 * jT_h2s.x[i][2][k] + 3.0 * jT_h2s.x[i][1][k];  // extrapolation
            jT_h2s.x[i][jm-1][k] = jT_h2s.x[i][jm-4][k] 
                - 3.0 * jT_h2s.x[i][jm-3][k] + 3.0 * jT_h2s.x[i][jm-2][k];  // extrapolation

            jT_nh3.x[i][0][k] = jT_nh3.x[i][3][k] 
                - 3.0 * jT_nh3.x[i][2][k] + 3.0 * jT_nh3.x[i][1][k];  // extrapolation
            jT_nh3.x[i][jm-1][k] = jT_nh3.x[i][jm-4][k] 
                - 3.0 * jT_nh3.x[i][jm-3][k] + 3.0 * jT_nh3.x[i][jm-2][k];  // extrapolation

            jT_nh4sh.x[i][0][k] = jT_nh4sh.x[i][3][k] 
                - 3.0 * jT_nh4sh.x[i][2][k] + 3.0 * jT_nh4sh.x[i][1][k];  // extrapolation
            jT_nh4sh.x[i][jm-1][k] = jT_nh4sh.x[i][jm-4][k] 
                - 3.0 * jT_nh4sh.x[i][jm-3][k] + 3.0 * jT_nh4sh.x[i][jm-2][k];  // extrapolation


            w_h2s.x[i][0][k] = w_h2s.x[i][3][k] 
                - 3.0 * w_h2s.x[i][2][k] + 3.0 * w_h2s.x[i][1][k];  // extrapolation
            w_h2s.x[i][jm-1][k] = w_h2s.x[i][jm-4][k] 
                - 3.0 * w_h2s.x[i][jm-3][k] + 3.0 * w_h2s.x[i][jm-2][k];  // extrapolation

            w_nh3.x[i][0][k] = w_nh3.x[i][3][k] 
                - 3.0 * w_nh3.x[i][2][k] + 3.0 * w_nh3.x[i][1][k];  // extrapolation
            w_nh3.x[i][jm-1][k] = w_nh3.x[i][jm-4][k] 
                - 3.0 * w_nh3.x[i][jm-3][k] + 3.0 * w_nh3.x[i][jm-2][k];  // extrapolation

            w_nh4sh.x[i][0][k] = w_nh4sh.x[i][3][k] 
                - 3.0 * w_nh4sh.x[i][2][k] + 3.0 * w_nh4sh.x[i][1][k];  // extrapolation
            w_nh4sh.x[i][jm-1][k] = w_nh4sh.x[i][jm-4][k] 
                - 3.0 * w_nh4sh.x[i][jm-3][k] + 3.0 * w_nh4sh.x[i][jm-2][k];  // extrapolation


            massflux_h2s.x[i][0][k] = massflux_h2s.x[i][3][k] 
                - 3.0 * massflux_h2s.x[i][2][k] + 3.0 * massflux_h2s.x[i][1][k];  // extrapolation
            massflux_h2s.x[i][jm-1][k] = massflux_h2s.x[i][jm-4][k] 
                - 3.0 * massflux_h2s.x[i][jm-3][k] + 3.0 * massflux_h2s.x[i][jm-2][k];  // extrapolation

            massflux_nh3.x[i][0][k] = massflux_nh3.x[i][3][k] 
                - 3.0 * massflux_nh3.x[i][2][k] + 3.0 * massflux_nh3.x[i][1][k];  // extrapolation
            massflux_nh3.x[i][jm-1][k] = massflux_nh3.x[i][jm-4][k] 
                - 3.0 * massflux_nh3.x[i][jm-3][k] + 3.0 * massflux_nh3.x[i][jm-2][k];  // extrapolation

            massflux_nh4sh.x[i][0][k] = massflux_nh4sh.x[i][3][k] 
                - 3.0 * massflux_nh4sh.x[i][2][k] + 3.0 * massflux_nh4sh.x[i][1][k];  // extrapolation
            massflux_nh4sh.x[i][jm-1][k] = massflux_nh4sh.x[i][jm-4][k] 
                - 3.0 * massflux_nh4sh.x[i][jm-3][k] + 3.0 * massflux_nh4sh.x[i][jm-2][k];  // extrapolation


            difflux_h2s.x[i][0][k] = difflux_h2s.x[i][3][k] 
                - 3.0 * difflux_h2s.x[i][2][k] + 3.0 * difflux_h2s.x[i][1][k];  // extrapolation
            difflux_h2s.x[i][jm-1][k] = difflux_h2s.x[i][jm-4][k] 
                - 3.0 * difflux_h2s.x[i][jm-3][k] + 3.0 * difflux_h2s.x[i][jm-2][k];  // extrapolation

            difflux_nh3.x[i][0][k] = difflux_nh3.x[i][3][k] 
                - 3.0 * difflux_nh3.x[i][2][k] + 3.0 * difflux_nh3.x[i][1][k];  // extrapolation
            difflux_nh3.x[i][jm-1][k] = difflux_nh3.x[i][jm-4][k] 
                - 3.0 * difflux_nh3.x[i][jm-3][k] + 3.0 * difflux_nh3.x[i][jm-2][k];  // extrapolation

            difflux_nh4sh.x[i][0][k] = difflux_nh4sh.x[i][3][k]
                - 3.0 * difflux_nh4sh.x[i][2][k] + 3.0 * difflux_nh4sh.x[i][1][k];  // extrapolation
            difflux_nh4sh.x[i][jm-1][k] = difflux_nh4sh.x[i][jm-4][k]
                - 3.0 * difflux_nh4sh.x[i][jm-3][k] + 3.0 * difflux_nh4sh.x[i][jm-2][k];  // extrapolation

            fluxlim_nh4sh.x[i][0][k] = 0.0;
            fluxlim_nh4sh.x[i][jm-1][k] = 0.0;


            thermalmassflux.x[i][0][k] = thermalmassflux.x[i][3][k] 
                - 3.0 * thermalmassflux.x[i][2][k] + 3.0 * thermalmassflux.x[i][1][k];  // extrapolation
            thermalmassflux.x[i][jm-1][k] = thermalmassflux.x[i][jm-4][k] 
                - 3.0 * thermalmassflux.x[i][jm-3][k] + 3.0 * thermalmassflux.x[i][jm-2][k];  // extrapolation


            CoriolisForce.x[i][0][k] = CoriolisForce.x[i][3][k] 
                - 3.0 * CoriolisForce.x[i][2][k] + 3.0 * CoriolisForce.x[i][1][k];  // extrapolation
            CoriolisForce.x[i][jm-1][k] = CoriolisForce.x[i][jm-4][k] 
                - 3.0 * CoriolisForce.x[i][jm-3][k] + 3.0 * CoriolisForce.x[i][jm-2][k];  // extrapolation

            CentrifugalForce.x[i][0][k] = CentrifugalForce.x[i][3][k] 
                - 3.0 * CentrifugalForce.x[i][2][k] + 3.0 * CentrifugalForce.x[i][1][k];  // extrapolation
            CentrifugalForce.x[i][jm-1][k] = CentrifugalForce.x[i][jm-4][k] 
                - 3.0 * CentrifugalForce.x[i][jm-3][k] + 3.0 * CentrifugalForce.x[i][jm-2][k];  // extrapolation

            BuoyancyForce.x[i][0][k] = BuoyancyForce.x[i][3][k] 
                - 3.0 * BuoyancyForce.x[i][2][k] + 3.0 * BuoyancyForce.x[i][1][k];  // extrapolation
            BuoyancyForce.x[i][jm-1][k] = BuoyancyForce.x[i][jm-4][k] 
                - 3.0 * BuoyancyForce.x[i][jm-3][k] + 3.0 * BuoyancyForce.x[i][jm-2][k];  // extrapolation

            PresGradForce.x[i][0][k] = PresGradForce.x[i][3][k] 
                - 3.0 * PresGradForce.x[i][2][k] + 3.0 * PresGradForce.x[i][1][k];  // extrapolation
            PresGradForce.x[i][jm-1][k] = PresGradForce.x[i][jm-4][k] 
                - 3.0 * PresGradForce.x[i][jm-3][k] + 3.0 * PresGradForce.x[i][jm-2][k];  // extrapolation

            Q_Latent.x[i][0][k] = Q_Latent.x[i][3][k] 
                - 3.0 * Q_Latent.x[i][2][k] + 3.0 * Q_Latent.x[i][1][k];  // extrapolation
            Q_Latent.x[i][jm-1][k] = Q_Latent.x[i][jm-4][k] 
                - 3.0 * Q_Latent.x[i][jm-3][k] + 3.0 * Q_Latent.x[i][jm-2][k];  // extrapolation

            Q_Sensible.x[i][0][k] = Q_Sensible.x[i][3][k] 
                - 3.0 * Q_Sensible.x[i][2][k] + 3.0 * Q_Sensible.x[i][1][k];  // extrapolation
            Q_Sensible.x[i][jm-1][k] = Q_Sensible.x[i][jm-4][k] 
                - 3.0 * Q_Sensible.x[i][jm-3][k] + 3.0 * Q_Sensible.x[i][jm-2][k];  // extrapolation


            // k*, dis*, nue* at the two poles (see the note at the top of this file).
            if(turb_bc){
                tke.x[i][0][k] = std::max(0.0,
                    c43 * tke.x[i][1][k] - c13 * tke.x[i][2][k]);
                dis.x[i][0][k] = std::max(bc_dis_min,
                    c43 * dis.x[i][1][k] - c13 * dis.x[i][2][k]);
                nue.x[i][0][k] = std::max(0.0,
                    c43 * nue.x[i][1][k] - c13 * nue.x[i][2][k]);

                tke.x[i][jm-1][k] = std::max(0.0,
                    c43 * tke.x[i][jm-2][k] - c13 * tke.x[i][jm-3][k]);
                dis.x[i][jm-1][k] = std::max(bc_dis_min,
                    c43 * dis.x[i][jm-2][k] - c13 * dis.x[i][jm-3][k]);
                nue.x[i][jm-1][k] = std::max(0.0,
                    c43 * nue.x[i][jm-2][k] - c13 * nue.x[i][jm-3][k]);
            }



/*
            h2o.x[i][0][k] = c43 * h2o.x[i][1][k] - c13 * h2o.x[i][2][k];
            h2o.x[i][jm-1][k] = c43 * h2o.x[i][jm-2][k] - c13 * h2o.x[i][jm-3][k];

            h2o_cloud.x[i][0][k] = c43 * h2o_cloud.x[i][1][k] - c13 * h2o_cloud.x[i][2][k];
            h2o_cloud.x[i][jm-1][k] = c43 * h2o_cloud.x[i][jm-2][k] - c13 * h2o_cloud.x[i][jm-3][k];

            h2o_ice.x[i][0][k] = c43 * h2o_ice.x[i][1][k] - c13 * h2o_ice.x[i][2][k];
            h2o_ice.x[i][jm-1][k] = c43 * h2o_ice.x[i][jm-2][k] - c13 * h2o_ice.x[i][jm-3][k];

            h2s.x[i][0][k] = c43 * h2s.x[i][1][k] - c13 * h2s.x[i][2][k];
            h2s.x[i][jm-1][k] = c43 * h2s.x[i][jm-2][k] - c13 * h2s.x[i][jm-3][k];

            nh3.x[i][0][k] = c43 * nh3.x[i][1][k] - c13 * nh3.x[i][2][k];
            nh3.x[i][jm-1][k] = c43 * nh3.x[i][jm-2][k] - c13 * nh3.x[i][jm-3][k];

            nh3_cloud.x[i][0][k] = c43 * nh3_cloud.x[i][1][k] - c13 * nh3_cloud.x[i][2][k];
            nh3_cloud.x[i][jm-1][k] = c43 * nh3_cloud.x[i][jm-2][k] - c13 * nh3_cloud.x[i][jm-3][k];

            nh3_ice.x[i][0][k] = c43 * nh3_ice.x[i][1][k] - c13 * nh3_ice.x[i][2][k];
            nh3_ice.x[i][jm-1][k] = c43 * nh3_ice.x[i][jm-2][k] - c13 * nh3_ice.x[i][jm-3][k];

            nh4sh.x[i][0][k] = c43 * nh4sh.x[i][1][k] - c13 * nh4sh.x[i][2][k];
            nh4sh.x[i][jm-1][k] = c43 * nh4sh.x[i][jm-2][k] - c13 * nh4sh.x[i][jm-3][k];


            j_h2s.x[i][0][k] = c43 * j_h2s.x[i][1][k] - c13 * j_h2s.x[i][2][k];
            j_h2s.x[i][jm-1][k] = c43 * j_h2s.x[i][jm-2][k] - c13 * j_h2s.x[i][jm-3][k];

            j_nh3.x[i][0][k] = c43 * j_nh3.x[i][1][k] - c13 * j_nh3.x[i][2][k];
            j_nh3.x[i][jm-1][k] = c43 * j_nh3.x[i][jm-2][k] - c13 * j_nh3.x[i][jm-3][k];

            j_nh4sh.x[i][0][k] = c43 * j_nh4sh.x[i][1][k] - c13 * j_nh4sh.x[i][2][k];
            j_nh4sh.x[i][jm-1][k] = c43 * j_nh4sh.x[i][jm-2][k] - c13 * j_nh4sh.x[i][jm-3][k];


            jT_h2s.x[i][0][k] = c43 * jT_h2s.x[i][1][k] - c13 * jT_h2s.x[i][2][k];
            jT_h2s.x[i][jm-1][k] = c43 * jT_h2s.x[i][jm-2][k] - c13 * jT_h2s.x[i][jm-3][k];

            jT_nh3.x[i][0][k] = c43 * jT_nh3.x[i][1][k] - c13 * jT_nh3.x[i][2][k];
            jT_nh3.x[i][jm-1][k] = c43 * jT_nh3.x[i][jm-2][k] - c13 * jT_nh3.x[i][jm-3][k];

            jT_nh4sh.x[i][0][k] = c43 * jT_nh4sh.x[i][1][k] - c13 * jT_nh4sh.x[i][2][k];
            jT_nh4sh.x[i][jm-1][k] = c43 * jT_nh4sh.x[i][jm-2][k] - c13 * jT_nh4sh.x[i][jm-3][k];


            w_h2s.x[i][0][k] = c43 * w_h2s.x[i][1][k] - c13 * w_h2s.x[i][2][k];
            w_h2s.x[i][jm-1][k] = c43 * w_h2s.x[i][jm-2][k] - c13 * w_h2s.x[i][jm-3][k];

            w_nh3.x[i][0][k] = c43 * w_nh3.x[i][1][k] - c13 * w_nh3.x[i][2][k];
            w_nh3.x[i][jm-1][k] = c43 * w_nh3.x[i][jm-2][k] - c13 * w_nh3.x[i][jm-3][k];

            w_nh4sh.x[i][0][k] = c43 * w_nh4sh.x[i][1][k] - c13 * w_nh4sh.x[i][2][k];
            w_nh4sh.x[i][jm-1][k] = c43 * w_nh4sh.x[i][jm-2][k] - c13 * w_nh4sh.x[i][jm-3][k];


            massflux_h2s.x[i][0][k] = c43 * massflux_h2s.x[i][1][k] - c13 * massflux_h2s.x[i][2][k];
            massflux_h2s.x[i][jm-1][k] = c43 * massflux_h2s.x[i][jm-2][k] - c13 * massflux_h2s.x[i][jm-3][k];

            massflux_nh3.x[i][0][k] = c43 * massflux_nh3.x[i][1][k] - c13 * massflux_nh3.x[i][2][k];
            massflux_nh3.x[i][jm-1][k] = c43 * massflux_nh3.x[i][jm-2][k] - c13 * massflux_nh3.x[i][jm-3][k];

            massflux_nh4sh.x[i][0][k] = c43 * massflux_nh4sh.x[i][1][k] - c13 * massflux_nh4sh.x[i][2][k];
            massflux_nh4sh.x[i][jm-1][k] = c43 * massflux_nh4sh.x[i][jm-2][k] - c13 * massflux_nh4sh.x[i][jm-3][k];

            difflux_h2s.x[i][0][k] = c43 * difflux_h2s.x[i][1][k] - c13 * difflux_h2s.x[i][2][k];
            difflux_h2s.x[i][jm-1][k] = c43 * difflux_h2s.x[i][jm-2][k] - c13 * difflux_h2s.x[i][jm-3][k];

            difflux_nh3.x[i][0][k] = c43 * difflux_nh3.x[i][1][k] - c13 * difflux_nh3.x[i][2][k];
            difflux_nh3.x[i][jm-1][k] = c43 * difflux_nh3.x[i][jm-2][k] - c13 * difflux_nh3.x[i][jm-3][k];

            difflux_nh4sh.x[i][0][k] = c43 * difflux_nh4sh.x[i][1][k] - c13 * difflux_nh4sh.x[i][2][k];
            difflux_nh4sh.x[i][jm-1][k] = c43 * difflux_nh4sh.x[i][jm-2][k] - c13 * difflux_nh4sh.x[i][jm-3][k];

            thermalmassflux.x[i][0][k] = c43 * thermalmassflux.x[i][1][k] - c13 * thermalmassflux.x[i][2][k];
            thermalmassflux.x[i][jm-1][k] = c43 * thermalmassflux.x[i][jm-2][k] - c13 * thermalmassflux.x[i][jm-3][k];

            CoriolisForce.x[i][0][k] = c43 * CoriolisForce.x[i][1][k] - c13 * CoriolisForce.x[i][2][k];
            CoriolisForce.x[i][jm-1][k] = c43 * CoriolisForce.x[i][jm-2][k] - c13 * CoriolisForce.x[i][jm-3][k];

            BuoyancyForce.x[i][0][k] = c43 * BuoyancyForce.x[i][1][k] - c13 * BuoyancyForce.x[i][2][k];
            BuoyancyForce.x[i][jm-1][k] = c43 * BuoyancyForce.x[i][jm-2][k] - c13 * BuoyancyForce.x[i][jm-3][k];

            PresGradForce.x[i][0][k] = c43 * PresGradForce.x[i][1][k] - c13 * PresGradForce.x[i][2][k];
            PresGradForce.x[i][jm-1][k] = c43 * PresGradForce.x[i][jm-2][k] - c13 * PresGradForce.x[i][jm-3][k];
*/
        }
    }

//    auto end = std::chrono::high_resolution_clock::now();
//    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
//    printf(" time measured: %.3f seconds for BC_theta\n", elapsed.count() * 1e-9)//;

//    cout << "      ATSAT: BC_theta ended" << endl;
    return;
}


void cSaturnModel::BC_phi(){
//    cout << endl << "      ATSAT: BC_phi" << endl;

//    auto begin = std::chrono::high_resolution_clock::now();

    const bool turb_bc = turb_active;

    #pragma omp parallel for
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            t.x[i][j][0] = c43 * t.x[i][j][1] - c13 * t.x[i][j][2];
            t.x[i][j][km-1] = c43 * t.x[i][j][km-2] - c13 * t.x[i][j][km-3];
            t.x[i][j][0] = t.x[i][j][km-1] = (t.x[i][j][0] + t.x[i][j][km-1])/2.0;

            u.x[i][j][0] = c43 * u.x[i][j][1] - c13 * u.x[i][j][2];
            u.x[i][j][km-1] = c43 * u.x[i][j][km-2] - c13 * u.x[i][j][km-3];
            u.x[i][j][0] = u.x[i][j][km-1] = (u.x[i][j][0] + u.x[i][j][km-1])/2.0;

            v.x[i][j][0] = c43 * v.x[i][j][1] - c13 * v.x[i][j][2];
            v.x[i][j][km-1] = c43 * v.x[i][j][km-2] - c13 * v.x[i][j][km-3];
            v.x[i][j][0] = v.x[i][j][km-1] = (v.x[i][j][0] + v.x[i][j][km-1])/2.0;

            w.x[i][j][0] = c43 * w.x[i][j][1] - c13 * w.x[i][j][2];
            w.x[i][j][km-1] = c43 * w.x[i][j][km-2] - c13 * w.x[i][j][km-3];
            w.x[i][j][0] = w.x[i][j][km-1] = (w.x[i][j][0] + w.x[i][j][km-1])/2.0;

            ch4.x[i][j][0] = c43 * ch4.x[i][j][1] - c13 * ch4.x[i][j][2];
            ch4.x[i][j][km-1] = c43 * ch4.x[i][j][km-2] - c13 * ch4.x[i][j][km-3];
            ch4.x[i][j][0] = ch4.x[i][j][km-1] = (ch4.x[i][j][0] + ch4.x[i][j][km-1])/2.0;

            ch4_cloud.x[i][j][0] = c43 * ch4_cloud.x[i][j][1] - c13 * ch4_cloud.x[i][j][2];
            ch4_cloud.x[i][j][km-1] = c43 * ch4_cloud.x[i][j][km-2] - c13 * ch4_cloud.x[i][j][km-3];
            ch4_cloud.x[i][j][0] = ch4_cloud.x[i][j][km-1] = (ch4_cloud.x[i][j][0] + ch4_cloud.x[i][j][km-1])/2.0;

            ch4_ice.x[i][j][0] = c43 * ch4_ice.x[i][j][1] - c13 * ch4_ice.x[i][j][2];
            ch4_ice.x[i][j][km-1] = c43 * ch4_ice.x[i][j][km-2] - c13 * ch4_ice.x[i][j][km-3];
            ch4_ice.x[i][j][0] = ch4_ice.x[i][j][km-1] = (ch4_ice.x[i][j][0] + ch4_ice.x[i][j][km-1])/2.0;

            h2o.x[i][j][0] = c43 * h2o.x[i][j][1] - c13 * h2o.x[i][j][2];
            h2o.x[i][j][km-1] = c43 * h2o.x[i][j][km-2] - c13 * h2o.x[i][j][km-3];
            h2o.x[i][j][0] = h2o.x[i][j][km-1] = (h2o.x[i][j][0] + h2o.x[i][j][km-1])/2.0;

            h2o_cloud.x[i][j][0] = c43 * h2o_cloud.x[i][j][1] - c13 * h2o_cloud.x[i][j][2];
            h2o_cloud.x[i][j][km-1] = c43 * h2o_cloud.x[i][j][km-2] - c13 * h2o_cloud.x[i][j][km-3];
            h2o_cloud.x[i][j][0] = h2o_cloud.x[i][j][km-1] = (h2o_cloud.x[i][j][0] + h2o_cloud.x[i][j][km-1])/2.0;

            h2o_ice.x[i][j][0] = c43 * h2o_ice.x[i][j][1] - c13 * h2o_ice.x[i][j][2];
            h2o_ice.x[i][j][km-1] = c43 * h2o_ice.x[i][j][km-2] - c13 * h2o_ice.x[i][j][km-3];
            h2o_ice.x[i][j][0] = h2o_ice.x[i][j][km-1] = (h2o_ice.x[i][j][0] + h2o_ice.x[i][j][km-1])/2.0;

            h2s.x[i][j][0] = c43 * h2s.x[i][j][1] - c13 * h2s.x[i][j][2];
            h2s.x[i][j][km-1] = c43 * h2s.x[i][j][km-2] - c13 * h2s.x[i][j][km-3];
            h2s.x[i][j][0] = h2s.x[i][j][km-1] = (h2s.x[i][j][0] + h2s.x[i][j][km-1])/2.0;

            nh3.x[i][j][0] = c43 * nh3.x[i][j][1] - c13 * nh3.x[i][j][2];
            nh3.x[i][j][km-1] = c43 * nh3.x[i][j][km-2] - c13 * nh3.x[i][j][km-3];
            nh3.x[i][j][0] = nh3.x[i][j][km-1] = (nh3.x[i][j][0] + nh3.x[i][j][km-1])/2.0;

            nh3_cloud.x[i][j][0] = c43 * nh3_cloud.x[i][j][1] - c13 * nh3_cloud.x[i][j][2];
            nh3_cloud.x[i][j][km-1] = c43 * nh3_cloud.x[i][j][km-2] - c13 * nh3_cloud.x[i][j][km-3];
            nh3_cloud.x[i][j][0] = nh3_cloud.x[i][j][km-1] = (nh3_cloud.x[i][j][0] + nh3_cloud.x[i][j][km-1])/2.0;

            nh3_ice.x[i][j][0] = c43 * nh3_ice.x[i][j][1] - c13 * nh3_ice.x[i][j][2];
            nh3_ice.x[i][j][km-1] = c43 * nh3_ice.x[i][j][km-2] - c13 * nh3_ice.x[i][j][km-3];
            nh3_ice.x[i][j][0] = nh3_ice.x[i][j][km-1] = (nh3_ice.x[i][j][0] + nh3_ice.x[i][j][km-1])/2.0;

            nh4sh.x[i][j][0] = c43 * nh4sh.x[i][j][1] - c13 * nh4sh.x[i][j][2];
            nh4sh.x[i][j][km-1] = c43 * nh4sh.x[i][j][km-2] - c13 * nh4sh.x[i][j][km-3];
            nh4sh.x[i][j][0] = nh4sh.x[i][j][km-1] = (nh4sh.x[i][j][0] + nh4sh.x[i][j][km-1])/2.0;

            j_h2s.x[i][j][0] = c43 * j_h2s.x[i][j][1] - c13 * j_h2s.x[i][j][2];
            j_h2s.x[i][j][km-1] = c43 * j_h2s.x[i][j][km-2] - c13 * j_h2s.x[i][j][km-3];
            j_h2s.x[i][j][0] = j_h2s.x[i][j][km-1] = (j_h2s.x[i][j][0] + j_h2s.x[i][j][km-1])/2.0;

            j_nh3.x[i][j][0] = c43 * j_nh3.x[i][j][1] - c13 * j_nh3.x[i][j][2];
            j_nh3.x[i][j][km-1] = c43 * j_nh3.x[i][j][km-2] - c13 * j_nh3.x[i][j][km-3];
            j_nh3.x[i][j][0] = j_nh3.x[i][j][km-1] = (j_nh3.x[i][j][0] + j_nh3.x[i][j][km-1])/2.0;

            j_nh4sh.x[i][j][0] = c43 * j_nh4sh.x[i][j][1] - c13 * j_nh4sh.x[i][j][2];
            j_nh4sh.x[i][j][km-1] = c43 * j_nh4sh.x[i][j][km-2] - c13 * j_nh4sh.x[i][j][km-3];
            j_nh4sh.x[i][j][0] = j_nh4sh.x[i][j][km-1] = (j_nh4sh.x[i][j][0] + j_nh4sh.x[i][j][km-1])/2.0;

            jT_h2s.x[i][j][0] = c43 * jT_h2s.x[i][j][1] - c13 * jT_h2s.x[i][j][2];
            jT_h2s.x[i][j][km-1] = c43 * jT_h2s.x[i][j][km-2] - c13 * jT_h2s.x[i][j][km-3];
            jT_h2s.x[i][j][0] = jT_h2s.x[i][j][km-1] = (jT_h2s.x[i][j][0] + jT_h2s.x[i][j][km-1])/2.0;

            jT_nh3.x[i][j][0] = c43 * jT_nh3.x[i][j][1] - c13 * jT_nh3.x[i][j][2];
            jT_nh3.x[i][j][km-1] = c43 * jT_nh3.x[i][j][km-2] - c13 * jT_nh3.x[i][j][km-3];
            jT_nh3.x[i][j][0] = jT_nh3.x[i][j][km-1] = (jT_nh3.x[i][j][0] + jT_nh3.x[i][j][km-1])/2.0;

            jT_nh4sh.x[i][j][0] = c43 * jT_nh4sh.x[i][j][1] - c13 * jT_nh4sh.x[i][j][2];
            jT_nh4sh.x[i][j][km-1] = c43 * jT_nh4sh.x[i][j][km-2] - c13 * jT_nh4sh.x[i][j][km-3];
            jT_nh4sh.x[i][j][0] = jT_nh4sh.x[i][j][km-1] = (jT_nh4sh.x[i][j][0] + jT_nh4sh.x[i][j][km-1])/2.0;

            w_h2s.x[i][j][0] = c43 * w_h2s.x[i][j][1] - c13 * w_h2s.x[i][j][2];
            w_h2s.x[i][j][km-1] = c43 * w_h2s.x[i][j][km-2] - c13 * w_h2s.x[i][j][km-3];
            w_h2s.x[i][j][0] = w_h2s.x[i][j][km-1] = (w_h2s.x[i][j][0] + w_h2s.x[i][j][km-1])/2.0;

            w_nh3.x[i][j][0] = c43 * w_nh3.x[i][j][1] - c13 * w_nh3.x[i][j][2];
            w_nh3.x[i][j][km-1] = c43 * w_nh3.x[i][j][km-2] - c13 * w_nh3.x[i][j][km-3];
            w_nh3.x[i][j][0] = w_nh3.x[i][j][km-1] = (w_nh3.x[i][j][0] + w_nh3.x[i][j][km-1])/2.0;

            w_nh4sh.x[i][j][0] = c43 * w_nh4sh.x[i][j][1] - c13 * w_nh4sh.x[i][j][2];
            w_nh4sh.x[i][j][km-1] = c43 * w_nh4sh.x[i][j][km-2] - c13 * w_nh4sh.x[i][j][km-3];
            w_nh4sh.x[i][j][0] = w_nh4sh.x[i][j][km-1] = (w_nh4sh.x[i][j][0] + w_nh4sh.x[i][j][km-1])/2.0;

            massflux_h2s.x[i][j][0] = c43 * massflux_h2s.x[i][j][1] - c13 * massflux_h2s.x[i][j][2];
            massflux_h2s.x[i][j][km-1] = c43 * massflux_h2s.x[i][j][km-2] - c13 * massflux_h2s.x[i][j][km-3];
            massflux_h2s.x[i][j][0] = massflux_h2s.x[i][j][km-1] = (massflux_h2s.x[i][j][0] + massflux_h2s.x[i][j][km-1])/2.0;

            massflux_nh3.x[i][j][0] = c43 * massflux_nh3.x[i][j][1] - c13 * massflux_nh3.x[i][j][2];
            massflux_nh3.x[i][j][km-1] = c43 * massflux_nh3.x[i][j][km-2] - c13 * massflux_nh3.x[i][j][km-3];
            massflux_nh3.x[i][j][0] = massflux_nh3.x[i][j][km-1] = (massflux_nh3.x[i][j][0] + massflux_nh3.x[i][j][km-1])/2.0;

            massflux_nh4sh.x[i][j][0] = c43 * massflux_nh4sh.x[i][j][1] - c13 * massflux_nh4sh.x[i][j][2];
            massflux_nh4sh.x[i][j][km-1] = c43 * massflux_nh4sh.x[i][j][km-2] - c13 * massflux_nh4sh.x[i][j][km-3];
            massflux_nh4sh.x[i][j][0] = massflux_nh4sh.x[i][j][km-1] = (massflux_nh4sh.x[i][j][0] + massflux_nh4sh.x[i][j][km-1])/2.0;

            difflux_h2s.x[i][j][0] = c43 * difflux_h2s.x[i][j][1] - c13 * difflux_h2s.x[i][j][2];
            difflux_h2s.x[i][j][km-1] = c43 * difflux_h2s.x[i][j][km-2] - c13 * difflux_h2s.x[i][j][km-3];
            difflux_h2s.x[i][j][0] = difflux_h2s.x[i][j][km-1] = (difflux_h2s.x[i][j][0] + difflux_h2s.x[i][j][km-1])/2.0;

            difflux_nh3.x[i][j][0] = c43 * difflux_nh3.x[i][j][1] - c13 * difflux_nh3.x[i][j][2];
            difflux_nh3.x[i][j][km-1] = c43 * difflux_nh3.x[i][j][km-2] - c13 * difflux_nh3.x[i][j][km-3];
            difflux_nh3.x[i][j][0] = difflux_nh3.x[i][j][km-1] = (difflux_nh3.x[i][j][0] + difflux_nh3.x[i][j][km-1])/2.0;

            difflux_nh4sh.x[i][j][0] = c43 * difflux_nh4sh.x[i][j][1] - c13 * difflux_nh4sh.x[i][j][2];
            difflux_nh4sh.x[i][j][km-1] = c43 * difflux_nh4sh.x[i][j][km-2] - c13 * difflux_nh4sh.x[i][j][km-3];
            difflux_nh4sh.x[i][j][0] = difflux_nh4sh.x[i][j][km-1] = (difflux_nh4sh.x[i][j][0] + difflux_nh4sh.x[i][j][km-1])/2.0;

            fluxlim_nh4sh.x[i][j][0] = c43 * fluxlim_nh4sh.x[i][j][1] - c13 * fluxlim_nh4sh.x[i][j][2];
            fluxlim_nh4sh.x[i][j][km-1] = c43 * fluxlim_nh4sh.x[i][j][km-2] - c13 * fluxlim_nh4sh.x[i][j][km-3];
            fluxlim_nh4sh.x[i][j][0] = fluxlim_nh4sh.x[i][j][km-1] = (fluxlim_nh4sh.x[i][j][0] + fluxlim_nh4sh.x[i][j][km-1])/2.0;

            thermalmassflux.x[i][j][0] = c43 * thermalmassflux.x[i][j][1] - c13 * thermalmassflux.x[i][j][2];
            thermalmassflux.x[i][j][km-1] = c43 * thermalmassflux.x[i][j][km-2] - c13 * thermalmassflux.x[i][j][km-3];
            thermalmassflux.x[i][j][0] = thermalmassflux.x[i][j][km-1] = (thermalmassflux.x[i][j][0] + thermalmassflux.x[i][j][km-1])/2.0;

            CoriolisForce.x[i][j][0] = c43 * CoriolisForce.x[i][j][1] - c13 * CoriolisForce.x[i][j][2];
            CoriolisForce.x[i][j][km-1] = c43 * CoriolisForce.x[i][j][km-2] - c13 * CoriolisForce.x[i][j][km-3];
            CoriolisForce.x[i][j][0] = CoriolisForce.x[i][j][km-1] = (CoriolisForce.x[i][j][0] + CoriolisForce.x[i][j][km-1])/2.0;

            CentrifugalForce.x[i][j][0] = c43 * CentrifugalForce.x[i][j][1] - c13 * CentrifugalForce.x[i][j][2];
            CentrifugalForce.x[i][j][km-1] = c43 * CentrifugalForce.x[i][j][km-2] - c13 * CentrifugalForce.x[i][j][km-3];
            CentrifugalForce.x[i][j][0] = CentrifugalForce.x[i][j][km-1] = (CentrifugalForce.x[i][j][0] + CentrifugalForce.x[i][j][km-1])/2.0;

            PresGradForce.x[i][j][0] = c43 * PresGradForce.x[i][j][1] - c13 * PresGradForce.x[i][j][2];
            PresGradForce.x[i][j][km-1] = c43 * PresGradForce.x[i][j][km-2] - c13 * PresGradForce.x[i][j][km-3];
            PresGradForce.x[i][j][0] = PresGradForce.x[i][j][km-1] = (PresGradForce.x[i][j][0] + PresGradForce.x[i][j][km-1])/2.0;

            BuoyancyForce.x[i][j][0] = c43 * BuoyancyForce.x[i][j][1] - c13 * BuoyancyForce.x[i][j][2];
            BuoyancyForce.x[i][j][km-1] = c43 * BuoyancyForce.x[i][j][km-2] - c13 * BuoyancyForce.x[i][j][km-3];
            BuoyancyForce.x[i][j][0] = BuoyancyForce.x[i][j][km-1] = (BuoyancyForce.x[i][j][0] + BuoyancyForce.x[i][j][km-1])/2.0;

            Q_Latent.x[i][j][0] = c43 * Q_Latent.x[i][j][1] - c13 * Q_Latent.x[i][j][2];
            Q_Latent.x[i][j][km-1] = c43 * Q_Latent.x[i][j][km-2] - c13 * Q_Latent.x[i][j][km-3];
            Q_Latent.x[i][j][0] = Q_Latent.x[i][j][km-1] = (Q_Latent.x[i][j][0] + Q_Latent.x[i][j][km-1])/2.0;

            Q_Sensible.x[i][j][0] = c43 * Q_Sensible.x[i][j][1] - c13 * Q_Sensible.x[i][j][2];
            Q_Sensible.x[i][j][km-1] = c43 * Q_Sensible.x[i][j][km-2] - c13 * Q_Sensible.x[i][j][km-3];
            Q_Sensible.x[i][j][0] = Q_Sensible.x[i][j][km-1] = (Q_Sensible.x[i][j][0] + Q_Sensible.x[i][j][km-1])/2.0;

            // k*, dis*, nue* across the phi seam (see the note at the top of this file). The two
            // faces are averaged and set equal, as every other field here is: 0 and km-1 are the
            // same meridian, so a jump between them is a discontinuity in the middle of the
            // domain, not a boundary.
            if(turb_bc){
                tke.x[i][j][0] = c43 * tke.x[i][j][1] - c13 * tke.x[i][j][2];
                tke.x[i][j][km-1] = c43 * tke.x[i][j][km-2] - c13 * tke.x[i][j][km-3];
                tke.x[i][j][0] = tke.x[i][j][km-1] =
                    std::max(0.0, (tke.x[i][j][0] + tke.x[i][j][km-1])/2.0);

                dis.x[i][j][0] = c43 * dis.x[i][j][1] - c13 * dis.x[i][j][2];
                dis.x[i][j][km-1] = c43 * dis.x[i][j][km-2] - c13 * dis.x[i][j][km-3];
                dis.x[i][j][0] = dis.x[i][j][km-1] =
                    std::max(bc_dis_min, (dis.x[i][j][0] + dis.x[i][j][km-1])/2.0);

                nue.x[i][j][0] = c43 * nue.x[i][j][1] - c13 * nue.x[i][j][2];
                nue.x[i][j][km-1] = c43 * nue.x[i][j][km-2] - c13 * nue.x[i][j][km-3];
                nue.x[i][j][0] = nue.x[i][j][km-1] =
                    std::max(0.0, (nue.x[i][j][0] + nue.x[i][j][km-1])/2.0);
            }
        }
    }

//    auto end = std::chrono::high_resolution_clock::now();
//    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
//    printf(" time measured: %.3f seconds for BC_phi\n", elapsed.count() * 1e-9);

//    cout << "      ATSAT: BC_phi ended" << endl;
    return;
}
/*
*
*/
/*
void cSaturnModel::BC_seamount(){
    cout << endl << "      ATSAT: BC_seamount" << endl;
    // seamount of elliptical cross section with height i_0 and half axes a and b

//    auto begin = std::chrono::high_resolution_clock::now();

    int j_ellipse = 0;
    int a_0 = 1;
    int b_0 = 1;
    int a = im-1;
    int b = im-1;
    int j_0 = 112; // center of Great Red Spot 22° south of equator
    int k_0 = 180;
    int i_0 = 37;
//    int i_0 = 30;

    for(int i = 0; i < i_0; i++){
        for(int k = k_0-a; k <= k_0+a; k++){  // elliptic contour
            for(int j = j_0-b; j <= j_0+b; j++){
                if(j <= j_0){
                    j_ellipse = j_0 - sqrt(pow(b,2) - pow((b * (k-k_0)/a),2));
                    if(j >= j_ellipse)  SeaMount.x[i][j][k] = 1.0;
                }else{
                    j_ellipse = j_0 + sqrt(pow(b,2) - pow((b * (k-k_0)/a),2));
                    if(j <= j_ellipse)  SeaMount.x[i][j][k] = 1.0;
                }
//if(i==37)cout << "   " << i << "   " << j << "   " << k << "   " << j_ellipse << "   " << SeaMount.x[i][j_ellipse][k] << "   " << SeaMount.x[i][j][k] << endl;
                if(i == i_0) Topography.y[j][k] = SeaMount.x[i_0][j][k]; // topography
            }
        }
        a = a_0 * ((im-1)-i);
        b = b_0 * ((im-1)-i);
    }


    for(int k = k_0-a; k <= k_0+a; k++){  // square contour
        for(int j = j_0-b; j <= j_0+b; j++){
            for(int i = 0; i < i_0; i++){
                SeaMount.x[i][j][k] = 1.0;
            }
        }
    }

//    auto end = std::chrono::high_resolution_clock::now();
//    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
//    printf(" time measured: %.3f seconds for BC_seamount\n", elapsed.count() * 1e-9);

    cout << "      ATSAT: BC_seamount ended" << endl;
    return;
}
*/
/*
*
*/
/*
void cSaturnModel::BC_solidground(){
    cout << endl << "      ATSAT: BC_solidground" << endl;

//    auto begin = std::chrono::high_resolution_clock::now();

//    #pragma omp parallel for
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            for(int i = 0; i < im; i++){
                if(SeaMount.x[i][j][k] == 1.0){
                    u.x[i][j][k] = 0.0;
                    v.x[i][j][k] = 0.0;
                    w.x[i][j][k] = 0.0;

                    t.x[i][j][k] = 1.0;
                    p_stat.x[i][j][k] = 1.0;
                    p_dyn.x[i][j][k] = 0.0;

                    h2o.x[i][j][k] = 0.0;
                    h2o_cloud.x[i][j][k] = 0.0;
                    h2o_ice.x[i][j][k] = 0.0;

                    h2s.x[i][j][k] = 0.0;
                    h2s_cloud.x[i][j][k] = 0.0;
                    h2s_ice.x[i][j][k] = 0.0;
                    j_h2s.x[i][j][k] = 0.0;
                    jT_h2s.x[i][j][k] = 0.0;

                    nh3.x[i][j][k] = 0.0;
                    nh3_cloud.x[i][j][k] = 0.0;
                    nh3.x[i][j][k] = 0.0;
                    j_nh3.x[i][j][k] = 0.0;
                    jT_nh3.x[i][j][k] = 0.0;

                    nh4sh.x[i][j][k] = 0.0;
//                    nh4sh_cloud.x[i][j][k] = 0.0;
//                    nh4sh_ice.x[i][j][k] = 0.0;
                    j_nh4sh.x[i][j][k] = 0.0;
                    jT_nh4sh.x[i][j][k] = 0.0;

                    massflux_h2s.x[i][j][k] = 0.0;
                    massflux_nh3.x[i][j][k] = 0.0;
                    massflux_nh4sh.x[i][j][k] = 0.0;
                    thermalmassflux.x[i][j][k] = 0.0;

                    w_h2s.x[i][j][k] = 0.0;
                    w_nh3.x[i][j][k] = 0.0;
                    w_nh4sh.x[i][j][k] = 0.0;
//                    w_nh4sh_tot.x[i][j][k] = 0.0;

                    difflux_h2s.x[i][j][k] = 0.0;
                    difflux_nh3.x[i][j][k] = 0.0;
                    difflux_nh4sh.x[i][j][k] = 0.0;

                    fluxlim_nh4sh.x[i][j][k] = 0.0;

                }

                if(h2s.x[i][j][k] < 0.0)  h2s.x[i][j][k] = 0.0;
            }
        }
    }


//    auto end = std::chrono::high_resolution_clock::now();
//    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
//    printf(" time measured: %.3f seconds for BC_solidground\n", elapsed.count() * 1e-9);

    cout << "      ATSAT: BC_solidground ended" << endl;
    return;
}
*/
/*
*
*/
void cSaturnModel::init_tropopause_layers(){
    tropopause_layers = std::vector<double>(jm, tropopause_pole);
    cout << endl << "      ATSAT: init_tropopause_layers" << endl;
// Versiera di Agnesi approach, two inflection points
    int i_max = im-1;
    int j_max = jm-1;
    int j_half = j_max/2;
//    double coeff_pole = 280.0;
    double coeff_pole = 285.0;

    for(int j=j_half; j>=0; j--){

        double x = coeff_pole * (1.0 - (double)(j_half-j)/(double)j_half);
        tropopause_layers[j] = AtomUtils::Agnesi(tropopause_equator, x); 
        tropopause_layers[j] = round(tropopause_layers[j] 
           /L_atm * (double)i_max); 

        tropopause_layers[j] = tropopause_equator/L_atm * (double)i_max;

/*
    cout << "   j = " << j << "   x = " << x 
        << "   Agnesi = " << AtomUtils::Agnesi(tropopause_equator, x) 
        << "   tropopause_equator = " << tropopause_equator 
        << "   tropopause_pole = " << tropopause_pole
        << "   tropopause_layers = " << round(tropopause_layers[j]) << endl;
*/

    }

    for(int j=j_max; j>j_half; j--){
        tropopause_layers[j] = tropopause_layers[j_max-j];
    }


    cout << "      ATSAT: init_tropopause_layers ended" << endl;
    return;
}
/*
*
*/
