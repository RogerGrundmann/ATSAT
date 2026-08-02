#include "cSaturnModel.h"
#include "BC_Sat.h"
#include "ATPhys.h"   // env_int, for the boundary-hardening knobs below

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
//
//
// ===== Why the three routines below are written as field lists =====
// Each of BC_radius, BC_theta and BC_phi applies ONE formula to every transported field. Written
// longhand — a named block per field per boundary — that fact was invisible: the file was 1121
// lines in which the same two lines of arithmetic appeared some 110 times with a different
// identifier pasted in, and the only way to see that (say) fluxlim_nh4sh is extrapolated at the
// radial walls but pinned to zero at the poles was to read all of it.
//
// The lists below say it directly. A field participates in a boundary if and only if it appears
// in that boundary's list, and changing the FORM of a boundary is now an edit to one expression
// rather than to ~38 blocks. That is the point: ATJUP already carries BC_TOP_TAPER, BC_POLE_COPY,
// BC_RADIUS_COPY and BC_T_LID_PIN as live branches on exactly these expressions (BC_Jup.h), and
// they could not be ported into the longhand form without tripling it. This commit changes no
// arithmetic; it only makes the next one possible.
//
// DELETED WITH THE LONGHAND: roughly 300 lines of commented-out alternative forms, all of them
// the 2-point Neumann f[s] = (4/3)f[a] - (1/3)f[b] that these routines used before the cubic.
// They are not lost — that form is what ATJUP_BC_RADIUS_COPY and friends select, so it returns
// as a live, testable branch rather than as a comment. Two of those dead blocks had typos that
// would have become bugs the moment anyone uncommented them, which is its own argument against
// keeping code in comments: BC_radius wrote nh3_cloud[0] from nh3[2] rather than nh3_cloud[2],
// and BC_theta wrote v[i][0][k] from v.x[3][3][k] rather than v.x[i][3][k].
namespace {
    constexpr double bc_dis_min = 1.0e-10;   // matches TurbulenceSat::dis_min and the RK4 floor

    // ===== Boundary-hardening knobs, ported from ATJUP's BCJupKnobs (BC_Jup.h) =====
    // ATSAT's convention, unlike ATJUP's: every one of these defaults to OFF, so the default
    // path is byte-identical to the code before the port and each knob is an independently
    // answerable question rather than a change smuggled in with its neighbours. ATJUP ships
    // top_taper and pole_copy ON, on Jupiter evidence; whether that evidence transfers to a
    // 500 km Saturn shell is exactly what these are here to measure.
    inline int knob(const char* name, int dflt){
        return ATPhys::env_int(cSaturnModel::planet_tag(), name, dflt);
    }

    // (2) Taper v,w to zero over the top three layers. See the note at the call site.
    inline int top_taper(){ static const int v = knob("BC_TOP_TAPER", 0); return v; }
}

void BC_Sat::bcRadius() { m.BC_radius(); }
void BC_Sat::bcTheta()  { m.BC_theta(); }
void BC_Sat::bcPhi()    { m.BC_phi(); }


// ---------------------------------------------------------------------------
// Deep boundary i=0 and model top i=im-1.
// ---------------------------------------------------------------------------
void cSaturnModel::BC_radius(){
    // 3-point cubic extrapolation, inward from each radial wall.
    Array* fields[] = {
        &t, &u, &v, &w,
        &ch4, &ch4_cloud, &ch4_ice,
        &h2o, &h2o_cloud, &h2o_ice,
        &h2s,
        &nh3, &nh3_cloud, &nh3_ice,
        &nh4sh,
        &j_h2s,  &j_nh3,  &j_nh4sh,
        &jT_h2s, &jT_nh3, &jT_nh4sh,
        &w_h2s,  &w_nh3,  &w_nh4sh,
        &massflux_h2s, &massflux_nh3, &massflux_nh4sh,
        &difflux_h2s,  &difflux_nh3,  &difflux_nh4sh,
        &fluxlim_nh4sh,
        &thermalmassflux,
        &CoriolisForce, &CentrifugalForce, &BuoyancyForce, &PresGradForce,
        &Q_Latent, &Q_Sensible
    };
    const int nf = (int)(sizeof(fields) / sizeof(fields[0]));

    // k*, dis*, nue* at the deep boundary and the model top: 2-point Neumann and a floor, for
    // the reasons in the note at the top of this file.
    Array* turb[] = { &tke, &dis, &nue };
    const double turb_floor[] = { 0.0, bc_dis_min, 0.0 };
    const int nt = (int)(sizeof(turb) / sizeof(turb[0]));
    const bool turb_bc = turb_active;

  #pragma omp parallel for
    for(int j = 1; j < jm-1; j++){
        for(int k = 1; k < km-1; k++){
            for(int f = 0; f < nf; f++){
                Array& F = *fields[f];
                F.x[0][j][k] = F.x[3][j][k]
                    - 3.0 * F.x[2][j][k] + 3.0 * F.x[1][j][k];  // extrapolation
                F.x[im-1][j][k] = F.x[im-4][j][k]
                    - 3.0 * F.x[im-3][j][k] + 3.0 * F.x[im-2][j][k];  // extrapolation
            }

            if(turb_bc){
                for(int f = 0; f < nt; f++){
                    Array& F = *turb[f];
                    F.x[0][j][k] = std::max(turb_floor[f],
                        c43 * F.x[1][j][k] - c13 * F.x[2][j][k]);
                    F.x[im-1][j][k] = std::max(turb_floor[f],
                        c43 * F.x[im-2][j][k] - c13 * F.x[im-3][j][k]);
                }
            }
        }
    }

    // --- (2) Taper the HORIZONTAL velocities to a quiet grid ceiling. ATSAT_BC_TOP_TAPER=1. ---
    //
    // init_v_or_w_above_tropopause() (InitVelocity_Sat.cpp) ramps v and w linearly to zero
    // between the tropopause and the model top, so the initial condition has v = w = 0 at
    // i = im-1 by construction. The extrapolation above then copies the interior value straight
    // back onto the lid and undoes it, dragging the zonal jet up to the ceiling. Ramping the top
    // three layers by 2/3, 1/3, 0 restores a quiet lid without the one-cell shear shock that a
    // hard zero at im-1 alone would create.
    //
    // TWO FORMS, because ATJUP's rationale for its own does not survive measurement:
    //
    //   =1  ATJUP parity. Multiply in place, exactly as BC_Jup.h does.
    //   =2  Re-derive from the first untapered level each call, so nothing compounds.
    //
    // ATJUP's comment claims the multiply cannot run away, "because RK4 integrates i = 1..im-2,
    // so v,w at im-2 and im-3 are recomputed from tendencies every iteration and the factor is
    // re-applied to a fresh value rather than to an already-tapered one". That is not what an
    // incremental integrator does. RK4 forms v = vn + dt*rhs, restoreVar then copies v back into
    // vn, so the value the factor multiplies at iteration n+1 is the value it already multiplied
    // at iteration n. The factor compounds geometrically and only the tendency replenishes it.
    //
    // MEASURED on ATSAT, 30 iterations, max|w| by level (=0 against =1):
    //     i=38  0.5717 -> 0.000056        (2/3)^30 = 5.2e-6
    //     i=39  0.4362 -> 0.0             (1/3)^30 = 5.2e-15
    //     i=40  0.1505 -> 0.0
    // So =1 is not a 2/3, 1/3, 0 ramp at all after the first few iterations: it is a hard zero
    // three layers deep, reached geometrically. It does do the job it was ported for — the |w|
    // maximum moves off the lid, i=38 (475 km) to i=37 (462 km) — but by emptying the top of the
    // shell rather than by grading it.
    //
    // =2 gets the intended profile and keeps it: v,w at im-2 and im-3 are set to 1/3 and 2/3 of
    // the first level the taper does not touch, i = im-4, which RK4 refreshes every iteration.
    // That is a genuine linear ramp from the interior to zero at the lid, re-derived from a fresh
    // value each call, which is what the ATJUP comment describes and its arithmetic does not do.
    //
    // i = im-1 is outside the RK4 range in both forms, so zeroing it is a clean Dirichlet
    // condition either way.
    //
    // WHAT DOES NOT TRANSFER FROM ATJUP. The three tapered layers are 3/40 of the shell in both
    // models, but ATSAT's shell is 500 km against ATJUP's 140, so this reaches ~37 km down
    // instead of ~10 km. It is also a far smaller fraction of the ramp it is meant to protect:
    // ATSAT's tropopause sits at level 10 of 40 (tropopause_equator 125 km of L_atm 500), so the
    // IC ramp already spans 30 levels and the taper touches only its top tenth.
    //
    // Its own loop, over the FULL j,k range, rather than inside the loop above: that loop runs
    // j = 1..jm-2, k = 1..km-2, and leaving the seam columns untapered would put a discontinuity
    // at k = 0/km-1 in the middle of the domain. bcRadius() is called last of the three
    // (cSaturnModel.cpp), so nothing overwrites this afterwards.
    const int taper = top_taper();
    if(taper != 0){
        const int iml = im - 1;
      #pragma omp parallel for
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                v.x[iml][j][k] = 0.0;
                w.x[iml][j][k] = 0.0;
                if(taper == 2){
                    v.x[iml-1][j][k] = (1.0/3.0) * v.x[iml-3][j][k];
                    w.x[iml-1][j][k] = (1.0/3.0) * w.x[iml-3][j][k];
                    v.x[iml-2][j][k] = (2.0/3.0) * v.x[iml-3][j][k];
                    w.x[iml-2][j][k] = (2.0/3.0) * w.x[iml-3][j][k];
                }else{
                    v.x[iml-1][j][k] *= (1.0/3.0);
                    w.x[iml-1][j][k] *= (1.0/3.0);
                    v.x[iml-2][j][k] *= (2.0/3.0);
                    w.x[iml-2][j][k] *= (2.0/3.0);
                }
            }
        }
    }

    return;
}
/*
*
*/
// ---------------------------------------------------------------------------
// The two poles, j=0 and j=jm-1.
// ---------------------------------------------------------------------------
void cSaturnModel::BC_theta(){
    // v, w and the nh4sh flux limiter are pinned to zero at both poles rather than extrapolated:
    // the meridional and zonal velocities have no meaning on the axis, and the limiter is a
    // correction to an advective flux that does not exist there.
    Array* zero_fields[] = { &v, &w, &fluxlim_nh4sh };
    const int nz = (int)(sizeof(zero_fields) / sizeof(zero_fields[0]));

    // Everything else takes the same 3-point cubic used at the radial walls. Note that u IS
    // extrapolated here and v, w are not — the radial component is tangential to the pole.
    Array* fields[] = {
        &t, &u,
        &ch4, &ch4_cloud, &ch4_ice,
        &h2o, &h2o_cloud, &h2o_ice,
        &h2s,
        &nh3, &nh3_cloud, &nh3_ice,
        &nh4sh,
        &j_h2s,  &j_nh3,  &j_nh4sh,
        &jT_h2s, &jT_nh3, &jT_nh4sh,
        &w_h2s,  &w_nh3,  &w_nh4sh,
        &massflux_h2s, &massflux_nh3, &massflux_nh4sh,
        &difflux_h2s,  &difflux_nh3,  &difflux_nh4sh,
        &thermalmassflux,
        &CoriolisForce, &CentrifugalForce, &BuoyancyForce, &PresGradForce,
        &Q_Latent, &Q_Sensible
    };
    const int nf = (int)(sizeof(fields) / sizeof(fields[0]));

    // k*, dis*, nue* at the two poles (see the note at the top of this file).
    Array* turb[] = { &tke, &dis, &nue };
    const double turb_floor[] = { 0.0, bc_dis_min, 0.0 };
    const int nt = (int)(sizeof(turb) / sizeof(turb[0]));
    const bool turb_bc = turb_active;

    #pragma omp parallel for
    for(int k = 1; k < km-1; k++){
        for(int i = 1; i < im-1; i++){
            for(int f = 0; f < nz; f++){
                Array& F = *zero_fields[f];
                F.x[i][0][k] = 0.0;
                F.x[i][jm-1][k] = 0.0;
            }

            for(int f = 0; f < nf; f++){
                Array& F = *fields[f];
                F.x[i][0][k] = F.x[i][3][k]
                    - 3.0 * F.x[i][2][k] + 3.0 * F.x[i][1][k];  // extrapolation
                F.x[i][jm-1][k] = F.x[i][jm-4][k]
                    - 3.0 * F.x[i][jm-3][k] + 3.0 * F.x[i][jm-2][k];  // extrapolation
            }

            if(turb_bc){
                for(int f = 0; f < nt; f++){
                    Array& F = *turb[f];
                    F.x[i][0][k] = std::max(turb_floor[f],
                        c43 * F.x[i][1][k] - c13 * F.x[i][2][k]);
                    F.x[i][jm-1][k] = std::max(turb_floor[f],
                        c43 * F.x[i][jm-2][k] - c13 * F.x[i][jm-3][k]);
                }
            }
        }
    }

    return;
}
/*
*
*/
// ---------------------------------------------------------------------------
// The phi seam, k=0 and k=km-1.
// ---------------------------------------------------------------------------
void cSaturnModel::BC_phi(){
    // k=0 and k=km-1 are the SAME meridian, so this is not a boundary condition in the sense of
    // the two above: each face is extrapolated from its own side and the two are then averaged
    // and set equal, because a jump between them is a discontinuity in the middle of the domain.
    Array* fields[] = {
        &t, &u, &v, &w,
        &ch4, &ch4_cloud, &ch4_ice,
        &h2o, &h2o_cloud, &h2o_ice,
        &h2s,
        &nh3, &nh3_cloud, &nh3_ice,
        &nh4sh,
        &j_h2s,  &j_nh3,  &j_nh4sh,
        &jT_h2s, &jT_nh3, &jT_nh4sh,
        &w_h2s,  &w_nh3,  &w_nh4sh,
        &massflux_h2s, &massflux_nh3, &massflux_nh4sh,
        &difflux_h2s,  &difflux_nh3,  &difflux_nh4sh,
        &fluxlim_nh4sh,
        &thermalmassflux,
        &CoriolisForce, &CentrifugalForce, &PresGradForce, &BuoyancyForce,
        &Q_Latent, &Q_Sensible
    };
    const int nf = (int)(sizeof(fields) / sizeof(fields[0]));

    // k*, dis*, nue* across the phi seam (see the note at the top of this file). The two faces
    // are averaged and set equal, as every other field here is; the floor is applied to the
    // average rather than to each face, so the two stay exactly equal.
    Array* turb[] = { &tke, &dis, &nue };
    const double turb_floor[] = { 0.0, bc_dis_min, 0.0 };
    const int nt = (int)(sizeof(turb) / sizeof(turb[0]));
    const bool turb_bc = turb_active;

    #pragma omp parallel for
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int f = 0; f < nf; f++){
                Array& F = *fields[f];
                F.x[i][j][0] = c43 * F.x[i][j][1] - c13 * F.x[i][j][2];
                F.x[i][j][km-1] = c43 * F.x[i][j][km-2] - c13 * F.x[i][j][km-3];
                F.x[i][j][0] = F.x[i][j][km-1] = (F.x[i][j][0] + F.x[i][j][km-1])/2.0;
            }

            if(turb_bc){
                for(int f = 0; f < nt; f++){
                    Array& F = *turb[f];
                    F.x[i][j][0] = c43 * F.x[i][j][1] - c13 * F.x[i][j][2];
                    F.x[i][j][km-1] = c43 * F.x[i][j][km-2] - c13 * F.x[i][j][km-3];
                    F.x[i][j][0] = F.x[i][j][km-1] =
                        std::max(turb_floor[f], (F.x[i][j][0] + F.x[i][j][km-1])/2.0);
                }
            }
        }
    }

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
