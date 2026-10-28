//
// ATSAT's binding of the SHARED BoundaryConditions<Planet> template.
//
// The three routines themselves — bcRadius, bcTheta, bcPhi, and the five hardening knobs — are in
// BoundaryConditions.h, which ATJUP instantiates too. What is here is the part that is ATSAT's
// alone: WHICH FIELDS each boundary treats. The lists are the whole of the difference between
// the two models' boundary conditions, which is the point of having moved them here.
//
// ===== BOUNDARY TREATMENT OF THE TURBULENCE FIELDS =====
//
// k* and dis* are prognostic variables of the RK4 system (RungeKutta_Sat_Turb.cpp), so their
// domain boundaries must be treated like every other transported field. ATJUP records the
// consequence of not doing so at BC_Jup.h:190 — the interior evolves away from the frozen faces,
// the resulting permanent Laplacian at i=im-1, the poles and the phi seam drains dis* to its
// floor within ~20 iterations, nue* = k*/omega* then saturates its ceiling, and with the coupling
// on that eddy viscosity destroys the momentum field. So this is a prerequisite for
// ATSAT_TURB_COUPLING, not an improvement on it.
//
// They are given their own list rather than being added to the ones below, for two reasons the
// shared file records in full: the 2-point Neumann form is used instead of the cubic, because the
// cubic amplifies an alternating error 7x per call and dis* appears in denominators throughout
// the closure; and the results are clamped to the bounds the RK4 stages enforce, k* >= 0 and
// dis* >= dis_min, because an extrapolation is free to produce a negative where the integration
// is not, and a negative dis* at one face is an infinite nue* at the next call of the closure.
//
// i=0 is also written, though TurbulenceSat::apply_wall_bc() re-imposes its own zero-gradient
// condition there later in the same iteration. The duplication is what keeps tken/disn — which
// restoreVar copies from tke/dis between the two — from carrying a stale deep boundary.
//
// Gated by turb_active. With the closure off tke, dis and nue are identically zero, so the
// extrapolation would be a no-op, and skipping it keeps the off path exactly as it was.

#include "cSaturnModel.h"
#include "BC_Sat.h"

using namespace std;

namespace {
    constexpr double bc_dis_min = 1.0e-10;   // matches TurbulenceSat::dis_min and the RK4 floor
}

// ATSAT_BC_SPECIES_LID_COPY (2026-10-11; DEFAULT 1 on the user's word, SAT-LID; 0 restores): the eleven species on the LID
// copy the cell below instead of taking the shared routine's form. ATSAT's radial form is the
// 3-point cubic f[s] = f[c] - 3 f[b] + 3 f[a], which is exact for a straight profile -- the right
// thing for the temperature at the bottom -- and overshoots wherever a field bends towards the
// boundary. MEASURED (satchk/giants/S1000, iteration 1000): ch4_cloud 0.093 / 0 / 0.70 g/m3 in the
// three cells below the lid at 2 S gives 2.2 on it (the run's "max ch4_cloud", growing without
// levelling), ch4_ice 5.66 / 12.98 / 23.26 gives 36.5, and ch4 116 / 80 / 40 gives -4, which the
// clamp then zeroes. A copy cannot overshoot. The bottom plane and every other field keep the
// shared form. The lid values of v and w are ATSAT_BC_TOP_TAPER's business, not this knob's.
// 224 iterations with this knob and ATSAT_BC_TOP_TAPER=2 against without (satchk/giants/LID3
// against SYNC2): max ch4_ice 93.2 g/m3 ON the lid -> 42.2 one cell below it, min w -16.5 m/s on
// the lid -> -12.8 at 412 km, the temperature untouched, 60 of 100 printed extrema identical.
// Changing the whole radial form instead (ATSAT_BC_RADIUS_COPY = 1 or 2, LID1 / LID2) removed the
// overshoot as well but flattened the bottom temperature by 10 to 13 K.
void BC_Sat::bcRadius() {
    BoundaryConditions<cSaturnModel>(m).bcRadius();

    static const int lid_copy = [](){
        const char* e = getenv("ATSAT_BC_SPECIES_LID_COPY"); return e ? atoi(e) : 1; }();
    // ATSAT_BC_SPECIES_WALL_COPY (2026-10-11; DEFAULT 1 on the user's word, 0 restores): the same
    // copy on the BOTTOM plane, ATJUP's bottom wall (ATJUP_BC_SPECIES_WALL_COPY). An extrapolated
    // wall value lets the centred diffusion stencil at i = 1 draw the species through a wall that
    // is meant to be closed. On ATSAT there is next to nothing for it to draw today: the level
    // means of total H2O, NH3, CH4 and H2S are flat over i = 0..5 (satchk/giants/S1000pf,
    // iteration 1000: q0 - q1 is 3e-8 to 2e-6 of the column, against 1e-3 on ATJUP), so this
    // closes the wall before a profile forms there, and no gain from it has been measured.
    // =2 (2026-10-28, row GIANTS-MCA-DEFAULT; DEFAULT 2 since 2026-10-28, user: "do (a) and (b) for ATSAT";
    // 1 restores the density copy): copy the MASS FRACTION instead of
    // the density, x[0] = x[1] * rho0 / rho1 with rho_at(). Why: =1 copies kg/m3, so under
    // ATSAT_LOCAL_RHO=1 the denser floor cell holds rho1/rho0 = 0.9426 of the mass fraction above
    // it (all gases alike, satchk/mcagiants_2026-10-28/dbg) -- lighter air under heavier, which the
    // Ledoux test of the MCA fires on in every column (65341 dry triggers in call 1; 0 on r_mix).
    // With ATSAT_LOCAL_RHO=0 rho_at() is r_mix on both cells and =2 is exactly =1.
    static const int wall_copy = [](){
        const char* e = getenv("ATSAT_BC_SPECIES_WALL_COPY"); return e ? atoi(e) : 2; }();
    if(lid_copy == 0 && wall_copy == 0) return;
    Array* species[] = {
        &m.h2o, &m.h2o_cloud, &m.h2o_ice,
        &m.h2s,
        &m.nh3, &m.nh3_cloud, &m.nh3_ice,
        &m.ch4, &m.ch4_cloud, &m.ch4_ice,
        &m.nh4sh };
    const int ns = (int)(sizeof(species) / sizeof(species[0]));
    const int im = m.im, jm = m.jm, km = m.km;
    #pragma omp parallel for
    for(int j = 0; j < jm; j++)
        for(int k = 0; k < km; k++){
            const double fr = (wall_copy == 2) ? m.rho_at(0, j, k) / m.rho_at(1, j, k) : 1.0;
            for(int f = 0; f < ns; f++){
                if(lid_copy  != 0) species[f]->x[im-1][j][k] = species[f]->x[im-2][j][k];
                if(wall_copy != 0) species[f]->x[0][j][k]    = species[f]->x[1][j][k] * fr;
            }
        }
}
void BC_Sat::bcTheta()  { BoundaryConditions<cSaturnModel>(m).bcTheta();  }
void BC_Sat::bcPhi()    { BoundaryConditions<cSaturnModel>(m).bcPhi();    }


// Every transported field is extrapolated at the deep boundary i=0 and the model top i=im-1.
std::vector<Array*> cSaturnModel::bc_fields_radius(){
    return {
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
}

// At the poles, u IS extrapolated — it is tangential there — while v and w are not.
// NOTE ATSAT extrapolates massflux_* and difflux_* here where ATJUP zeroes them; see the note in
// BoundaryConditions.h. That difference is preserved, not resolved.
std::vector<Array*> cSaturnModel::bc_fields_theta_extrap(){
    return {
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
}

// v, w and the nh4sh flux limiter are pinned to zero at both poles: the meridional and zonal
// velocities have no meaning on the axis, and the limiter corrects an advective flux that does
// not exist there.
std::vector<Array*> cSaturnModel::bc_fields_theta_zero(){
    return { &v, &w, &fluxlim_nh4sh };
}

std::vector<Array*> cSaturnModel::bc_fields_phi(){
    return {
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
}

std::vector<Array*> cSaturnModel::bc_turb_fields(){ return { &tke, &dis, &nue }; }
std::vector<double> cSaturnModel::bc_turb_floors(){ return { 0.0, bc_dis_min, 0.0 }; }
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
