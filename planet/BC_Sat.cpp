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

void BC_Sat::bcRadius() { BoundaryConditions<cSaturnModel>(m).bcRadius(); }
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
