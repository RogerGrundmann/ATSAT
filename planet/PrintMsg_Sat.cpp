#include "cSaturnModel.h"
#include "Reporting.h"

using namespace std;

/*
* The min/max report machinery and steadyQuery now live in the SHARED Reporting.h, as
* Reporting<Planet>. What stayed here is what is genuinely Saturn's: which fields printMinMax
* lists, in what unit, and the welcome/final text. See Reporting.h for why the split falls there.
*/
/*
*
*/
/*
 * UNSETTLED, and left as it stands rather than changed under a refactor: the species units.
 *
 * ATJUP prints every species with coefficient 1.0 and labels it kg/m3, having found that its own
 * 1e3*r_mix factor made the printed number 1.2844x too large — c is already a mass density in
 * kg/m3. The lines below still carry 1e3*r_mix and the label g/m3, so if ATJUP's reading is right
 * then every species number in ATSAT's log is 1284x its true value in those units.
 *
 * The temperature line has the same shape of problem. It scales t by 273.15 and then subtracts
 * 273.15, but t is nondimensional and t_ref (134.0) is its scale, so the printed degC is not the
 * model's temperature. ATJUP corrected exactly this in its own copy (its note: "it printed the
 * offset from t_ref while labelling it degC — every temperature in the log read ~108 K too
 * high"). ATSAT_TRACE, which uses t*t_ref, reports a plausible 30-400 K against this block.
 *
 * Both are OUTPUT-ONLY — no prognostic field depends on them — which is why correcting them can
 * wait for its own commit and its own measurement instead of riding along with the sharing.
 */
void cSaturnModel::printMinMax(){

    cout << endl << endl << " Courant time step   dt = " << dt << endl << endl;

    cout << endl << endl << " Temperatures " << endl;
    searchMinMax_3D(" max 3D temperature ", " min 3D temperature ", " degC",
        t, t_ref, [](double i)->double{return i - 273.15;}, true);
    searchMinMax_3D(" max 3D thermalflux ", " min 3D thermalflux ", " W/m3", thermalmassflux, 1.0);
    cout << endl;

    cout << endl << " Velocities " << endl;
    searchMinMax_3D(" max 3D u-component ", " min 3D u-component ", " m/s", u, u_0);
    searchMinMax_3D(" max 3D v-component ", " min 3D v-component ", " m/s", v, u_0);
    searchMinMax_3D(" max 3D w-component ", " min 3D w-component ", " m/s", w, u_0);
    cout << endl;

    cout << endl << " Pressures " << endl;
    searchMinMax_3D(" max 3D pressure dynamic ", " min 3D pressure dynamic ", " bar", p_dyn, p_dyn_to_bar());
    searchMinMax_3D(" max 3D pressure static ", " min 3D pressure static ", " bar", p_stat, 1.0);
//    searchMinMax_3D(" max 3D density ", " min 3D density ", " kg/m³", rho, 1.0);
    cout << endl;

    cout << endl << " Water " << endl;
    searchMinMax_3D(" max 3D h2o ",  " min 3D h2o ", " kg/m3", h2o, 1.0);
    searchMinMax_3D(" max 3D h2o_cloud ", " min 3D h2o_cloud ", " kg/m3", h2o_cloud, 1.0);
    searchMinMax_3D(" max 3D h2o_ice ", " min 3D h2o_ice ", " kg/m3", h2o_ice, 1.0);
//    searchMinMax_3D(" max 3D cloudiness_h2o ", " min 3D cloudiness_h2o ", "[/]", cloudiness_h2o, 1.0);
    cout << endl;

    // CH4 was never printed at all — the fields existed, condensed and (until t_00_ch4 was
    // corrected) accumulated without a sink, and none of it appeared in any log. That is how a
    // 45.6 g/m3 methane ice deck with zero methane snow stayed invisible until someone read the
    // .vtk by hand.
    cout << endl << " Methane " << endl;
    searchMinMax_3D(" max 3D ch4 ",  " min 3D ch4 ", " kg/m3", ch4, 1.0);
    searchMinMax_3D(" max 3D ch4_cloud ", " min 3D ch4_cloud ", " kg/m3", ch4_cloud, 1.0);
    searchMinMax_3D(" max 3D ch4_ice ", " min 3D ch4_ice ", " kg/m3", ch4_ice, 1.0);
    cout << endl;

    cout << endl << " Hydrogen Sulfide " << endl;
    searchMinMax_3D(" max 3D h2s ",  " min 3D h2s ", " kg/m3", h2s, 1.0);
    searchMinMax_3D(" max 3D w_h2s ", " min 3D w_h2s ", " kg/(m3s)", w_h2s, 1.0);
    searchMinMax_3D(" max 3D j_h2s ", " min 3D j_h2s ", " kg/m4", j_h2s, 1.0);
    searchMinMax_3D(" max 3D jT_h2s ", " min 3D jT_h2s ", " kg/m4", jT_h2s, 1.0);
    searchMinMax_3D(" max 3D massflux_h2s ", " min 3D massflux_h2s ", " kg/(m3s)", massflux_h2s, 1.0);
    searchMinMax_3D(" max 3D diff_h2s ", " min 3D diff_h2s ", " kg/(m3s)", difflux_h2s, 1.0);
    cout << endl;

    cout << endl << " Ammonia " << endl;
    searchMinMax_3D(" max 3D nh3 ",  " min 3D nh3 ", " kg/m3", nh3, 1.0);
    searchMinMax_3D(" max 3D nh3_cloud ", " min 3D nh3_cloud ", " kg/m3", nh3_cloud, 1.0);
    searchMinMax_3D(" max 3D nh3_ice ", " min 3D nh3_ice ", " kg/m3", nh3_ice, 1.0);
    searchMinMax_3D(" max 3D w_nh3 ", " min 3D w_nh3 ", " kg/(m3s)", w_nh3, 1.0);
    searchMinMax_3D(" max 3D j_nh3 ", " min 3D j_nh3 ", " kg/m4", j_nh3, 1.0);
    searchMinMax_3D(" max 3D jT_nh3 ", " min 3D jT_nh3 ", " kg/m4", jT_nh3, 1.0);
    searchMinMax_3D(" max 3D massflux_nh3 ", " min 3D massflux_nh3 ", " kg/(m3s)", massflux_nh3, 1.0);
    searchMinMax_3D(" max 3D diff_nh3 ", " min 3D diff_nh3 ", " kg/(m3s)", difflux_nh3, 1.0);
//    searchMinMax_3D(" max 3D cloudiness_nh3 ", " min 3D cloudiness_nh3 ", "[/]", cloudiness_nh3, 1.0);
    cout << endl;

    cout << endl << " Ammonia Hydrosufide " << endl;
    searchMinMax_3D(" max 3D nh4sh ",  " min 3D nh4sh ", " mg/m3", nh4sh, 1e6);
    searchMinMax_3D(" max 3D w_nh4sh ", " min 3D w_nh4sh ", " mg/(m3s)", w_nh4sh, 1e6);
    searchMinMax_3D(" max 3D j_nh4sh ", " min 3D j_nh4sh ", " mg/m4", j_nh4sh, 1e6);
    searchMinMax_3D(" max 3D jT_nh4sh ", " min 3D jT_nh4sh ", " mg/m4", jT_nh4sh, 1e6);
    searchMinMax_3D(" max 3D massflux_nh4sh ", " min 3D massflux_nh4sh ", " mg/(m3s)", massflux_nh4sh, 1e6);
    searchMinMax_3D(" max 3D diff_nh4sh ", " min 3D diff_nh4sh ", " mg/(m3s)", difflux_nh4sh, 1e6);
    cout << endl;

    cout << endl << " Forces " << endl;
    searchMinMax_3D(" max 3D Coriolis force ", " min 3D Coriolis force ", " mN/m3", CoriolisForce, 1e3);
    searchMinMax_3D(" max 3D centrifugal force ", " min 3D centrifugal force ", " mN/m3", CentrifugalForce, 1e3);
    searchMinMax_3D(" max 3D buoyancy force ", " min 3D buoyancy force ", " N/m3", BuoyancyForce, 1.0);
    searchMinMax_3D(" max 3D presgrad force ", " min 3D presgrad force ", " N/m3", PresGradForce, 1.0);

    cout << endl << " Energies " << endl;
    searchMinMax_3D(" max 3D sensible heat ", " min 3D sensible heat ", " W/m3", Q_Sensible, 1.0);
    searchMinMax_3D(" max 3D latent heat ", " min 3D latent heat ", " W/m3", Q_Latent, 1.0);

    // Radiation and turbulence diagnostics (ported classes). All identically zero unless the
    // corresponding knob is set, so this block costs nothing when they are off — but without it
    // the ported fields are invisible and the ports cannot be judged.
    searchMinMax_3D(" max 3D net radiation ", " min 3D net radiation ", "W/m2", radiation, 1.0);
    searchMinMax_3D(" max 3D Q_rad ", " min 3D Q_rad ", "W/m3", Q_rad, 1.0);
    searchMinMax_3D(" max 3D emissivity ", " min 3D emissivity ", "/", epsilon, 1.0);
    searchMinMax_3D(" max 3D tke ", " min 3D tke ", "/", tke, 1.0);
    searchMinMax_3D(" max 3D dis ", " min 3D dis ", "/", dis, 1.0);
    searchMinMax_3D(" max 3D nue ", " min 3D nue ", "/", nue, 1.0);
    searchMinMax_3D(" max 3D prod ", " min 3D prod ", "/", prod, 1.0);
    searchMinMax_3D(" max 3D P_rain ", " min 3D P_rain ", "kg/m2/s", P_rain, 1.0);
    searchMinMax_3D(" max 3D P_nh3_rain ", " min 3D P_nh3_rain ", "kg/m2/s", P_nh3_rain, 1.0);
    searchMinMax_3D(" max 3D Q_precip ", " min 3D Q_precip ", "W/m3", Q_precip, 1.0);
    searchMinMax_3D(" max 3D P_ch4_rain ", " min 3D P_ch4_rain ", "kg/m2/s", P_ch4_rain, 1.0);
    cout << endl << endl;

    reportClampBudget();
}
/*
*
*/
/*
*
*/
// Forwarders to the shared Reporting<cSaturnModel>. The bodies used to be here in full and were
// the highest-overlap pair between the two models (83.8 % line-for-line).
void cSaturnModel::searchMinMax_3D(string name_maxValue, string name_minValue,
    string name_unitValue, Array &value_D, double coeff,
    std::function< double(double) > lambda, bool print_heading){
    Reporting<cSaturnModel>(*this).searchMinMax_3D(name_maxValue, name_minValue,
        name_unitValue, value_D, coeff, lambda, print_heading);
}
/*
*
*/
void cSaturnModel::searchMinMax_2D(string name_maxValue, string name_minValue,
    string name_unitValue, Array_2D &value, double coeff){
    Reporting<cSaturnModel>(*this).searchMinMax_2D(name_maxValue, name_minValue,
        name_unitValue, value, coeff);
}
/*
*
*/
void cSaturnModel::print_welcome_msg(){
    if(verbose){
        cout << endl << endl << endl;
        cout << "***** Atmosphere Saturn General Circulation Model(ATSAT) applied to laminar flow" << endl;
        cout << "***** program for the computation of Saturn-atmospherical circulating flows in a spherical shell" << endl;
        cout << "***** finite difference scheme for the solution of the 3D Navier-Stokes equations" << endl;
        cout << "***** with 6 additional transport equations to describe the water vapour, cloud water, cloud ice and nh3 vapour, nh3 cloud and nh3 ice" << endl;
        cout << "***** 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop" << endl;
        cout << "***** Poisson equation for the pressure solution in an outer iterational loop" << endl;
        cout << "***** temperature distribution given as a parabolic distribution from pole to pole, zonally constant" << endl;
        cout << "***** water and nh3 vapour distribution given by Clausius-Claperon equation for the partial pressure" << endl;
        cout << "***** water vapour is part of the Boussinesq approximation and the absorptivity in the radiation model" << endl;
        cout << "***** two category ice scheme for cold clouds applying parameterization schemes provided by the COSMO code(German Weather Forecast)" << endl;
        cout << "***** rain and snow precipitation solved by column equilibrium applying the diagnostic equations" << endl;
        cout << "***** code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz(roger.grundmann@web.de)" << endl << endl;
        cout << "***** original program name:  " << __FILE__ << endl;
        cout << "***** compiled:  " << __DATE__  << "  at time:  " << __TIME__ << endl << endl;
        has_welcome_msg_printed = true;
    }
    return;
}
/*
*
*/
void cSaturnModel::initMsg(){
    cout << "  present state of the computation " << endl << "  current number of iterations nm = " << nm << "    iter_n = " << iter_n << endl;
    return;
}
/*
*
*/
void cSaturnModel::print_final_msg(){
    cout << endl << "***** end of the SaturnAtmosphere General Circulation Modell(ATSAT) *****" << endl << endl;
    if(iter_n == nm)   cout <<  "***** number of artificial time steps      iter_n = " << iter_n << ", end of program reached because of limit of maximum artificial time steps ***** \n\n" << endl;
}
/*
*
*/
// See Reporting.h for what this reports and the three defects that were fixed when it was
// revived. It must run BEFORE restoreVar(), or every number it prints is identically zero.
void cSaturnModel::steadyQuery(){
    Reporting<cSaturnModel>(*this).steadyQuery();
}
