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
    // The 39 fields every model shares are the SHARED Reporting<Planet>; see
    // print_minmax_common() there for why they could not be shared until the units were
    // settled. Below is only what is ATSAT's own.
    Reporting<cSaturnModel>(*this).print_minmax_common();

    cout << endl;

    cout << endl << " Pressures " << endl;
//    searchMinMax_3D(" max 3D density ", " min 3D density ", " kg/m³", rho, 1.0);
    cout << endl;

    cout << endl << " Water " << endl;
//    searchMinMax_3D(" max 3D cloudiness_h2o ", " min 3D cloudiness_h2o ", "[/]", cloudiness_h2o, 1.0);
    cout << endl;

    // CH4 was never printed at all — the fields existed, condensed and (until t_00_ch4 was
    // corrected) accumulated without a sink, and none of it appeared in any log. That is how a
    // 45.6 g/m3 methane ice deck with zero methane snow stayed invisible until someone read the
    // .vtk by hand.
    cout << endl;

    cout << endl << " Ammonia " << endl;
//    searchMinMax_3D(" max 3D cloudiness_nh3 ", " min 3D cloudiness_nh3 ", "[/]", cloudiness_nh3, 1.0);
    cout << endl;

    cout << endl;

    cout << endl << " Energies " << endl;

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

    // Equatorial column profile (j=jm/2, k=km/2), top -> bottom, for a direct check of the
    // radiation / Q_rad fields against the actual T(p). ATSAT was the last model without it.
    //
    // It is the table that caught ATNEPT's 607x pressure defect after the fact: that model's
    // photosphere sat at 17.2287 bar and was read as an opacity failure across two commits, when
    // what was actually wrong was init_PressureStatic putting the TOP of the domain at 15 bar. A
    // single column of p[bar] shows that at a glance and no column-MEAN diagnostic can.
    //
    // SATURN HAS THE OPEN QUESTION THIS IS MOST LIKELY TO BEAR ON. Its photosphere sits at
    // 0.0660 bar against Jupiter's 0.3263 with the same Jupiter-tuned opacity, and it UNDER-emits
    // — 1.306 W/m2 against the 4.450 that enters it, with a photosphere 31.0 K colder than
    // T_eff(in). Whether that is Saturn or a symptom is the question Radiation.h's banner has
    // carried since the shared port, and eps against p[bar] down the column is the direct view of
    // the opacity that question is about. The p[bar] column is worth having with the radiation
    // knob off as well, which is why the table is not gated on it.
    {
        const int j0 = jm / 2, k0 = km / 2;
        cout << endl << " Equatorial column  (j=" << j0 << ", k=" << k0
             << ")   top -> bottom" << endl;
        printf("   %3s  %10s  %8s  %8s  %12s  %14s\n",
               "i", "p[bar]", "T[K]", "eps", "netRad[W/m2]", "Q_rad[W/m3]");
        for(int i = im - 1; i >= 0; i--){
            printf("   %3d  %10.4f  %8.2f  %8.4f  %12.4f  %14.4e\n",
                   i, p_stat.x[i][j0][k0], t.x[i][j0][k0] * t_ref,
                   epsilon.x[i][j0][k0], radiation.x[i][j0][k0], Q_rad.x[i][j0][k0]);
        }
    }
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
