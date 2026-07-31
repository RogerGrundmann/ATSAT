/*
 * Saturn Atmosphere Circulation Model (JACM) applied to laminar flow
 * program for the computation of saturn-atmospherical circulating flows in a spherical shell
 * modeling of the atmosphere with gases, their cloud and ice formation: H2O, H2S, NH3 and NH4SH
 * finite difference scheme for the solution of the 3D Navier-Stokes equations
 * with 1 additional transport equations to describe the salinity
 * 4th order Runge-Kutta scheme to solve 2nd order differential equations inside an inner iterational loop
 * Poisson equation for the pressure solution in an outer iterational loop
 * temperature distribution given as a parabolic distribution from pole to pole, zonaly constant
 * code developed by Roger Grundmann, Zum Marktsteig 1, D-01728 Bannewitz(roger.grundmann@web.de)
*/

#include "cSaturnModel.h"
#include "ChemistrySat.h"
#include "PressureSolverSat.h"
#include "SaturationAdjustmentSat.h"
#include "BC_Sat.h"
#include "VelocityInitializerSat.h"
#include "ConvectiveAdjustmentSat.h"
#include "RadiationSat.h"
#include "ThermalWindDiagSat.h"
#include "PrecipitationSat.h"
#include "TurbulenceSat.h"

using namespace std;
using namespace tinyxml2;
using namespace AtomUtils;

cSaturnModel* cSaturnModel::m_model = NULL;

const double cSaturnModel::pi180 = 180.0/M_PI;      // pi180 = 57.3

const double cSaturnModel::the_degree = 1.0;         // compares to 1° step size laterally
const double cSaturnModel::phi_degree = 1.0;         // compares to 1° step size longitudinally

const double cSaturnModel::the0 = 0.0;             // North Pole
const double cSaturnModel::phi0 = 0.0;             // zero meridian in Greenwich

const double cSaturnModel::r0 = 1.0; // value much too small, Saturn radius 72000km

const double cSaturnModel::dr = 0.025;    // 0.025 x 40 = 1.0 compares to 16 km : 40 = 400 m for 1 radial step
const double cSaturnModel::dthe = the_degree/pi180; 
const double cSaturnModel::dphi = phi_degree/pi180;


cSaturnModel::cSaturnModel():
    j_ellipse(std::vector<std::vector<int> >(jm, std::vector<int>(km, 0))),
    has_welcome_msg_printed(false){
//    if(PythonStream::is_enable()){
//        backup = std::cout.rdbuf();
//        std::cout.rdbuf(&ps);
//    }
    // If Ctrl-C is pressed, quit
    signal(SIGINT, exit);
    // set default configuration
    SetDefaultConfig();
    m_model = this;
    rad.initArray_1D(im, 0); // radial coordinate direction
    the.initArray_1D(jm, 0); // lateral coordinate direction
    phi.initArray_1D(km, 0); // longitudinal coordinate direction
    rad.Coordinates(im, r0, dr);
    the.Coordinates(jm, the0, dthe);
    phi.Coordinates(km, phi0, dphi);
    init_layer_heights();
}

cSaturnModel::~cSaturnModel(){
    delete m_chem;
    m_chem = nullptr;
    delete m_pressure;
    m_pressure = nullptr;
    if(PythonStream::is_enable()){
        std::cout.rdbuf(backup);
    }
    m_model = NULL;
}


// Dry convective adjustment (ConvectiveAdjustmentSat), ported from ATJUP. Default OFF, so every
// existing ATSAT run stays bit-identical. It restores any superadiabatic column to the dry adiabat
// while conserving the column's mass-weighted enthalpy; nothing else in the model does that.
// Whether ATSAT develops such columns at all is unmeasured — switching this on and reading the
// per-iteration report (columns touched, layers mixed, max dT, enthalpy drift) is how to find out.
static int conv_adj_enabled(){
    static const int v = [](){ const char* e = getenv("ATSAT_CONV_ADJ"); return e ? atoi(e) : 0; }();
    return v;
}


// Grey multi-layer radiation (RadiationSat), ported from ATJUP. Default OFF, bit-identical when
// off. It fills the DIAGNOSTIC arrays radiation / epsilon / Q_rad and touches neither t nor any
// rhs — wiring the heating into the temperature equation is a separate step. Runs on even
// iterations, alongside the rest of the physics block.
static int radiation_enabled(){
    static const int v = [](){ const char* e = getenv("ATSAT_RADIATION"); return e ? atoi(e) : 0; }();
    return v;
}


// Thermal-wind residual (ThermalWindDiagSat), a MEASUREMENT only — it modifies nothing.
// ATSAT_TW_DIAG, default 0.
static int tw_diag_enabled(){
    static const int v = [](){ const char* e = getenv("ATSAT_TW_DIAG"); return e ? atoi(e) : 0; }();
    return v;
}


// Timestep (ATSAT_DT, nondimensional; unset keeps the formula below unchanged).
//
// The default is dt = 2.8284 * dr/u_0 * 0.2, and 2.8284 = 2*sqrt(2) is RK4's stability limit on
// the imaginary axis with 0.2 as a safety factor — so the intent is a proper CFL condition. The
// evaluation is not one: dr is DIMENSIONLESS (0.025) while u_0 is in m/s (470), and a
// dimensionless length divided by a dimensional velocity is not a Courant number. It gives
// dt = 3.01e-5, and with the time unit L/u_0 = 500 km / 470 m/s = 1064 s that is 0.032 s of
// Saturn time per iteration: the configured 224 iterations span SEVEN SECONDS. Reaching even a
// ten-minute window would take ~19700 iterations and one Saturn rotation about 1.2 million, so
// nothing horizontal can develop in any run of practical length and every measurement made on
// such a run is a statement about the initial state.
//
// Read as an actual CFL condition, with the NONDIMENSIONAL velocity (max|w|/u_0 = 268/470 = 0.57)
// in place of u_0, the same formula gives
//     dt = 2.8284 * 0.2 * dr / 0.57 = 0.025,
// which is 26 s per iteration and makes 224 iterations 1.6 hours — a factor of 825. The identical
// unit mix was found in ATJUP's dt on 2026-07-29, and there the corresponding value ran cleanly.
//
// The default is left alone so every existing run is bit-identical; this knob is how a longer
// horizon gets measured. Each run prints what it is actually covering.
static double timestep_override(){
    static const double v = [](){ const char* e = getenv("ATSAT_DT"); return e ? atof(e) : 0.0; }();
    return v;
}


// Precipitation microphysics (PrecipitationSat), ported from ATJUP. Default OFF, bit-identical.
// It DOES feed back: the condensate it converts is removed from the cloud/ice fields in place, so
// it must run AFTER the SaturationAdjustmentSat calls or the adjustment would simply undo it.
//
// WHY IT STAYS OFF, measured 2026-07-31 over 50 iterations (config_m50, 12 threads). With it on
// the scheme runs and produces fluxes, but they are the SAFETY CAP and not physics. At the
// i=20 level (250 km), reading the .vtk directly:
//
//     P_snow    8.6400 mm/day in 65337 of 65341 cells  <- P_max_flux exactly, everywhere
//     P_graupel 2.1243 mm/day max     P_rain 0 (too cold there)     P_nh3_rain 0.0131 max
//     Precip_total at the surface: 0 — nothing arrives at the bottom
//
// 8.64 mm/day is 1.0e-4 kg/m2/s, the P_max_flux clamp in PrecipitationSat::column. ATJUP's note
// on that constant says it is set "well above the energy-budget scale, so the cap guards against
// a numerical runaway instead of silently becoming the operative limiter as it used to". On
// Saturn it IS the operative limiter, over an entire model level.
//
// The cause is not the cap but the coefficients around it: c_c_au and the rest are ATOM's
// terrestrial numbers, rescaled ONCE to Jupiter's energy budget and then carried here unchanged,
// while ATSAT's condensate loading is far higher (max h2o_cloud 113 g/m3 against ATJUP's ~30).
// Switching this on before they are recalibrated to Saturn's own budget would put a saturated
// clamp into the moisture equations and call it precipitation. Recalibrate first: the observable
// is the same one ATJUP used, Lv*P against Saturn's emitted flux.
static int precip_enabled(){
    static const int v = [](){ const char* e = getenv("ATSAT_PRECIP"); return e ? atoi(e) : 0; }();
    return v;
}


// Shapiro filtering of the velocity fields (ATJUP's dampVelocities). ATSAT's Utils has the
// 1-2-1 filter but not ATJUP's higher-order damp_wiggles_ho, so only the 2nd-order form is
// available here; ATJUP_VEL_SHAPIRO_ORDER=4 has no ATSAT counterpart.
static int shapiro_vel_inloop(){
    static const int v = [](){
        const char* e = getenv("ATSAT_VEL_SHAPIRO_INLOOP"); return e ? atoi(e) : 0; }();
    return v;
}
static double shapiro_strength(){
    static const double v = [](){
        const char* e = getenv("ATSAT_SHAPIRO_STRENGTH"); return e ? atof(e) : 1.0; }();
    return v;
}
static int nancheck_on(){
    static const int v = [](){
        const char* e = getenv("ATSAT_NANCHECK"); return e ? atoi(e) : 0; }();
    return v;
}


// Turbulence closure (TurbulenceSat), ported from ATJUP. Default OFF, bit-identical: with it off
// tke/dis/nue and the rest stay identically zero.
//
// Which model runs is now a CONFIGURATION choice, turb_model in config_atsat.xml
// (none | k_epsilon | k_omega | k_omega_SST, declared in param.py, default k_omega_SST), with
// ATSAT_TURB_MODEL overriding it at runtime. Both are folded into cSaturnModel::turb_active where
// the override is applied, and that flag is the only gate the rest of the model reads.
//
// nue* reaches the momentum and scalar equations when ATSAT_TURB_COUPLING is also set
// (RHS_Sat_Turb.cpp). The closure remains UNVALIDATED against anything observed on Saturn.
static int turb_env_enabled(){
    static const int v = [](){ const char* e = getenv("ATSAT_TURB"); return e ? atoi(e) : 0; }();
    return v;
}

#include "cSaturnDefaults.cpp.inc"
/*
*
*/
void cSaturnModel::dampVelocities(){
    const double s = shapiro_strength();
    AtomUtils::damp_wiggles(u, nullptr, true, true, true, s);
    AtomUtils::damp_wiggles(v, nullptr, true, true, true, s);
    AtomUtils::damp_wiggles(w, nullptr, true, true, true, s);
}

ChemistrySat& cSaturnModel::getChemistry(){
    if(!m_chem)
        m_chem = new ChemistrySat(*this);
    return *m_chem;
}

PressureSolverSat& cSaturnModel::getPressureSolver(){
    if(!m_pressure)
        m_pressure = new PressureSolverSat(*this);
    return *m_pressure;
}

void cSaturnModel::ChemMassRateSat()    { getChemistry().ChemMassRateSat(); }
void cSaturnModel::DiffMassFluxSat()    { getChemistry().DiffMassFluxSat(); }
void cSaturnModel::ThermalPropertiesSat(){ getChemistry().ThermalPropertiesSat(); }
/*
*
*/
void cSaturnModel::LoadConfig(const char *filename){
    XMLDocument doc;
    XMLError err = doc.LoadFile(filename);
    if(err){
        doc.PrintError();
        throw std::invalid_argument("   couldn't load config file inside cSaturnModel");
    }
    XMLElement *atsat = doc.FirstChildElement("atsat");
    if(!atsat){
        return;
    }
    XMLElement* elem_common = doc.FirstChildElement("atsat")->FirstChildElement("common");
    if(!elem_common){
        return;
    }
    XMLElement* elem_saturn = doc.FirstChildElement("atsat")->FirstChildElement("saturn");
    if(!elem_saturn){
        return;
    }
#include "SaturnLoadConfig.cpp.inc"
}
/*
*
*/
void cSaturnModel::Run(){
    // ATSAT_FPE=1 turns the first invalid floating-point operation into a SIGFPE instead of a
    // silently propagating NaN. Run the CLI under gdb to get the exact line:
    //     OMP_NUM_THREADS=1 ATSAT_FPE=1 gdb -batch -ex run -ex "bt 6" -ex "info locals"
    //         --args cli/sat . config_ca.xml   (one line; split here only for width)
    // (build with -g -O0 for line numbers). NaN is INVISIBLE to printMinMax, whose
    // searchMinMax_3D compares with a bare > and so skips every non-finite cell — which is why
    // a domain-wide NaN can look like a perfectly quiet run. Off by default: trapping would
    // abort on the first harmless inf in a diagnostic field.
    if(getenv("ATSAT_FPE")) feenableexcept(FE_INVALID | FE_DIVBYZERO);


    #ifdef _OPENMP
        printf("\n\n   number of processors: %d\n\n", omp_get_num_procs());

    #pragma omp parallel
        {
            printf("   thread %d of %d in \"Desktop-Dell-XPS 8960\"\n", 
                omp_get_thread_num(), omp_get_num_threads());
        }
    #else
        printf("   OpenMP is not supported\n");
    #endif
        printf("   ended\n\n");

    output_path = output_path + "-Saturn";
    mkdir(output_path.c_str(), 0777);
    m_model = this;
    cout.precision(6);
    cout.setf(ios::fixed);
    if(!has_welcome_msg_printed)
        print_welcome_msg();
    initMsg();

    iter_n = 0;
    panorama_cnt = 0;

    resetArrays();

    dt = 2.8284 * dr/u_0 * 0.2;
    if(timestep_override() > 0.0) dt = timestep_override();
    // Resolve the turbulence switch. Order matters: the configuration value is already loaded,
    // ATSAT_TURB_MODEL overrides it here, and only then is turb_active fixed — everything
    // downstream (TurbulenceSat, RHSSat, RungeKuttaSat, the BC routines) reads that flag and never
    // the string. "none" now means what it says; before, it fell through to k-omega SST.
    if(const char* tm = getenv("ATSAT_TURB_MODEL")) turb_model = tm;
    turb_active = (turb_env_enabled() != 0) && (turb_model != "none");
    if(turb_env_enabled() != 0 && !turb_active)
        cout << "      ATSAT: turbulence requested but turb_model = \"none\" — closure off"
             << endl;
    printf("      ATSAT: dt = %.6g nondimensional = %.4g s of Saturn time per iteration"
           " (%d iterations = %.4g s = %.3f %% of a rotation)\n",
           dt, dt * L_atm * 1.0e3 / u_0, nm, nm * dt * L_atm * 1.0e3 / u_0,
           100.0 * nm * dt * L_atm * 1.0e3 / u_0 / 38018.0);

    init_layer_heights();
    TropopauseLocation();
    init_tropopause_layers();
    VelocityInitializerSat(*this).compute();

    // ===== The two latitude bands the integrator never reaches (ATSAT_GHOST_BANDS) =====
    //
    // RungeKuttaSat runs j = 2 .. jm-3, and BC_theta writes only j = 0 and j = jm-1. So j = 1
    // and j = jm-2 are advanced by nothing and constrained by nothing: whatever the velocity
    // initialiser leaves there stays, bit for bit, to the end of the run. They are not initial
    // data — they are a permanent prescribed forcing on the edge of the computed domain, and no
    // amount of iterating can relax them.
    //
    // ATJUP has the same defect with six bands rather than two (its loop runs j = 3 .. jm-4) and
    // measured, at iteration 450, that its frozen rows held radial velocities five times larger
    // than anything in the integrated range. The report below is so that ATSAT's version of that
    // number is visible rather than inferred.
    //
    //   ATSAT_GHOST_BANDS=1   report what the frozen bands hold, change nothing
    //   ATSAT_GHOST_BANDS=2   report, then zero u,v,w there
    //
    // Zeroing is a DIAGNOSTIC, not a repair. The honest repair is to integrate the rows or to
    // give them a real polar boundary condition; this only measures what the prescribed forcing
    // is worth. Default 0 = every existing run bit-identical.
    {
        static const int ghost = [](){
            const char* e = getenv("ATSAT_GHOST_BANDS"); return e ? atoi(e) : 0; }();
        if(ghost > 0){
            const int rows[2] = {1, jm-2};
            double band_max = 0.0, inner_max = 0.0;
            for(int r = 0; r < 2; r++)
                for(int i = 0; i < im; i++)
                    for(int k = 0; k < km; k++){
                        band_max = std::max(band_max, std::fabs(u.x[i][rows[r]][k]));
                        band_max = std::max(band_max, std::fabs(v.x[i][rows[r]][k]));
                        band_max = std::max(band_max, std::fabs(w.x[i][rows[r]][k]));
                    }
            for(int i = 0; i < im; i++)
                for(int j = 2; j < jm-2; j++)
                    for(int k = 0; k < km; k++){
                        inner_max = std::max(inner_max, std::fabs(u.x[i][j][k]));
                        inner_max = std::max(inner_max, std::fabs(v.x[i][j][k]));
                        inner_max = std::max(inner_max, std::fabs(w.x[i][j][k]));
                    }
            printf("      ATSAT: GHOST BANDS j=1 and j=%d are never integrated."
                   " max|vel| there = %.6g m/s, against %.6g m/s in the integrated range"
                   " (j=2..%d)\n",
                   jm-2, band_max * u_0, inner_max * u_0, jm-3);
            if(ghost >= 2){
                printf("      ATSAT: ATSAT_GHOST_BANDS=2 - zeroing u,v,w in j=1 and j=%d\n", jm-2);
                for(int r = 0; r < 2; r++)
                    for(int i = 0; i < im; i++)
                        for(int k = 0; k < km; k++)
                            u.x[i][rows[r]][k] = v.x[i][rows[r]][k] = w.x[i][rows[r]][k] = 0.0;
            }
        }
    }

    AtomUtils::damp_wiggles(u, nullptr, true, true, true);
    AtomUtils::damp_wiggles(v, nullptr, true, true, true);
    AtomUtils::damp_wiggles(w, nullptr, true, true, true);

//    goto Printout;

    ChemistrySat(*this).ThermalPropertiesSat();
    init_temperature();
    init_PressureStatic();
    init_PressureDynamic();

    init_vapour_cloud_ice("CH4", ch4_tropopause, coeff_ch4_A, coeff_ch4_B,
        coeff_ch4_A_i, coeff_ch4_B_i, t_0_ch4, t_00_ch4,
        ep_ch4, r_ch4, m_ch4,
        C_ch4, L0_ch4, R_ch4, del_alf_ch4, del_bet_ch4,
        ch4, ch4_cloud, ch4_ice, cloudiness_h2o);

    AtomUtils::damp_wiggles(ch4,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(ch4_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(ch4_ice,   nullptr, true, true, true);

    init_vapour_cloud_ice("H2O", h2o_tropopause, coeff_h2o_A, coeff_h2o_B,
        coeff_h2o_A_i, coeff_h2o_B_i, t_0_h2o, t_00_h2o,
        ep_h2o, r_h2o, m_h2o,
        C_h2o, L0_h2o, R_h2o, del_alf_h2o, del_bet_h2o,
        h2o, h2o_cloud, h2o_ice, cloudiness_h2o);

    AtomUtils::damp_wiggles(h2o,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(h2o_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(h2o_ice,   nullptr, true, true, true);

    init_h2s("H2S", h2s_tropopause, coeff_h2s_A, coeff_h2s_B,
        coeff_h2s_A_i, coeff_h2s_B_i, t_0_h2s, t_00_h2s,
        ep_h2s, r_h2s, m_h2s,
        C_h2s, L0_h2s, R_h2s, del_alf_h2s, del_bet_h2s, h2s);

    AtomUtils::damp_wiggles(h2s, nullptr, true, true, true);

    init_vapour_cloud_ice("NH3", nh3_tropopause, coeff_nh3_A, coeff_nh3_B,
        coeff_nh3_A_i, coeff_nh3_B_i, t_0_nh3, t_00_nh3,
        ep_nh3, r_nh3, m_nh3,
        C_nh3, L0_nh3, R_nh3, del_alf_nh3, del_bet_nh3,
        nh3, nh3_cloud, nh3_ice, cloudiness_nh3);

    AtomUtils::damp_wiggles(nh3,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(nh3_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(nh3_ice,   nullptr, true, true, true);

    init_nh4sh("NH4SH", nh4sh_tropopause, coeff_nh4sh_A, coeff_nh4sh_B,
        ep_nh4sh, r_nh4sh, m_nh4sh,
        C_nh4sh, L0_nh4sh, R_nh4sh,
        del_alf_nh4sh, del_bet_nh4sh, nh4sh);

    AtomUtils::damp_wiggles(nh4sh, nullptr, true, true, true);

//    goto Printout;

    SaturationAdjustmentSat(*this).run("CH4",
        coeff_ch4_A, coeff_ch4_B, coeff_ch4_A_i, coeff_ch4_B_i,
        t_0_ch4, t_00_ch4,
        ep_ch4, lv_ch4, ls_ch4, cp_ch4, r_ch4,
        C_ch4, L0_ch4, R_ch4, del_alf_ch4, del_bet_ch4, m_ch4,
        ch4, ch4_cloud, ch4_ice);

    AtomUtils::damp_wiggles(ch4,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(ch4_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(ch4_ice,   nullptr, true, true, true);

    SaturationAdjustmentSat(*this).run("H2O",
        coeff_h2o_A, coeff_h2o_B, coeff_h2o_A_i, coeff_h2o_B_i,
        t_0_h2o, t_00_h2o,
        ep_h2o, lv_h2o, ls_h2o, cp_h2o, r_h2o,
        C_h2o, L0_h2o, R_h2o, del_alf_h2o, del_bet_h2o, m_h2o,
        h2o, h2o_cloud, h2o_ice);

    AtomUtils::damp_wiggles(h2o,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(h2o_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(h2o_ice,   nullptr, true, true, true);

    SaturationAdjustmentSat(*this).run("NH3",
        coeff_nh3_A, coeff_nh3_B, coeff_nh3_A_i, coeff_nh3_B_i,
        t_0_nh3, t_00_nh3,
        ep_nh3, lv_nh3, ls_nh3, cp_nh3, r_nh3,
        C_nh3, L0_nh3, R_nh3, del_alf_nh3, del_bet_nh3, m_nh3,
        nh3, nh3_cloud, nh3_ice);

    AtomUtils::damp_wiggles(nh3,       nullptr, true, true, true);
    AtomUtils::damp_wiggles(nh3_cloud, nullptr, true, true, true);
    AtomUtils::damp_wiggles(nh3_ice,   nullptr, true, true, true);

    ChemistrySat(*this).DiffMassFluxSat();

    AtomUtils::damp_wiggles(difflux_h2s,     nullptr, true, true, true);
    AtomUtils::damp_wiggles(difflux_nh3,     nullptr, true, true, true);
    AtomUtils::damp_wiggles(difflux_nh4sh,   nullptr, true, true, true);
    AtomUtils::damp_wiggles(thermalmassflux, nullptr, true, true, true);

    ChemistrySat(*this).ChemMassRateSat();
    ChemistrySat(*this).FluxLimiterNH4SH();

    AtomUtils::damp_wiggles(massflux_h2s,   nullptr, true, true, true);
    AtomUtils::damp_wiggles(massflux_nh3,   nullptr, true, true, true);
    AtomUtils::damp_wiggles(massflux_nh4sh, nullptr, true, true, true);
    AtomUtils::damp_wiggles(fluxlim_nh4sh,  nullptr, true, true, true);

    Forces();
    Latent_Heat();

    BC_Sat(*this).bcTheta();
    BC_Sat(*this).bcPhi();
    BC_Sat(*this).bcRadius();

    restoreVar(1.0);

//    goto Printout;

    if(turb_active) TurbulenceSat(*this).init();

    // ---- Restart, if one was asked for and can be read ----
    // The loop resumes at restart_from_iter+1, so the PARITY of iter_n is preserved and the
    // physics block (which runs on even iterations) lands where it would have in one long run.
    // A missing or mismatched file is not fatal: load_state says so and the run starts from
    // scratch, because that is almost always what one wants at the end of a long queue.
    int iter_start = 1;
    if(restart_from_iter >= 0 && load_state(restart_from_iter)){
        iter_start = restart_from_iter + 1;
    }

    for(iter_n = iter_start; iter_n <= nm; iter_n++){

        auto begin = std::chrono::high_resolution_clock::now();

        cout << endl << endl;
        cout << " >>>>>>>>>>>>>>>>>>>>>>>>>>>>>    3D    <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << endl;
        cout << " 3D Saturn iterational process" << endl;
        cout << " present state of the computation " << endl << endl
             << " ======================== iteration number n = " << iter_n << " ========================= " << endl << endl
             << "    max total iteration number nm = " << nm << endl
             << "    checkpoint when to write 3D-panorama = " << checkpoint << endl
             << "    panorama_print = " << panorama_print << endl << endl;

        if(iter_n % 2 == 0){

            // One density for the whole physics block, before anything reads it. ATJUP does
            // the same at the same point; everything downstream either goes through rho_at()
            // (gated by ATSAT_LOCAL_RHO) or reads rho_mix directly when it must be local.
            computeMixtureDensity();

            getPressureSolver().run();
            AtomUtils::damp_wiggles(p_dyn, nullptr, true, true, true);

            SaturationAdjustmentSat(*this).run("CH4",
                coeff_ch4_A, coeff_ch4_B, coeff_ch4_A_i, coeff_ch4_B_i,
                t_0_ch4, t_00_ch4,
                ep_ch4, lv_ch4, ls_ch4, cp_ch4, r_ch4,
                C_ch4, L0_ch4, R_ch4, del_alf_ch4, del_bet_ch4, m_ch4,
                ch4, ch4_cloud, ch4_ice);

            SaturationAdjustmentSat(*this).run("H2O",
                coeff_h2o_A, coeff_h2o_B, coeff_h2o_A_i, coeff_h2o_B_i,
                t_0_h2o, t_00_h2o,
                ep_h2o, lv_h2o, ls_h2o, cp_h2o, r_h2o,
                C_h2o, L0_h2o, R_h2o, del_alf_h2o, del_bet_h2o, m_h2o,
                h2o, h2o_cloud, h2o_ice);

            SaturationAdjustmentSat(*this).run("NH3",
                coeff_nh3_A, coeff_nh3_B, coeff_nh3_A_i, coeff_nh3_B_i,
                t_0_nh3, t_00_nh3,
                ep_nh3, lv_nh3, ls_nh3, cp_nh3, r_nh3,
                C_nh3, L0_nh3, R_nh3, del_alf_nh3, del_bet_nh3, m_nh3,
                nh3, nh3_cloud, nh3_ice);

            // After the saturation adjustments, whose condensate it consumes; before the
            // chemistry, as in ATJUP.
            if(precip_enabled()) PrecipitationSat(*this).run();

            ChemistrySat(*this).DiffMassFluxSat();

            // Smooth difflux and thermalmassflux BEFORE difflux enters massflux assembly.
            AtomUtils::damp_wiggles(difflux_h2s,     nullptr, true, true, true);
            AtomUtils::damp_wiggles(difflux_nh3,     nullptr, true, true, true);
            AtomUtils::damp_wiggles(difflux_nh4sh,   nullptr, true, true, true);
            AtomUtils::damp_wiggles(thermalmassflux, nullptr, true, true, true);

            ChemistrySat(*this).ChemMassRateSat();
            ChemistrySat(*this).FluxLimiterNH4SH();

            // Smooth the assembled massflux.
            AtomUtils::damp_wiggles(massflux_h2s,   nullptr, true, true, true);
            AtomUtils::damp_wiggles(massflux_nh3,   nullptr, true, true, true);
            AtomUtils::damp_wiggles(massflux_nh4sh, nullptr, true, true, true);
            AtomUtils::damp_wiggles(fluxlim_nh4sh,  nullptr, true, true, true);

            Forces();
            Latent_Heat();

            if(radiation_enabled()) RadiationSat(*this).run();
            if(tw_diag_enabled())   ThermalWindDiagSat(*this).run();

        }  // if iter_n % 2

        RungeKuttaSat();

        BC_Sat(*this).bcTheta();
        BC_Sat(*this).bcPhi();
        BC_Sat(*this).bcRadius();

        // Floor the species at zero before the n-level copies are refreshed, so the clamped
        // values are what the next Runge-Kutta step starts from. After the boundary conditions,
        // because the (4/3,-1/3) extrapolation at the radial planes is one of the two sources of
        // the undershoot. See FileIO_Sat.cpp for the measurement.
        clampNegativeSpecies();

        restoreVar(1.0);

        // After the state has been advanced and the boundaries applied: put any superadiabatic
        // column back on the dry adiabat. Off by default (ATSAT_CONV_ADJ).
        // Reads the velocity field left by RK4 and the BCs, so it runs after them.
        if(turb_active) TurbulenceSat(*this).run();

        if(conv_adj_enabled()) ConvectiveAdjustmentSat(*this).run();

        panorama_cnt++;

        if(iter_n % checkpoint == 0){
            printMinMax();
            writeData();
        }

        // Shapiro filter on the velocities, opt-in. ATJUP applies one at initialisation and
        // optionally n passes per iteration; ATSAT already filters p_dyn and the mass fluxes
        // this way but has never filtered u, v, w. Grid-scale checkerboard in the velocity is
        // what the pressure projection cannot see and what feeds the polar cells.
        if(shapiro_vel_inloop() > 0)
            for(int n = 0; n < shapiro_vel_inloop(); n++) dampVelocities();

        // Non-finite census, opt-in. Runs at the checkpoint cadence so it costs nothing on the
        // iterations in between.
        if(nancheck_on() && iter_n % checkpoint == 0) nan_watch(iter_n);

        // ---- Binary restart checkpoints ----
        // One explicit dump at checkpoint_save_iter, plus a periodic one every
        // restart_save_stride iterations, which is ATJUP's arrangement and its stride of 100.
        // The periodic dump is written ONLY when the state is clean, so a diverged run can never
        // overwrite a good restart point — being able to resume from it is the file's whole
        // value. Written at the END of the iteration, after restoreVar has run, so the stored
        // n-level copies are consistent with the fields they were built from.
        if(checkpoint_save_iter >= 0 && iter_n == checkpoint_save_iter)
            save_state(iter_n);

        {
            constexpr int restart_save_stride = 100;
            if(restart_save_stride > 0 && iter_n > 0 && iter_n % restart_save_stride == 0
               && iter_n != checkpoint_save_iter){
                if(restart_state_is_clean())
                    save_state(iter_n);
                else
                    cout << "      ATSAT: restart checkpoint SKIPPED at iter " << iter_n
                         << " - non-finite cell present (state not clean)" << endl;
            }
        }

        if(panorama_cnt == panorama_print) panorama_cnt = 1;

        auto end = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for one time step\n", elapsed.count() * 1e-9);

    } // end iter_n

    cout << endl << "      Saturn: run_3D_loop atm ended ..........................." << endl;

    print_final_msg();

    return;
}
/*
*
*/
void cSaturnModel::resetArrays(){
    cout << endl << "      ATSAT: resetArrays" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    Topography.initArray_2D(jm, km, 0.0); // topography
    LatentHeat.initArray_2D(jm, km, 0.0);            // areas of higher latent heat
    precip_srf_h2o.initArray_2D(jm, km, 0.0);
    precip_srf_nh3.initArray_2D(jm, km, 0.0);
    precip_srf_ch4.initArray_2D(jm, km, 0.0);
    precip_srf_nh4sh.initArray_2D(jm, km, 0.0);
    precip_srf_total.initArray_2D(jm, km, 0.0);
    vel_star.initArray_2D(jm, km, 0.0);
    Precipitation.initArray_2D(jm, km, 0.0);         // areas of higher precipitation
    precipitable_water.initArray_2D(jm, km, 0.0);    // areas of precipitable water in the air

    aux_2D_v.initArray_2D(jm, km, 0.0);              // auxilliar field v
    aux_2D_w.initArray_2D(jm, km, 0.0);              // auxilliar field w
    tropopause_height.initArray_2D(jm, km, 0.0); // local height of the tropopause

    t.initArray(im, jm, km, ta);                    // temperature
    u.initArray(im, jm, km, ua);                    // u-component velocity component in r-direction
    v.initArray(im, jm, km, va);                    // v-component velocity component in theta-direction
    w.initArray(im, jm, km, wa);                    // w-component velocity component in phi-direction

    ch4.initArray(im, jm, km, 0.0);                 // ch4-vapour
    ch4_cloud.initArray(im, jm, km, 0.0);           // ch4-cloud
    ch4_ice.initArray(im, jm, km, 0.0);             // ch4-ice

    h2o.initArray(im, jm, km, 0.0);                    // water vapour
    h2o_cloud.initArray(im, jm, km, 0.0);                // cloud water
    h2o_ice.initArray(im, jm, km, 0.0);                  // cloud ice

    h2s.initArray(im, jm, km, 0.0);                    // water vapour

    nh3.initArray(im, jm, km, 0.0);                 // nh3-vapour
    nh3_cloud.initArray(im, jm, km, 0.0);           // nh3-cloud
    nh3_ice.initArray(im, jm, km, 0.0);             // nh3-ice

    nh4sh.initArray(im, jm, km, 0.0);                 // nh4sh

    tn.initArray(im, jm, km, ta);                    // temperature new
    un.initArray(im, jm, km, ua);                    // u-velocity component in r-direction new
    vn.initArray(im, jm, km, va);                    // v-velocity component in theta-direction new
    wn.initArray(im, jm, km, wa);                    // w-velocity component in phi-direction new

    ch4n.initArray(im, jm, km, 0.0);                // ch4 new
    ch4_cloudn.initArray(im, jm, km, 0.0);          // ch4_cloud new
    ch4_icen.initArray(im, jm, km, 0.0);            // ch4_ice new

    h2on.initArray(im, jm, km, 0.0);                    // water vapour new
    h2o_cloudn.initArray(im, jm, km, 0.0);                // cloud water new
    h2o_icen.initArray(im, jm, km, 0.0);                    // cloud ice new

    h2sn.initArray(im, jm, km, 0.0);                    // water vapour new

    nh3n.initArray(im, jm, km, 0.0);                // nh3 new
    nh3_cloudn.initArray(im, jm, km, 0.0);            // nh3_cloud new
    nh3_icen.initArray(im, jm, km, 0.0);            // nh3_ice new

    nh4shn.initArray(im, jm, km, 0.0);                // nh4sh new

    massflux_h2s.initArray(im, jm, km, 0.0);   // mass flux h2s
    massflux_nh3.initArray(im, jm, km, 0.0);   // mass flux nh3
    massflux_nh4sh.initArray(im, jm, km, 0.0);   // mass flux nh4sh

    difflux_h2s.initArray(im, jm, km, 0.0);   // diffusive flux h2s
    difflux_nh3.initArray(im, jm, km, 0.0);   // diffusive flux nh3
    difflux_nh4sh.initArray(im, jm, km, 0.0);   // diffusive flux nh4sh

    fluxlim_nh4sh.initArray(im, jm, km, 0.0);   // TVD flux-limiter correction for nh4sh advection

    thermalmassflux.initArray(im, jm, km, 0.0);   // thermal massflux_h2s

    p_dyn.initArray(im, jm, km, pa);                // dynamic pressure
    p_stat.initArray(im, jm, km, 1.0);                // static pressure
    rho_mix.initArray(im, jm, km, 0.0);               // mixture density, computeMixtureDensity()
    radiation.initArray(im, jm, km, 0.0);             // net thermal radiative flux [W/m2]
    epsilon.initArray(im, jm, km, 0.0);               // layer emissivity
    Q_rad.initArray(im, jm, km, 0.0);                 // radiative heating rate [W/m3]
    P_rain.initArray(im, jm, km, 0.0);
    P_snow.initArray(im, jm, km, 0.0);
    P_graupel.initArray(im, jm, km, 0.0);
    P_nh3_rain.initArray(im, jm, km, 0.0);
    P_nh3_snow.initArray(im, jm, km, 0.0);
    P_nh3_graupel.initArray(im, jm, km, 0.0);
    P_ch4_rain.initArray(im, jm, km, 0.0);
    P_ch4_snow.initArray(im, jm, km, 0.0);
    P_ch4_graupel.initArray(im, jm, km, 0.0);
    P_nh4sh.initArray(im, jm, km, 0.0);
    Q_precip.initArray(im, jm, km, 0.0);
    tke.initArray(im, jm, km, 0.0);
    dis.initArray(im, jm, km, 0.0);
    tken.initArray(im, jm, km, 0.0);
    disn.initArray(im, jm, km, 0.0);
    rhs_tke.initArray(im, jm, km, 0.0);
    rhs_dis.initArray(im, jm, km, 0.0);
    nue.initArray(im, jm, km, 0.0);
    nue_t.initArray(im, jm, km, 0.0);
    prod.initArray(im, jm, km, 0.0);
    tke_source.initArray(im, jm, km, 0.0);
    dis_source.initArray(im, jm, km, 0.0);
    wall_nue.initArray(im, jm, km, 0.0);
//    rho.initArray(im, jm, km, 1.0);                // density

    rhs_t.initArray(im, jm, km, 0.0);                // auxilliar field RHS temperature
    rhs_u.initArray(im, jm, km, 0.0);                // auxilliar field RHS u-velocity component
    rhs_v.initArray(im, jm, km, 0.0);                // auxilliar field RHS v-velocity component
    rhs_w.initArray(im, jm, km, 0.0);                // auxilliar field RHS w-velocity component

    rhs_ch4.initArray(im, jm, km, 0.0);                // auxilliar field RHS ch4
    rhs_ch4_cloud.initArray(im, jm, km, 0.0);          // auxilliar field RHS ch4_cloud
    rhs_ch4_ice.initArray(im, jm, km, 0.0);            // auxilliar field RHS ch4_ice

    rhs_h2o.initArray(im, jm, km, 0.0);                // auxilliar field RHS water vapour
    rhs_h2o_cloud.initArray(im, jm, km, 0.0);            // auxilliar field RHS cloud water
    rhs_h2o_ice.initArray(im, jm, km, 0.0);                // auxilliar field RHS cloud ice

    rhs_h2s.initArray(im, jm, km, 0.0);                // auxilliar field RHS water vapour

    rhs_nh3.initArray(im, jm, km, 0.0);                // auxilliar field RHS nh3
    rhs_nh3_cloud.initArray(im, jm, km, 0.0);        // auxilliar field RHS nh3_cloud
    rhs_nh3_ice.initArray(im, jm, km, 0.0);            // auxilliar field RHS nh3_ice

    rhs_nh4sh.initArray(im, jm, km, 0.0);                // auxilliar field RHS nh4sh

    aux.initArray(im, jm, km, 0.0);                // auxilliar field u-velocity component
    aux_u.initArray(im, jm, km, 0.0);                // auxilliar field u-velocity component
    aux_v.initArray(im, jm, km, 0.0);                // auxilliar field v-velocity component
    aux_w.initArray(im, jm, km, 0.0);                // auxilliar field w-velocity component

    Q_Latent.initArray(im, jm, km, 0.0);                // latent heat
    Q_Sensible.initArray(im, jm, km, 0.0);            // sensible heat
    CoriolisForce.initArray(im, jm, km, 0.0);        // Coriolis force
    CentrifugalForce.initArray(im, jm, km, 0.0);             // centrifugal force
    BuoyancyForce.initArray(im, jm, km, 0.0);        // buoyancy force, Boussinesque approximation
    PresGradForce.initArray(im, jm, km, 0.0);// pressure gradient force
    SeaMount.initArray(im, jm, km, 0.0);             // sea mount contour

    cloudiness_h2o.initArray(im, jm, km, 0.0); // cloudiness, N in literature
    cloudiness_nh3.initArray(im, jm, km, 0.0); // cloudiness, N in literature

    w_nh3.initArray(im, jm, km, 0.0);                // chemical reaction rate of nh3
    w_h2s.initArray(im, jm, km, 0.0);                // chemical reaction rate of h2s
    w_nh4sh.initArray(im, jm, km, 0.0);                // chemical reaction rate of nh4sh

    j_nh3.initArray(im, jm, km, 0.0);                // ordinary-diffusion mass flux of nh3
    j_h2s.initArray(im, jm, km, 0.0);                // ordinary-diffusion mass flux of h2s
    j_nh4sh.initArray(im, jm, km, 0.0);                // ordinary-diffusion mass flux of nh4sh

    jT_nh3.initArray(im, jm, km, 0.0);                // thermo-diffusion mass flux of nh3
    jT_h2s.initArray(im, jm, km, 0.0);                // thermo-diffusion mass flux of h2s
    jT_nh4sh.initArray(im, jm, km, 0.0);                // thermo-diffusion mass flux of nh4sh

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for resetArrays\n", elapsed.count() * 1e-9);

    cout << "      ATSAT: resetArrays ended" << endl;
    return;
}
/*
*
*/
void cSaturnModel::restoreVar(double coeff){
//    cout << endl << "      ATSAT: restoreVar" << endl;

//    auto begin = std::chrono::high_resolution_clock::now();

    #pragma omp parallel for
    for(int i = 0; i < im; i++){
        for(int j = 0; j < jm; j++){
            for(int k = 0; k < km; k++){
                tn.x[i][j][k] = coeff * t.x[i][j][k];
                // k* and dis* keep the same start-of-step copies the other prognostic fields get;
                // without them the RK4 stages would integrate from a moving base. coeff is the
                // restore factor the other fields use, so they stay in step with it.
                tken.x[i][j][k] = coeff * tke.x[i][j][k];
                disn.x[i][j][k] = coeff * dis.x[i][j][k];
                un.x[i][j][k] = coeff * u.x[i][j][k];
                vn.x[i][j][k] = coeff * v.x[i][j][k];
                wn.x[i][j][k] = coeff * w.x[i][j][k];


//                u.x[i][j][k] = un.x[i][j][k] = 0.0;
//                v.x[i][j][k] = vn.x[i][j][k] = 0.0;


                ch4n.x[i][j][k] = coeff * ch4.x[i][j][k];
                ch4_cloudn.x[i][j][k] = coeff * ch4_cloud.x[i][j][k];
                ch4_icen.x[i][j][k] = coeff * ch4_ice.x[i][j][k];
                h2on.x[i][j][k] = coeff * h2o.x[i][j][k];
                h2o_cloudn.x[i][j][k] = coeff * h2o_cloud.x[i][j][k];
                h2o_icen.x[i][j][k] = coeff * h2o_ice.x[i][j][k];
                h2sn.x[i][j][k] = coeff * h2s.x[i][j][k];
                nh3n.x[i][j][k] = coeff * nh3.x[i][j][k];
                nh3_cloudn.x[i][j][k] = coeff * nh3_cloud.x[i][j][k];
                nh3_icen.x[i][j][k] = coeff * nh3_ice.x[i][j][k];
                nh4shn.x[i][j][k] = coeff * nh4sh.x[i][j][k];
            }
        }
    }

//    auto end = std::chrono::high_resolution_clock::now();
//    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
//    printf(" time measured: %.3f seconds for restoreVar\n", elapsed.count() * 1e-9);

//    cout << "      ATSAT: restoreVar ended" << endl;

    return;
}
/*
*
*/
/*
    cout << endl << "      ATSAT: thermodynamic fft_gaussian_filter_3d in Run begin ......................." << endl;

    fft_gaussian_filter_3d(nh3,1);
    fft_gaussian_filter_3d(nh3_cloud,1);
    fft_gaussian_filter_3d(nh3_ice,1);

    fft_gaussian_filter_3d(nh4sh,1);
//    fft_gaussian_filter_3d(nh4sh_cloud,1);
//    fft_gaussian_filter_3d(nh4sh_ice,1);

    fft_gaussian_filter_3d(h2s,1);
    fft_gaussian_filter_3d(h2s_cloud,1);
    fft_gaussian_filter_3d(h2s_ice,1);

    fft_gaussian_filter_3d(h2o,1);
    fft_gaussian_filter_3d(h2o_cloud,1);
    fft_gaussian_filter_3d(h2o_ice,1);

    cout << endl << "      ATSAT: thermodynamic fft_gaussian_filter_3d in Run end ......................." << endl;
*/


