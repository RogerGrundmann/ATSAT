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
#include "TbSat.h"
#include "CwbSat.h"         // ATSAT_CWB_DIAG, print-only
#include "UbSat.h"           // ATSAT_UB_DIAG, print-only
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


// Grey multi-layer radiation (RadiationSat), ported from ATJUP. DEFAULT ON since 2026-10-11 on the
// user's word (ATSAT_RADIATION=0 restores the run without it; not byte-checked). It fills the DIAGNOSTIC arrays radiation / epsilon / Q_rad and touches neither t nor any
// rhs — wiring the heating into the temperature equation is a separate step. Runs on even
// iterations, alongside the rest of the physics block.
static int radiation_enabled(){
    static const int v = [](){ const char* e = getenv("ATSAT_RADIATION"); return e ? atoi(e) : 1; }();
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
// DEFAULT 0.001 SINCE 2026-10-09 (ATJUP's value; 1.06 s of Saturn time per iteration, 224
// iterations = 4 minutes). It needed two other defaults to move with it, both in cSaturnModel.h:
// the metric radius (the unit sphere's polar spacing was what limited the step) and the pole
// copy (the extrapolated pole cooled without bound once the step was long enough to see it).
// MEASURED with both, 224 iterations, 8 threads (ATJUP/satchk/giants/G1, G5, G25): dt = 0.001
// and 0.005 run clean, min T flat at -243.1 degC, max|w| 267 -> 259 and 263 -> 251 m/s;
// dt = 0.025 goes non-finite after iteration 160. ATSAT_DT=<value> overrides; ATSAT_DT=-1
// gives the old formula.
static double timestep_override(){
    static const double v = [](){ const char* e = getenv("ATSAT_DT"); return e ? atof(e) : 0.001; }();
    return v;
}


// Precipitation microphysics (PrecipitationSat), ported from ATJUP. DEFAULT ON since 2026-07-31,
// by decision, with the caveat below standing. It DOES feed back: the condensate it converts is
// removed from the cloud/ice fields in place, so it must run AFTER the SaturationAdjustmentSat
// calls or the adjustment would simply undo it. ATSAT_PRECIP=0 restores the old behaviour.
//
// THE CAVEAT, measured 2026-07-31 over 50 iterations (config_m50, 12 threads), and not repaired
// by switching the knob: the scheme runs and produces fluxes, but where it is strongest they are
// the SAFETY CAP and not physics. At the i=20 level (250 km), reading the .vtk directly:
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
// So where P_snow reads 8.6400 the number is the clamp speaking, not the microphysics, and the
// recalibration is still owed: the observable is the one ATJUP used, Lv*P against Saturn's
// emitted flux. Until then, read the fluxes as an upper bound wherever they sit at 8.6400 and as
// physics wherever they do not.
//
// RECALIBRATED 2026-10-09: the five rate coefficients carry Planet::precip_rate_scale() = 0.006
// on Saturn (cSaturnModel.h has the measurement). No cell sits on the cap any more; the caveat
// above describes the runs before that.
static int precip_enabled(){
    static const int v = [](){ const char* e = getenv("ATSAT_PRECIP"); return e ? atoi(e) : 1; }();
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
// ATSAT_TURB: -1 = not set. turb_model IS the switch since 2026-10-11 (the user's word), as on ATJUP
// since 2026-10-09: the closure runs unless turb_model (or ATSAT_TURB_MODEL, which overrides it) is
// "none" / "laminar" / empty, and param.py's default is k_omega_SST, so a stock run has it ON.
// Until then the closure needed ATSAT_TURB=1 as well and the configured model alone did nothing.
// ATSAT_TURB=0 still switches it off; any other value leaves the decision to turb_model.
static int turb_env_enabled(){
    static const int v = [](){ const char* e = getenv("ATSAT_TURB"); return e ? atoi(e) : -1; }();
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
// ATSAT_INIT_SPECIES_LOCAL_RHO: see the note at rho_init() in cSaturnModel.h.
void cSaturnModel::fill_rho_mix_init(){
    if(!init_species_local_rho()) return;
    double rmin = 1.0e30, rmax = 0.0;
    for(int i = 0; i < im; i++)
        for(int j = 0; j < jm; j++)
            for(int k = 0; k < km; k++){
                const double r = rho_init(i, j, k);
                rho_mix.x[i][j][k] = r;
                rmin = std::min(rmin, r); rmax = std::max(rmax, r);
            }
    printf("      %s: INIT_SPECIES_LOCAL_RHO: R_mix = %.6g, r_mix = %.6g, set-up density %.6g .. %.6g kg/m3\n",
           "ATSAT", R_mix, r_mix, rmin, rmax);
}

void cSaturnModel::damp_species_init(Array& f){
    if(!init_species_local_rho()){ AtomUtils::damp_wiggles(f, nullptr, true, true, true); return; }
    for(int i = 0; i < im; i++)
        for(int j = 0; j < jm; j++)
            for(int k = 0; k < km; k++) f.x[i][j][k] /= rho_init(i, j, k);
    AtomUtils::damp_wiggles(f, nullptr, true, true, true);
    for(int i = 0; i < im; i++)
        for(int j = 0; j < jm; j++)
            for(int k = 0; k < km; k++) f.x[i][j][k] *= rho_init(i, j, k);
}

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

    // ATSAT_BC_PHI_PERIODIC (2026-10-13; DEFAULT 1 on the user's word, 0 restores byte-identically;
    // ATPhys.h has the description): the seam
    // meridian k = 0 = km-1 is integrated instead of being reconstructed from its neighbours.
    // The knob is read by every loop that differences in phi; the Shapiro filter lives in lib/
    // and is told here.
    AtomUtils::set_phi_seam_periodic(ATPhys::phi_periodic<cSaturnModel>());
    if(ATPhys::phi_periodic<cSaturnModel>())
        printf("      ATSAT: the phi seam is periodic: k = 0 is integrated, k = km-1 is its copy (ATSAT_BC_PHI_PERIODIC=0 restores the reconstructed seam)\n");
    cout.precision(6);
    cout.setf(ios::fixed);
    if(!has_welcome_msg_printed)
        print_welcome_msg();
    initMsg();

    iter_n = 0;
    panorama_cnt = 0;

    resetArrays();

    // ATSAT_NH3_ALF overrides the exponent of the NH3 saturation curves, liquid and ice together
    // (cSaturnModel.h has the note); -1.2 is what every run before 2026-10-09 was made with.
    if(const char* e = getenv("ATSAT_NH3_ALF")){ del_alf_nh3 = atof(e); del_alf_nh3_ice = atof(e); }

    dt = 2.8284 * dr/u_0 * 0.2;
    if(timestep_override() > 0.0) dt = timestep_override();
    // Resolve the turbulence switch. Order matters: the configuration value is already loaded,
    // ATSAT_TURB_MODEL overrides it here, and only then is turb_active fixed — everything
    // downstream (TurbulenceSat, RHSSat, RungeKuttaSat, the BC routines) reads that flag and never
    // the string. "none" now means what it says; before, it fell through to k-omega SST.
    if(const char* tm = getenv("ATSAT_TURB_MODEL")) turb_model = tm;
    turb_active = (turb_env_enabled() != 0) && !(turb_model.empty() || turb_model == "none" || turb_model == "laminar");
    if(turb_env_enabled() > 0 && !turb_active)
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
    fill_rho_mix_init();
    init_PressureDynamic();

    init_vapour_cloud_ice("CH4", ch4_tropopause, coeff_ch4_A, coeff_ch4_B,
        coeff_ch4_A_i, coeff_ch4_B_i, t_0_ch4, t_00_ch4,
        ep_ch4, r_ch4, m_ch4,
        C_ch4, L0_ch4, R_ch4, del_alf_ch4, del_bet_ch4,
        ch4, ch4_cloud, ch4_ice, cloudiness_h2o);

    damp_species_init(ch4);
    damp_species_init(ch4_cloud);
    damp_species_init(ch4_ice);

    init_vapour_cloud_ice("H2O", h2o_tropopause, coeff_h2o_A, coeff_h2o_B,
        coeff_h2o_A_i, coeff_h2o_B_i, t_0_h2o, t_00_h2o,
        ep_h2o, r_h2o, m_h2o,
        C_h2o, L0_h2o, R_h2o, del_alf_h2o, del_bet_h2o,
        h2o, h2o_cloud, h2o_ice, cloudiness_h2o);

    damp_species_init(h2o);
    damp_species_init(h2o_cloud);
    damp_species_init(h2o_ice);

    init_h2s("H2S", h2s_tropopause, coeff_h2s_A, coeff_h2s_B,
        coeff_h2s_A_i, coeff_h2s_B_i, t_0_h2s, t_00_h2s,
        ep_h2s, r_h2s, m_h2s,
        C_h2s, L0_h2s, R_h2s, del_alf_h2s, del_bet_h2s, h2s);

    damp_species_init(h2s);

    init_vapour_cloud_ice("NH3", nh3_tropopause, coeff_nh3_A, coeff_nh3_B,
        coeff_nh3_A_i, coeff_nh3_B_i, t_0_nh3, t_00_nh3,
        ep_nh3, r_nh3, m_nh3,
        C_nh3, L0_nh3, R_nh3, del_alf_nh3, del_bet_nh3,
        nh3, nh3_cloud, nh3_ice, cloudiness_nh3);

    damp_species_init(nh3);
    damp_species_init(nh3_cloud);
    damp_species_init(nh3_ice);

    init_nh4sh("NH4SH", nh4sh_tropopause, coeff_nh4sh_A, coeff_nh4sh_B,
        ep_nh4sh, r_nh4sh, m_nh4sh,
        C_nh4sh, L0_nh4sh, R_nh4sh,
        del_alf_nh4sh, del_bet_nh4sh, nh4sh);

    damp_species_init(nh4sh);

//    goto Printout;

    SaturationAdjustmentSat(*this).run("CH4",
        coeff_ch4_A, coeff_ch4_B, coeff_ch4_A_i, coeff_ch4_B_i,
        t_0_ch4, t_00_ch4,
        ep_ch4, lv_ch4, ls_ch4, cp_ch4, r_ch4,
        C_ch4, L0_ch4, R_ch4, del_alf_ch4, del_bet_ch4, m_ch4,
        ch4, ch4_cloud, ch4_ice);

    damp_species_init(ch4);
    damp_species_init(ch4_cloud);
    damp_species_init(ch4_ice);

    SaturationAdjustmentSat(*this).run("H2O",
        coeff_h2o_A, coeff_h2o_B, coeff_h2o_A_i, coeff_h2o_B_i,
        t_0_h2o, t_00_h2o,
        ep_h2o, lv_h2o, ls_h2o, cp_h2o, r_h2o,
        C_h2o, L0_h2o, R_h2o, del_alf_h2o, del_bet_h2o, m_h2o,
        h2o, h2o_cloud, h2o_ice);

    damp_species_init(h2o);
    damp_species_init(h2o_cloud);
    damp_species_init(h2o_ice);

    SaturationAdjustmentSat(*this).run("NH3",
        coeff_nh3_A, coeff_nh3_B, coeff_nh3_A_i, coeff_nh3_B_i,
        t_0_nh3, t_00_nh3,
        ep_nh3, lv_nh3, ls_nh3, cp_nh3, r_nh3,
        C_nh3, L0_nh3, R_nh3, del_alf_nh3, del_bet_nh3, m_nh3,
        nh3, nh3_cloud, nh3_ice);

    damp_species_init(nh3);
    damp_species_init(nh3_cloud);
    damp_species_init(nh3_ice);

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

    // ATSAT_GEO_BALANCE (cSaturnModel.h): the balanced start, after the boundary passes. Its
    // solver half is the shared PressureSolver's ATSAT_PRESS_INCR, whose default is 0 there; the
    // balanced start sets it unless the knob is given.
    if(geo_balance()) setenv("ATSAT_PRESS_INCR", "1", 0);
    initCellHold();                                                     // ATSAT_CELL_HOLD, cSaturnModel.h
    initGeoBalance();

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

    // ATSAT_TB_DIAG=<n> (print-only, see TbSat.h): the level-mean temperature after every stage of
    // the loop that can write t, charged to that stage, in K; printed every n iterations.
    static const int tb_every = [](){ const char* e = getenv("ATSAT_TB_DIAG"); return e ? atoi(e) : 0; }();
    const bool tb_on = tb_every > 0 && im <= SatTb::MAXLEV;
    SatTb::st().on = tb_on;
    enum { TB_PRE = 0, TB_SATADJ, TB_PRECIP, TB_PHYS, TB_RK_RESET, TB_RK_REST, TB_BC, TB_CLAMP, TB_AFTER, TB_N };
    static const char* tb_name[TB_N] = { "pressure", "sat.adjust", "precip", "chem/latent", "RK reset", "RK rest",
        "BC", "clamp", "turb/conv" };
    static const char* tb_term[SatTb::NTERM] = { "RK press.work", "RK transp r", "RK transp the", "RK transp phi",
        "RK diffusion", "RK th.massflux", "RK radiation", "RK precip" };
    std::vector<std::vector<double>> tb_sum(TB_N, std::vector<double>(im, 0.0)), tb_tsum(SatTb::NTERM, std::vector<double>(im, 0.0));
    std::vector<double> tb_last(im, 0.0), tb_first(im, 0.0);
    double tb_norm = 0.0;
    auto tb_lev = [&](Array& a) -> std::vector<double> {
        std::vector<double> m_(im, 0.0);
        #pragma omp parallel for schedule(static)
        for(int i = 0; i < im; i++){
            double s_ = 0.0;
            for(int j = 0; j < jm; j++){
                const double sj = std::sin(the.z[j]);
                for(int k = 0; k < km; k++) s_ += sj * a.x[i][j][k];
            }
            m_[i] = s_ * t_ref / tb_norm;
        }
        return m_; };
    auto tb_mark = [&](int stage){
        if(!tb_on) return;
        const std::vector<double> w_ = tb_lev(t);
        for(int i = 0; i < im; i++){ tb_sum[stage][i] += w_[i] - tb_last[i]; tb_last[i] = w_[i]; } };
    if(tb_on){
        for(int j = 0; j < jm; j++) tb_norm += std::sin(the.z[j]) * km;
        tb_last = tb_first = tb_lev(t);
    }

    // ATSAT_UB_DIAG=<n> (print-only, see UbSat.h): the level mean of the radial wind after every
    // stage of the loop that can write u, charged to that stage, in m/s; the Runge-Kutta step split
    // into the terms of rhs_u; and the radial pressure force with what the pressure solver and the
    // filter on p_dyn did to it. For the whole level and for the equatorial band (5 S .. 5 N by
    // default, ATSAT_UB_BAND=<degrees>). Printed every n iterations.
    static const int ub_every = [](){ const char* e = getenv("ATSAT_UB_DIAG"); return e ? atoi(e) : 0; }();
    static const int ub_band  = [](){ const char* e = getenv("ATSAT_UB_BAND"); return e ? atoi(e) : 5; }();
    const bool ub_on = ub_every > 0 && im <= SatUb::MAXLEV;
    SatUb::st().on = ub_on;
    SatUb::st().j_lo = (jm - 1) / 2 - ub_band;
    SatUb::st().j_hi = (jm - 1) / 2 + ub_band;
    enum { UB_RK_RESET = 0, UB_RK_REST, UB_FILTER, UB_BC, UB_N };
    static const char* ub_name[UB_N] = { "RK reset", "RK rest", "filter on u", "BC" };
    static const char* ub_term[SatUb::NTERM] = { "RK -dp/dr", "RK +(v2+w2)/r", "RK -advection", "RK diffusion",
        "RK Coriolis", "RK centrifugal", "RK buoyancy", "RK sponge" };
    static const char* ub_reg[SatUb::NREG] = { "whole level", "equatorial band" };
    typedef std::vector<std::vector<double>> ub_tab;                       // [region][level]
    std::vector<ub_tab> ub_sum(UB_N, ub_tab(SatUb::NREG, std::vector<double>(im, 0.0)));
    std::vector<ub_tab> ub_tsum(SatUb::NTERM, ub_tab(SatUb::NREG, std::vector<double>(im, 0.0)));
    ub_tab ub_last(SatUb::NREG, std::vector<double>(im, 0.0)), ub_first = ub_last;
    ub_tab ub_f_first = ub_last, ub_m_first = ub_last, ub_f_solver = ub_last, ub_f_filter = ub_last;
    double ub_norm[SatUb::NREG] = { 0.0, 0.0 };
    // the cells the Runge-Kutta step integrates in phi: k = 0..km-2 with the periodic seam
    const int ub_k0 = ATPhys::phi_first(ATPhys::phi_periodic<cSaturnModel>());
    auto ub_in = [&](int r_, int j) -> bool {
        return r_ == SatUb::R_ALL || (j >= SatUb::st().j_lo && j <= SatUb::st().j_hi); };
    // level mean of a field [m/s]
    auto ub_lev = [&](Array& a) -> ub_tab {
        ub_tab m_(SatUb::NREG, std::vector<double>(im, 0.0));
        #pragma omp parallel for schedule(static)
        for(int i = 0; i < im; i++){
            double s_[SatUb::NREG] = { 0.0, 0.0 };
            for(int j = 1; j < jm-1; j++){
                const double sj = std::sin(the.z[j]);
                double row_ = 0.0;
                for(int k = ub_k0; k < km-1; k++) row_ += a.x[i][j][k];
                for(int r_ = 0; r_ < SatUb::NREG; r_++) if(ub_in(r_, j)) s_[r_] += sj * row_;
            }
            for(int r_ = 0; r_ < SatUb::NREG; r_++) m_[r_][i] = s_[r_] * u_0 / ub_norm[r_];
        }
        return m_; };
    // level mean of the radial pressure force -dp_dyn/dr and of the metric term (v2 + w2)/r as
    // rhs_u forms them, in m/s per iteration
    auto ub_force = [&](bool metric) -> ub_tab {
        ub_tab m_(SatUb::NREG, std::vector<double>(im, 0.0));
        const double scale_p_ = 1.0e5 / (r_mix * u_0 * u_0);
        #pragma omp parallel for schedule(static)
        for(int i = 1; i < im-1; i++){
            double s_[SatUb::NREG] = { 0.0, 0.0 };
            const double inv_rm_ = 1.0 / metricRadius(rad.z[i]);
            for(int j = 1; j < jm-1; j++){
                const double sj = std::sin(the.z[j]);
                double row_ = 0.0;
                for(int k = ub_k0; k < km-1; k++)
                    row_ += metric ? (v.x[i][j][k] * v.x[i][j][k] + w.x[i][j][k] * w.x[i][j][k]) * inv_rm_
                                   : -(p_dyn.x[i+1][j][k] - p_dyn.x[i-1][j][k]) / (2.0 * dr) * (coord_stretching ? inv_rm_ / (1.0 + inv_rm_) : 1.0) + 0.0 * scale_p_;
                for(int r_ = 0; r_ < SatUb::NREG; r_++) if(ub_in(r_, j)) s_[r_] += sj * row_;
            }
            for(int r_ = 0; r_ < SatUb::NREG; r_++) m_[r_][i] = s_[r_] * dt * u_0 / ub_norm[r_];
        }
        return m_; };
    auto ub_mark = [&](int stage){
        if(!ub_on) return;
        const ub_tab w_ = ub_lev(u);
        for(int r_ = 0; r_ < SatUb::NREG; r_++) for(int i = 0; i < im; i++){
            ub_sum[stage][r_][i] += w_[r_][i] - ub_last[r_][i]; ub_last[r_][i] = w_[r_][i]; } };
    if(ub_on){
        for(int j = 1; j < jm-1; j++) for(int r_ = 0; r_ < SatUb::NREG; r_++)
            if(ub_in(r_, j)) ub_norm[r_] += std::sin(the.z[j]) * (km - 1 - ub_k0);
        ub_last = ub_first = ub_lev(u);
        ub_f_first = ub_force(false);
        ub_m_first = ub_force(true);
    }

    // ATSAT_CWB_DIAG=<n> (print-only, see CwbSat.h): the H2O, NH3 and CH4 columns after every
    // stage of the loop that can write them, charged to that stage, in g/m2; printed every n
    // iterations. The stages are the ones ATSAT_TB_DIAG marks.
    static const int cwb_every = [](){ const char* e = getenv("ATSAT_CWB_DIAG"); return e ? atoi(e) : 0; }();
    const bool cwb_on = cwb_every > 0;
    SatCwb::st().on = cwb_on;
    static const char* cw_gas[SatCwb::NGAS] = { "H2O", "NH3", "CH4" };
    static const char* cw_name[TB_N] = { "before the adjustment (pressure)", "SaturationAdjustment (all three calls)",
        "Precipitation (in place)", "scalar sync, chemistry, forces, radiation", "RungeKutta: reset to the n-level copies",
        "RungeKutta: floor, limiters and remainder", "bcTheta, bcPhi, bcRadius", "clampNegativeSpecies", "turbulence, convective adjustment" };
    static const char* cw_term[SatCwb::NTERM] = { "RungeKutta: radial transport   -u dq/dr", "RungeKutta: meridional transport",
        "RungeKutta: zonal transport", "RungeKutta: diffusion", "RungeKutta: chemistry source", "RungeKutta: precipitation source",
        "RungeKutta: flux-form term  -s q div(u)", "(reference, not a term) q*div(u)" };
    Array* cw_f[SatCwb::NGAS][3]  = { { &h2o,  &h2o_cloud,  &h2o_ice  }, { &nh3,  &nh3_cloud,  &nh3_ice  }, { &ch4,  &ch4_cloud,  &ch4_ice  } };
    Array* cw_fn[SatCwb::NGAS][3] = { { &h2on, &h2o_cloudn, &h2o_icen }, { &nh3n, &nh3_cloudn, &nh3_icen }, { &ch4n, &ch4_cloudn, &ch4_icen } };
    double cw_sum[SatCwb::NGAS][TB_N] = {{0.0}}, cw_tsum[SatCwb::NGAS][SatCwb::NTERM] = {{0.0}};
    double cw_norm = 0.0, cw_last[SatCwb::NGAS] = {0.0}, cw_first[SatCwb::NGAS] = {0.0}, cw_first0[SatCwb::NGAS] = {0.0};
    // i_lo = 0: every level with a thickness; i_hi = 1: the plane i = 0 alone.
    auto cwb_col = [&](Array* f[3], int i_hi) -> double {
        double s_ = 0.0;
        const int ie = i_hi > 0 ? i_hi : im;
        #pragma omp parallel for collapse(2) schedule(static) reduction(+:s_)
        for(int k = 0; k < km; k++){
            for(int j = 0; j < jm; j++){
                const double sj = std::sin(the.z[j]);
                for(int i = 0; i < ie; i++){
                    const double dz = layer_thickness_m(i);
                    if(dz > 0.0) s_ += sj * dz * (f[0]->x[i][j][k] + f[1]->x[i][j][k] + f[2]->x[i][j][k]);
                }
            }
        }
        return s_ * 1.0e3 / cw_norm; };
    auto cwb_mark = [&](int stage){
        if(!cwb_on) return;
        for(int g_ = 0; g_ < SatCwb::NGAS; g_++){
            const double w_ = cwb_col(cw_f[g_], 0);
            cw_sum[g_][stage] += w_ - cw_last[g_]; cw_last[g_] = w_;
        } };
    if(cwb_on){
        for(int j = 0; j < jm; j++) cw_norm += std::sin(the.z[j]) * km;
        for(int g_ = 0; g_ < SatCwb::NGAS; g_++){
            cw_last[g_] = cw_first[g_] = cwb_col(cw_f[g_], 0);
            cw_first0[g_] = cwb_col(cw_f[g_], 1);
        }
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

            ub_tab ub_f0, ub_f1;
            if(ub_on) ub_f0 = ub_force(false);
            getPressureSolver().run();
            if(ub_on) ub_f1 = ub_force(false);
            // ATSAT_DAMP_PDYN=0 (default 0 with ATSAT_GEO_BALANCE, else 1): no Shapiro pass on p_dyn;
            // on ATNEPT it halved the balancing meridional gradient in 224 iterations (DAMP_P_HORIZ)
            static const bool damp_pdyn = [](){ const char* e = getenv("ATSAT_DAMP_PDYN"); return e ? atoi(e) != 0 : !geo_balance(); }();
            if(damp_pdyn) AtomUtils::damp_wiggles(p_dyn, nullptr, true, true, true);
            if(ub_on){
                const ub_tab ub_f2 = ub_force(false);
                for(int r_ = 0; r_ < SatUb::NREG; r_++) for(int i = 0; i < im; i++){
                    ub_f_solver[r_][i] += ub_f1[r_][i] - ub_f0[r_][i];
                    ub_f_filter[r_][i] += ub_f2[r_][i] - ub_f1[r_][i];
                }
            }
            tb_mark(TB_PRE); cwb_mark(TB_PRE);

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
            tb_mark(TB_SATADJ); cwb_mark(TB_SATADJ);
            if(precip_enabled()) PrecipitationSat(*this).run();
            tb_mark(TB_PRECIP); cwb_mark(TB_PRECIP);

            // ATSAT_RK_SCALAR_SYNC=<0|1|2> (2026-10-11; DEFAULT 2 on the user's word, 0 restores):
            // ATJUP's repair of 2026-10-08, ported. RungeKuttaSat integrates every scalar from
            // its n-level copy (h2on, ..., tn), and restoreVar refreshed those at the END of the
            // previous iteration -- before the saturation adjustment above ran. So everything the
            // adjustment writes directly (the condensate it forms, the latent heat it releases)
            // was overwritten by the step that follows. MEASURED with ATSAT_TB_DIAG (satchk/
            // giants/TB1, 56 iterations, level means): the adjustment +7.2 / +14.9 / +4.5 K at
            // levels 8 / 12 / 16 and +12.4 K at level 39, "RK reset" the same amounts with the
            // opposite sign, to the digit -- the adjustment re-condensed the same cells at every
            // call and its heat never reached the temperature.
            // 1 copies the nine condensable fields into their n-level copies, so the step starts
            // from the adjusted state; 2 copies the temperature too, so the latent heat goes with
            // the condensate. The in-place precipitation above is kept by it as well.
            static const int rk_scalar_sync = [](){ const char* e = getenv("ATSAT_RK_SCALAR_SYNC"); return e ? atoi(e) : 2; }();
            if(rk_scalar_sync >= 1){
                #pragma omp parallel for collapse(2) schedule(static)
                for(int i = 0; i < im; i++){
                    for(int j = 0; j < jm; j++){
                        for(int k = 0; k < km; k++){
                            h2on.x[i][j][k]       = h2o.x[i][j][k];
                            h2o_cloudn.x[i][j][k] = h2o_cloud.x[i][j][k];
                            h2o_icen.x[i][j][k]   = h2o_ice.x[i][j][k];
                            nh3n.x[i][j][k]       = nh3.x[i][j][k];
                            nh3_cloudn.x[i][j][k] = nh3_cloud.x[i][j][k];
                            nh3_icen.x[i][j][k]   = nh3_ice.x[i][j][k];
                            ch4n.x[i][j][k]       = ch4.x[i][j][k];
                            ch4_cloudn.x[i][j][k] = ch4_cloud.x[i][j][k];
                            ch4_icen.x[i][j][k]   = ch4_ice.x[i][j][k];
                            if(rk_scalar_sync >= 2) tn.x[i][j][k] = t.x[i][j][k];
                        }
                    }
                }
            }

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
            tb_mark(TB_PHYS); cwb_mark(TB_PHYS);

        }  // if iter_n % 2

        std::vector<double> tb_n;
        if(tb_on){
            tb_n = tb_lev(tn);                                          // what the step starts from
            for(int t_ = 0; t_ < SatTb::MAXTHR; t_++) for(int q_ = 0; q_ < SatTb::NTERM; q_++)
                for(int i = 0; i < im; i++) SatTb::st().acc[t_][q_][i] = 0.0;
        }

        double cw_n[SatCwb::NGAS] = {0.0};
        if(cwb_on){
            for(int g_ = 0; g_ < SatCwb::NGAS; g_++) cw_n[g_] = cwb_col(cw_fn[g_], 0);     // what the step starts from
            for(int t_ = 0; t_ < SatCwb::MAXTHR; t_++) for(int g_ = 0; g_ < SatCwb::NGAS; g_++)
                for(int q_ = 0; q_ < 8; q_++) SatCwb::st().acc[t_][g_][q_] = 0.0;
        }

        ub_tab ub_n;
        if(ub_on){
            ub_n = ub_lev(un);                                          // what the step starts from
            for(int t_ = 0; t_ < SatUb::MAXTHR; t_++) for(int r_ = 0; r_ < SatUb::NREG; r_++)
                for(int q_ = 0; q_ < SatUb::NTERM; q_++) for(int i = 0; i < im; i++) SatUb::st().acc[t_][r_][q_][i] = 0.0;
        }

        RungeKuttaSat();

        if(ub_on){
            const ub_tab w_ = ub_lev(u);
            for(int r_ = 0; r_ < SatUb::NREG; r_++) for(int i = 0; i < im; i++){
                double step_terms = 0.0;
                for(int q_ = 0; q_ < SatUb::NTERM; q_++){
                    double a_ = 0.0;
                    for(int t_ = 0; t_ < SatUb::MAXTHR; t_++) a_ += SatUb::st().acc[t_][r_][q_][i];
                    a_ *= dt / 6.0 * u_0 / ub_norm[r_];
                    ub_tsum[q_][r_][i] += a_; step_terms += a_;
                }
                ub_sum[UB_RK_RESET][r_][i] += ub_n[r_][i] - ub_last[r_][i];
                ub_sum[UB_RK_REST][r_][i]  += w_[r_][i] - ub_n[r_][i] - step_terms;
                ub_last[r_][i] = w_[r_][i];
            }
        }


        if(cwb_on){
            for(int g_ = 0; g_ < SatCwb::NGAS; g_++){
                const double w_ = cwb_col(cw_f[g_], 0);
                double step_terms = 0.0;
                for(int q_ = 0; q_ < SatCwb::NTERM; q_++){
                    double a_ = 0.0;
                    for(int t_ = 0; t_ < SatCwb::MAXTHR; t_++) a_ += SatCwb::st().acc[t_][g_][q_];
                    a_ *= dt / 6.0 * 1.0e3 / cw_norm;
                    cw_tsum[g_][q_] += a_;
                    if(q_ != SatCwb::T_QDIV) step_terms += a_;
                }
                cw_sum[g_][TB_RK_RESET] += cw_n[g_] - cw_last[g_];
                cw_sum[g_][TB_RK_REST]  += w_ - cw_n[g_] - step_terms;
                cw_last[g_] = w_;
            }
        }

        if(tb_on){
            const std::vector<double> w_ = tb_lev(t);
            for(int i = 0; i < im; i++){
                double step_terms = 0.0;
                for(int q_ = 0; q_ < SatTb::NTERM; q_++){
                    double a_ = 0.0;
                    for(int t_ = 0; t_ < SatTb::MAXTHR; t_++) a_ += SatTb::st().acc[t_][q_][i];
                    a_ *= dt / 6.0 * t_ref / tb_norm;
                    tb_tsum[q_][i] += a_; step_terms += a_;
                }
                tb_sum[TB_RK_RESET][i] += tb_n[i] - tb_last[i];
                tb_sum[TB_RK_REST][i]  += w_[i] - tb_n[i] - step_terms;
                tb_last[i] = w_[i];
            }
        }

        BC_Sat(*this).bcTheta();
        BC_Sat(*this).bcPhi();
        BC_Sat(*this).bcRadius();
        ub_mark(UB_BC);

        // Floor the species at zero before the n-level copies are refreshed, so the clamped
        // values are what the next Runge-Kutta step starts from. After the boundary conditions,
        // because the (4/3,-1/3) extrapolation at the radial planes is one of the two sources of
        // the undershoot. See FileIO_Sat.cpp for the measurement.
        tb_mark(TB_BC); cwb_mark(TB_BC);
        clampNegativeSpecies();
        tb_mark(TB_CLAMP); cwb_mark(TB_CLAMP);

        // Steady-state query: max|f - f_n| per field with its location, plus the continuity
        // residual. BEFORE restoreVar, which is what makes the differences non-zero — see the
        // note on the function. Same cadence as printMinMax; ATSAT_STEADY=0 switches it off.
        static const int steady_on = [](){
            const char* e = getenv("ATSAT_STEADY"); return e ? atoi(e) : 1; }();
        if(steady_on && iter_n % checkpoint == 0) steadyQuery();

        restoreVar(1.0);

        // After the state has been advanced and the boundaries applied: put any superadiabatic
        // column back on the dry adiabat. Off by default (ATSAT_CONV_ADJ).
        // Reads the velocity field left by RK4 and the BCs, so it runs after them.
        if(turb_active) TurbulenceSat(*this).run();

        if(conv_adj_enabled()) ConvectiveAdjustmentSat(*this).run();
        tb_mark(TB_AFTER); cwb_mark(TB_AFTER);
        if(tb_on && iter_n % tb_every == 0){
            static const int lv_[] = { 0, 1, 2, 3, 5, 8, 12, 16, 20, 24, 28, 32, 35, 37, 38, 39, 40 };
            const int nl_ = (int)(sizeof(lv_) / sizeof(lv_[0]));
            printf("\n      ATSAT: [TB] level-mean temperature, cumulative over iterations %d..%d, K\n      ATSAT: [TB] %-15s", iter_start, iter_n, "level");
            for(int l = 0; l < nl_; l++) if(lv_[l] < im) printf(" %8d", lv_[l]);
            printf("\n      ATSAT: [TB] %-15s", "T now");
            for(int l = 0; l < nl_; l++) if(lv_[l] < im) printf(" %8.2f", tb_last[lv_[l]]);
            printf("\n      ATSAT: [TB] %-15s", "change");
            for(int l = 0; l < nl_; l++) if(lv_[l] < im) printf(" %+8.3f", tb_last[lv_[l]] - tb_first[lv_[l]]);
            for(int q_ = 0; q_ < TB_N; q_++){
                printf("\n      ATSAT: [TB] %-15s", tb_name[q_]);
                for(int l = 0; l < nl_; l++) if(lv_[l] < im) printf(" %+8.3f", tb_sum[q_][lv_[l]]);
                if(q_ == TB_RK_RESET) for(int r_ = 0; r_ < SatTb::NTERM; r_++){
                    printf("\n      ATSAT: [TB] %-15s", tb_term[r_]);
                    for(int l = 0; l < nl_; l++) if(lv_[l] < im) printf(" %+8.3f", tb_tsum[r_][lv_[l]]);
                }
            }
            printf("\n");
        }

        if(ub_on && iter_n % ub_every == 0){
            static const int lv_[] = { 0, 1, 2, 3, 5, 8, 12, 16, 20, 24, 28, 32, 35, 37, 38, 39, 40 };
            const int nl_ = (int)(sizeof(lv_) / sizeof(lv_[0]));
            const ub_tab f_now = ub_force(false), m_now = ub_force(true);
            for(int r_ = 0; r_ < SatUb::NREG; r_++){
                auto row_ = [&](const char* name, const std::vector<double>& x_, const char* fmt){
                    printf("\n      ATSAT: [UB] %-16s", name);
                    for(int l = 0; l < nl_; l++) if(lv_[l] < im) printf(fmt, x_[lv_[l]]); };
                printf("\n      ATSAT: [UB] radial wind, %s (rows %d..%d), cumulative over iterations %d..%d, m/s\n      ATSAT: [UB] %-16s",
                       ub_reg[r_], r_ == SatUb::R_ALL ? 1 : SatUb::st().j_lo, r_ == SatUb::R_ALL ? jm-2 : SatUb::st().j_hi, 1, iter_n, "level");
                for(int l = 0; l < nl_; l++) if(lv_[l] < im) printf(" %9d", lv_[l]);
                row_("u now", ub_last[r_], " %+9.5f");
                std::vector<double> d_(im, 0.0);
                for(int i = 0; i < im; i++) d_[i] = ub_last[r_][i] - ub_first[r_][i];
                row_("change", d_, " %+9.5f");
                for(int q_ = 0; q_ < UB_N; q_++){
                    row_(ub_name[q_], ub_sum[q_][r_], " %+9.5f");
                    if(q_ == UB_RK_RESET) for(int t_ = 0; t_ < SatUb::NTERM; t_++) row_(ub_term[t_], ub_tsum[t_][r_], " %+9.5f");
                }
                printf("\n      ATSAT: [UB] the radial force in the same region, m/s per iteration (x 1e3)");
                std::vector<double> x_(im, 0.0);
                for(int i = 0; i < im; i++) x_[i] = 1.0e3 * f_now[r_][i];                 row_("-dp/dr now", x_, " %+9.4f");
                for(int i = 0; i < im; i++) x_[i] = 1.0e3 * m_now[r_][i];                 row_("(v2+w2)/r now", x_, " %+9.4f");
                for(int i = 0; i < im; i++) x_[i] = 1.0e3 * (f_now[r_][i] + m_now[r_][i]); row_("sum of the two", x_, " %+9.4f");
                for(int i = 0; i < im; i++) x_[i] = 1.0e3 * ub_m_first[r_][i];            row_("(v2+w2)/r start", x_, " %+9.4f");
                for(int i = 0; i < im; i++) x_[i] = 1.0e3 * ub_f_first[r_][i];            row_("-dp/dr at start", x_, " %+9.4f");
                for(int i = 0; i < im; i++) x_[i] = 1.0e3 * ub_f_solver[r_][i];           row_("  solver", x_, " %+9.4f");
                for(int i = 0; i < im; i++) x_[i] = 1.0e3 * ub_f_filter[r_][i];           row_("  filter on p", x_, " %+9.4f");
                printf("\n");
            }
        }

        if(cwb_on && iter_n % cwb_every == 0){
            for(int g_ = 0; g_ < SatCwb::NGAS; g_++){
                double sum_ = 0.0;
                printf("\n      ATSAT: [CWB] %s column (vapour+cloud+ice), cumulative over iterations %d..%d, g/m2:  start %.9e  now %.9e  change %+.6e (%+.4f %%)\n",
                       cw_gas[g_], iter_start, iter_n, cw_first[g_], cw_last[g_], cw_last[g_] - cw_first[g_], 100.0 * (cw_last[g_] / cw_first[g_] - 1.0));
                const double p0_ = cwb_col(cw_f[g_], 1) - cw_first0[g_], in_ = (cw_last[g_] - cw_first[g_]) - p0_;
                printf("      ATSAT: [CWB] %s   of the change: plane i=0 (set by bcRadius) %+.6e, integrated cells i>=1 %+.6e (%+.4f %% of their start)\n",
                       cw_gas[g_], p0_, in_, 100.0 * in_ / (cw_first[g_] - cw_first0[g_]));
                for(int q_ = 0; q_ < TB_N; q_++){
                    printf("      ATSAT: [CWB] %s   %-46s %+.6e\n", cw_gas[g_], cw_name[q_], cw_sum[g_][q_]); sum_ += cw_sum[g_][q_];
                    if(q_ == TB_RK_RESET)
                        for(int r_ = 0; r_ < SatCwb::NTERM; r_++){
                            printf("      ATSAT: [CWB] %s   %-46s %+.6e\n", cw_gas[g_], cw_term[r_], cw_tsum[g_][r_]);
                            if(r_ != SatCwb::T_QDIV) sum_ += cw_tsum[g_][r_];
                        }
                }
                printf("      ATSAT: [CWB] %s   %-46s %+.6e\n", cw_gas[g_], "unattributed", cw_last[g_] - cw_first[g_] - sum_);
            }
        }

        panorama_cnt++;

        // ATSAT_VTK_STRIDE=<n> (2026-10-09, mirrors ATJUP's): write the VTK files every n
        // iterations instead of at every checkpoint, so the printed extrema can stay dense while
        // the files do not. Default 0 = follow `checkpoint`, as before.
        static const int vtk_stride = [](){
            const char* e = getenv("ATSAT_VTK_STRIDE"); return e ? atoi(e) : 0; }();
        if(iter_n % checkpoint == 0) printMinMax();
        if(vtk_stride > 0 ? (iter_n % vtk_stride == 0) : (iter_n % checkpoint == 0)) writeData();

        // Shapiro filter on the velocities, opt-in. ATJUP applies one at initialisation and
        // optionally n passes per iteration; ATSAT already filters p_dyn and the mass fluxes
        // this way but has never filtered u, v, w. Grid-scale checkerboard in the velocity is
        // what the pressure projection cannot see and what feeds the polar cells.
        if(shapiro_vel_inloop() > 0)
            for(int n = 0; n < shapiro_vel_inloop(); n++) dampVelocities();

        // Non-finite census, opt-in. Runs at the checkpoint cadence so it costs nothing on the
        // iterations in between.
        if(nancheck_on() && iter_n % checkpoint == 0) nan_watch(iter_n);

        // Per-iteration trace, opt-in and read-only. ATSAT_TRACE=n emits a line every n
        // iterations (1 = every one); 0, the default, emits nothing and costs one integer test.
        // Placed here, at the end of the iteration alongside the other diagnostics and after
        // restoreVar, so the line describes the state the NEXT iteration starts from — the same
        // state save_state would serialise.
        static const int trace_every = [](){
            const char* e = getenv("ATSAT_TRACE"); return e ? atoi(e) : 0; }();
        if(trace_every > 0 && iter_n % trace_every == 0) trace_line(iter_n);

        // ---- Binary restart checkpoints ----
        // One explicit dump at checkpoint_save_iter, plus a periodic one every
        // restart_save_stride iterations, which is ATJUP's arrangement and its stride of 100.
        // The periodic dump is written ONLY when the state is clean, so a diverged run can never
        // overwrite a good restart point — being able to resume from it is the file's whole
        // value. Written at the END of the iteration, after restoreVar has run, so the stored
        // n-level copies are consistent with the fields they were built from.
        if(checkpoint_save_iter >= 0 && iter_n == checkpoint_save_iter)
            save_state(iter_n);

        // ATSAT_RESTART_STRIDE=<n> (2026-10-09, mirrors ATJUP's): default 100 = the stride this
        // always had; 0 = NO periodic dump. checkpoint_save_iter = -1 switches off only the
        // explicit dump above. Each sat_restart_<iter>.bin is 756 MB, so a comparison arm that
        // will never be resumed should say 0.
        {
            static const int restart_save_stride = [](){
                const char* e = getenv("ATSAT_RESTART_STRIDE"); return e ? atoi(e) : 100; }();
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

    p_dyn.initArray(im, jm, km, pa);
    p_metric.initArray(im, jm, km, 0.0);            // ATSAT_METRIC_SPLIT
    p_dynn.initArray(im, jm, km, pa);                // dynamic pressure, previous iteration
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
    S_precip_h2o.initArray(im, jm, km, 0.0);        // precipitation source terms, moisture eqs
    S_precip_h2o_cloud.initArray(im, jm, km, 0.0);
    S_precip_h2o_ice.initArray(im, jm, km, 0.0);
    S_precip_nh3.initArray(im, jm, km, 0.0);
    S_precip_nh3_cloud.initArray(im, jm, km, 0.0);
    S_precip_nh3_ice.initArray(im, jm, km, 0.0);
    S_precip_ch4.initArray(im, jm, km, 0.0);
    S_precip_ch4_cloud.initArray(im, jm, km, 0.0);
    S_precip_ch4_ice.initArray(im, jm, km, 0.0);
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
    acc_t.initArray(im, jm, km, 0.0);   // RK4 stage accumulator
    acc_u.initArray(im, jm, km, 0.0);   // RK4 stage accumulator
    acc_v.initArray(im, jm, km, 0.0);   // RK4 stage accumulator
    acc_w.initArray(im, jm, km, 0.0);   // RK4 stage accumulator
    acc_h2o.initArray(im, jm, km, 0.0);   // RK4 stage accumulator
    acc_h2o_cloud.initArray(im, jm, km, 0.0);   // RK4 stage accumulator
    acc_h2o_ice.initArray(im, jm, km, 0.0);   // RK4 stage accumulator
    acc_ch4.initArray(im, jm, km, 0.0);   // RK4 stage accumulator
    acc_ch4_cloud.initArray(im, jm, km, 0.0);   // RK4 stage accumulator
    acc_ch4_ice.initArray(im, jm, km, 0.0);   // RK4 stage accumulator
    acc_h2s.initArray(im, jm, km, 0.0);   // RK4 stage accumulator
    acc_nh3.initArray(im, jm, km, 0.0);   // RK4 stage accumulator
    acc_nh3_cloud.initArray(im, jm, km, 0.0);   // RK4 stage accumulator
    acc_nh3_ice.initArray(im, jm, km, 0.0);   // RK4 stage accumulator
    acc_nh4sh.initArray(im, jm, km, 0.0);   // RK4 stage accumulator
    acc_tke.initArray(im, jm, km, 0.0);   // RK4 stage accumulator
    acc_dis.initArray(im, jm, km, 0.0);   // RK4 stage accumulator
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
                // p_dynn is the previous-iteration dynamic pressure. Nothing maintained it and
                // nothing even allocated it, while steadyQuery pretended to difference against
                // it. It belongs with the other n-copies.
                p_dynn.x[i][j][k] = coeff * p_dyn.x[i][j][k];
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


