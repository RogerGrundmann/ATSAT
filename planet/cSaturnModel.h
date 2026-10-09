#ifndef _CSATURNMODEL_H
#define _CSATURNMODEL_H

// Forward declarations — full definitions included at the bottom of files that
// instantiate these classes, after cSaturnModel is complete.
class ChemistrySat;
class SaturationAdjustmentSat;
template<class M> class SaturationAdjustment;
class BC_Sat;
template<class M> class BoundaryConditions;
template<class M> class FluxLimiter;
class VelocityInitializerSat;

#include <fenv.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <ctime>    
#include <cmath>
#include <map>
#include <set>
#include <limits>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <stdexcept>
#include <sys/stat.h>
#include <sys/types.h>

#include "Array.h"
#include "Array_1D.h"
#include "Array_2D.h"
#include "tinyxml2.h"
#include "PythonStream.h"
#include "Utils.h"
#include "Config.h"
#include "BoundaryConditions.h"   // BCForm, and the shared BC template


#ifdef _OPENMP
#include <omp.h>
#endif


using namespace std;

namespace{
    std::function<double(double)> default_lambda=[](double i)->double{return i;};
}

// PressureSolverSat is a typedef of the SHARED PressureSolver template, not a class of its own,
// so it cannot be forward-declared as one. The template is declared here and the alias formed
// below; a pointer to a specialization needs neither complete.
template<class M> class PressureSolver;
class cSaturnModel;
typedef PressureSolver<cSaturnModel> PressureSolverSat;

class cSaturnModel{

    friend class ChemistrySat;
    friend class SaturationAdjustmentSat;
    friend class BC_Sat;
    friend class VelocityInitializerSat;
    friend class ThermalWindDiagSat;
    template<class M> friend class PressureSolver;
    template<class M> friend class SaturationAdjustment;
    template<class M> friend class BoundaryConditions;
    template<class M> friend class FluxLimiter;
    template<class M> friend class Reporting;
    template<class M> friend class ParaViewWriter;
    template<class M> friend class ConvectiveAdjustment;
    template<class M> friend class Radiation;
    template<class M> friend class Precipitation;
    template<class M> friend class Turbulence;

public:

    const char *filename;

    cSaturnModel();
    ~cSaturnModel();

    cSaturnModel(const cSaturnModel&) = delete;
    cSaturnModel& operator=(const cSaturnModel&) = delete;

    static cSaturnModel* get_model(){
        if(!m_model){
            m_model = new cSaturnModel();
        }
        return m_model;
    }

    void LoadConfig(const char *filename);
    void Run();

    #include "SaturnParams.h.inc"

    static const double pi180, the_degree, phi_degree, dthe, dphi, dr;
    double dt = 0.0;
    static const double the0, phi0, r0;


    int iter_n, panorama_cnt;
    int i_res, j_res, k_res;

    std::vector<double> tropopause_layers; // keep the tropopause layer index
    std::vector<std::vector<int> > i_topography;
    std::vector<double> u_trans;
    std::vector<double> v_trans;
    std::vector<double> w_trans;

    /*
     * This function must be called after init_layer_heights()
     * Given a layer index i, return the height of this layer
    */
    float get_layer_height(int i){
        if(0>i || i>im-1){
            return -1;
        }
        return m_layer_heights[i];
    }
    /*
     * Thickness of layer i in METRES. m_layer_heights is built from L_atm, which is specified in
     * km, so get_layer_height() returns KILOMETRES; anything forming a per-metre physical quantity
     * (a volumetric heating rate W/m3, say) must use this accessor and not the raw difference.
     */
    /*
     * Metric radius (ATSAT_METRIC_RADIUS, in km; default 0 = off, every run bit-identical).
     *
     * rad.z runs 1..2, so every 1/r factor in the equations works on a sphere of about one length
     * unit instead of on Saturn — the curvature and divergence terms are then too large by the
     * ratio of the planetary radius to the shell thickness, here 58232/500 = 116. ATJUP had
     * exactly this and it mattered: with the radius switched on there, a growing mode that killed
     * every long run simply stopped growing (2026-07-29, ATJUP_METRIC_RADIUS).
     *
     * Only the 1/r FACTORS are shifted. rad.z itself must not be moved: it is also the stretched
     * radial coordinate that the layer heights and any exp(zeta*(rad.z-1)) stretching sit on, and
     * shifting it would divide every radial derivative by ~116 and overflow the stretching. This
     * is the same construction ATOM uses (metricRadius at its five geometry sites) and the reason
     * is recorded there too.
     *
     * Pass the radius in km, e.g. ATSAT_METRIC_RADIUS=58232.
     */
    double metricRadius(double rm){
        static const double R_km = [](){
            const char* e = getenv("ATSAT_METRIC_RADIUS"); return e ? atof(e) : 0.0; }();
        if(!(R_km > 0.0)) return rm;
        return rm + (R_km / L_atm - 1.0);
    }

    /*
     * THE one gate on which density the model uses, mirroring ATJUP's local_rho(). Default
     * FALSE, i.e. the constant reference density r_mix — which is ATJUP's default too, and the
     * reason it is the default is that the local density is a strong feedback in the schemes
     * that saturate on it. Set ATSAT_LOCAL_RHO=1 to move every gated caller at once.
     */
    /*
     * Polar metric floor. sin(theta) is held at this value instead of going to zero, so that the
     * 1/sin(theta) and 1/sin^2(theta) metric terms stay bounded near the poles.
     *
     * ATSAT HAS NEVER HAD ONE, and that is the difference from ATJUP worth knowing here rather
     * than the knob itself. ATJUP floors at 0.55 — theta = 33.4 deg, i.e. active poleward of
     * 56.6 degrees latitude, 16.5% of the sphere — and its comment records that the value was
     * RAISED from 0.4 to stop a long-run polar blow-up, so it is a deliberate stability trade
     * and not an oversight. ATSAT instead relies on its Runge-Kutta loop stopping two rows short
     * of each pole: the innermost integrated row is j=2, where sin(theta) = 0.035 and 1/sin^2 is
     * 820. Bounded, but three orders of magnitude above what ATJUP allows itself.
     *
     * Default 0.0 = no floor, i.e. exactly what ATSAT does today and every run bit-identical.
     * Set ATSAT_SINTHE_MIN=0.55 to run ATJUP's polar metric and find out whether ATSAT's polar
     * behaviour is the same problem.
     */
    // The model's own name, used by the SHARED physics headers for their log prefix and to
    // build their environment-variable names (ATSAT_CONV_ADJ_LAPSE and so on). It is the
    // only thing those files know about which planet they are running on.
    static const char* planet_tag(){ return "ATSAT"; }

    // p_dyn is stored as the NONDIMENSIONAL kinematic pressure, so displaying it in bar needs
    // r_mix*u_0^2*1e-5. Default OFF (returns 1.0) so the printed number is unchanged; ATSAT_PDYN_UNITS=1
    // makes it actually bar. Mirrors ATJUP's accessor of the same name, which is what lets the
    // pressure row read the same in all three models.
    double p_dyn_to_bar() const {
        static const bool on = [](){ const char* e = getenv("ATSAT_PDYN_UNITS"); return e && atoi(e) != 0; }();
        return on ? r_mix * u_0 * u_0 * 1.0e-5 : 1.0;
    }


    // ---- ATSAT_QHEAT_SCALE / the latent+sensible heating scale ----
    //
    // Returns exactly 1.0 when off, so the default path is bit-identical.
    //
    // Thermo_*.cpp forms the heating as  lv * velocity_av * grad(rho) / (L_atm * L_atm), where
    // velocity_av is built from the RAW NONDIMENSIONAL u,v,w and grad() is per nondimensional
    // length. The physical volumetric rate is lv * (v . grad rho) in W/m3, which needs u_0 to
    // make the velocity a speed and ONE division by the length scale IN METRES. The code divides
    // by the length scale twice, in kilometres. The ratio is u_0 * L_atm / 1e3.
    //
    // Q_Latent and Q_Sensible are OUTPUT-ONLY — nothing in RHS_*/RungeKutta_* reads them — so
    // this changes no physics, but it does change the .vtk/.vts values, which is why it is a knob
    // and not a silent repair.
    double qheat_fix() const {
        static const int on = [](){ const char* e = getenv("ATSAT_QHEAT_SCALE"); return e ? atoi(e) : 0; }();
        return on ? u_0 * L_atm / 1.0e3 : 1.0;
    }

    // Hooks for the shared ParaViewWriter<Planet> (ParaViewWriter.h). planet_name() is the
    // word in an output FILE name ("Saturn_radial_20_1.vtk", "PlotData_Saturn.xyz");
    // planet_short() is the abbreviation inside a .vtk title line
    // ("Radial_Data_Sat_Circulation"). Both models carried both spellings by hand.
    static const char* planet_name(){ return "Saturn"; }
    static const char* planet_short(){ return "Sat"; }
    // What the panorama .vts prints in its "Temperature" array. ATSAT writes KELVIN/10 here;
    // ATJUP writes degrees Celsius (t*t_ref - 273.15). Two different quantities under one
    // array name, so a ParaView state file coloured for one model mis-scales the other. This
    // hook names the divergence; settling it belongs with the rest of the units question.
    double paraview_temperature(double t_nd) const { return t_nd * t_ref / 10.0; }

    // ---- The surface of a column, for the SHARED physics headers ----
    //
    // ATSAT HAS NO OBSTACLE — BC_seamount is never called and i_topography is declared but never
    // sized, so it must not be read. Every column starts at i = 0 and no cell is solid.
    //
    // Every shared header reaches the topography through these two and never touches SeaMount or
    // i_topography directly, which is what lets one implementation serve a model with an obstacle
    // and one without. They were extracted when Turbulence was shared: all six real differences
    // between TurbulenceJup.h and TurbulenceSat.h were this one concept, spelled out inline.
    int  surface_index(int, int) const { return 0; }
    bool is_solid(int, int, int) const { return false; }

    // ---- What the SHARED BoundaryConditions.h asks of this model ----
    //
    // The field lists, the loop margin, the default extrapolation form and each knob's default.
    // See BoundaryConditions.h for why every one of these is a model fact rather than a variant
    // of the algorithm. ATSAT works the INTERIOR rows (margin 1) and defaults to the 3-point
    // cubic; every hardening knob is OFF, because on Saturn none of them is measured yet.
    static int bc_margin(){ return 1; }
    static int bc_default_form(){ return BCForm::CUBIC; }
    static int bc_default_rigid_lid(){ return 0; }
    static int bc_default_top_taper(){ return 0; }
    static int bc_default_pole_copy(){ return 0; }
    static int bc_default_radius_copy(){ return 0; }

    std::vector<Array*> bc_fields_radius();
    std::vector<Array*> bc_fields_theta_extrap();
    std::vector<Array*> bc_fields_theta_zero();
    std::vector<Array*> bc_fields_phi();
    std::vector<Array*> bc_turb_fields();
    std::vector<double> bc_turb_floors();
    bool bc_turb_active() const { return turb_active; }

    // ---- What the SHARED SaturationAdjustment.h asks of this model ----
    //
    // ATSAT's saturation adjustment ends each cell by rewriting p_stat from the adjusted
    // temperature; ATJUP does not touch p_stat there at all. That is a real disagreement about
    // where the hydrostatic pressure may respond to latent heating, it was flagged and left
    // unsettled when the mirrored routine was written, and sharing the file was not the moment to
    // settle it. Answering true here keeps ATSAT doing exactly what it did.
    static bool satadj_updates_pstat(){ return true; }
    // Defaults of ATSAT_SATADJ_NEWTON and ATSAT_SATADJ_CONSERVE (shared SaturationAdjustment.h): ON since
    // 2026-10-09. They act only when the shared routine is selected (ATSAT_SATADJ=1; the inherited
    // routine is still this model's default). Run then for 16 iterations: without them the shared
    // routine deleted ice above the melting point by the per cent of a gas's column and left a
    // fifth to a half of its cell-calls unconverged; with them both are zero.
    static bool satadj_default_newton(){ return true; }
    static bool satadj_default_conserve(){ return true; }

    // ---- What the SHARED PressureSolver.h asks of this model ----
    //
    // has_obstacle() is the same fact is_solid() states cell by cell, asked once: ATSAT contains
    // no solid body, so the solver's wall condition on p_dyn defaults off. With nothing solid the
    // walled and unwalled stencils are the same sum anyway.
    static bool has_obstacle(){ return false; }

    // Rigid radial walls on aux_u, default OFF here where ATJUP has them ON. ATSAT's i=0 is the
    // deep interior of a gas giant rather than a floor, and whether a lid belongs there at all is
    // a modelling question the port did not get to settle. ATSAT_BC_RIGID_LID=1 turns it on.
    static bool press_rigid_lid(){ return false; }

    // The model's own boundary conditions on aux_*/rhs_* before the projection takes their
    // divergence. Defined in Pressure_Sat.cpp; these are computePressure()'s conditions, kept
    // because they are ATSAT's and the shared solver reads the fields they write.
    void prepareProjectionBoundaries(bool rigid_lid);

    /*
     * ---- Saturn's radiative constants, for the SHARED Radiation.h ----
     *
     * These five are MEASURED PROPERTIES OF SATURN, which is why they live here and not in the
     * shared physics file: a mechanical copy of the Jupiter radiation would have produced
     * Jupiter's radiation budget on Saturn's grid.
     *
     *                            ATJUP        ATSAT      source
     *     solar constant         50.5         14.83      1361 / a^2, a = 5.20 vs 9.58 AU
     *     Bond albedo            0.343        0.342      measured
     *     intrinsic flux F_int   5.4          2.01       measured
     *     x_H2 / x_He            0.863/0.134  0.96/0.032 Saturn's upper atmosphere is He-poor
     *
     * The grey-opacity CALIBRATION (C_cia, the kappa values, opac_cal) is deliberately NOT here —
     * it is one shared set, tuned on Jupiter, and Radiation.h says why. Saturn's lever on it is
     * ATSAT_CIA_STRENGTH / ATSAT_OPACITY_STRENGTH, and the first Saturn measurement put the
     * thermal photosphere at 0.054 bar against Jupiter's 0.5 bar target, so that lever has a
     * question waiting on it.
     *
     * NOT ACCOUNTED FOR: the shared file's insolation is an ANNUAL MEAN. Saturn's obliquity is
     * 26.7 degrees against Jupiter's 3.1, so a real seasonal cycle is being averaged away here in
     * a way it is not on Jupiter.
     */
    static double rad_F_int()       { return 2.01;  }  // Saturn intrinsic heat flux [W/m2]
    static double rad_S_solar()     { return 14.83; }  // solar constant at 9.58 AU [W/m2]
    static double rad_albedo_bond() { return 0.342; }  // Saturn Bond albedo (so S*(1-A) is ABSORBED)
    static double rad_x_H2()        { return 0.96;  }  // H2 mole fraction
    static double rad_x_He()        { return 0.032; }  // He mole fraction (depleted vs Jupiter)

    static double sinthe_min(){
        static const double v = [](){
            const char* e = getenv("ATSAT_SINTHE_MIN");
            const double x = e ? atof(e) : 0.0;
            return (x >= 0.0 && x < 1.0) ? x : 0.0;
        }();
        return v;
    }

    // ---- Hooks for the shared Reporting<Planet> (Reporting.h) ----
    // Each is a divergence between ATSAT and ATJUP that the shared reporting code would
    // otherwise have had to choose between. Every default here is ATSAT's existing behaviour,
    // so adopting the shared header changes no output.

    // ATJUP flips cos(theta) in the southern hemisphere for the continuity residual; ATSAT never
    // has. Same knob shape as ATJUP's so the two models answer one question rather than differ
    // by a missing line. Default off = unchanged.
    static bool costhe_abs(){
        static const bool v = [](){ const char* e = getenv("ATSAT_COSTHE_ABS"); return e && atoi(e) != 0; }();
        return v;
    }

    // Column layout of the min/max report. ATJUP widened its unit column to 12 because its unit
    // strings are longer (" kg/(m3s)") and separates the max and min halves with three spaces;
    // ATSAT has always used 6 and ten spaces. Purely cosmetic, and preserved rather than unified
    // so that sharing the machinery changes no log line. Unifying is a separate decision.
    static int minmax_unit_width()      { return 6; }
    static const char *minmax_separator(){ return "          "; }

    static const char *steady_heading(){
        return " 3D iterational process for the surface boundary conditions\n printout of maximum and minimum absolute and relative errors of the computed values at their locations: level, latitude, longitude";
    }

    // The iteration line of the steady-state header, including its trailing newlines. ATSAT
    // counts with iter_n and prints one; ATJUP counts with n and prints two.
    std::string steady_iter_line() const {
        return "      iter_n = " + std::to_string(iter_n) + "\n";
    }

    static bool local_rho(){
        static const bool v = [](){
            const char* e = getenv("ATSAT_LOCAL_RHO"); return e && atoi(e) != 0; }();
        return v;
    }

    /*
     * Mixture density for the cell, in kg/m3. Reads the rho_mix array that
     * computeMixtureDensity() fills once per physics block; falls back to r_mix wherever that
     * array holds nothing usable, so no caller can divide by zero. Exactly ATJUP's rho_at().
     *
     * Callers that must ALWAYS see the local field — the radiative heating, the thermal-wind
     * residual — read rho_mix directly and bypass this gate, as ATJUP's do.
     */
    double rho_at(int i, int j, int k){
        if(!local_rho()) return r_mix;
        const double rho = rho_mix.x[i][j][k];
        return (rho > 0.0 && std::isfinite(rho)) ? rho : r_mix;
    }

    // Fills rho_mix from the ideal gas law with the model's own R_mix (Pressure_Sat.cpp).
    void computeMixtureDensity();

    // Horizontal (area-weighted) mean of the buoyancy expression at each level, recomputed once
    // per RK4 step by computeBuoyancyRefLevel() in RungeKutta_Sat_Turb.cpp. Subtracting it is
    // what makes the buoyancy an ANOMALY: zero mean at every height, so only horizontal density
    // contrasts drive vertical motion and hydrostatic balance carries the mean.
    std::vector<double> buoy_ref_level;
    void computeBuoyancyRefLevel();

    // ---- Numerical safety nets, ported from ATJUP. All off by default. ----
    // ATSAT_NANCHECK: census of non-finite cells, per field, with the index extent.
    bool nan_watch(int iter);
    // ATSAT_TRACE: one read-only line per iteration, for dating the onset of things the
    // checkpoint-cadence printMinMax is too coarse to see. See FileIO_Sat.cpp.
    void trace_line(int iter);
    // ATSAT_VEL_SHAPIRO_INLOOP / ATSAT_SHAPIRO_STRENGTH: 1-2-1 filter on u, v, w.
    void dampVelocities();

    // Zero floor for the condensable species, with accounting — see FileIO_Sat.cpp.
    void clampNegativeSpecies();
    void reportClampBudget();
    std::vector<double> clamp_added;        // cumulative mass added by the floor, per field
    std::vector<long>   clamp_cells;        // cumulative number of clipped cells, per field
    std::vector<double> clamp_added_bnd;    // ... of which on the two radial boundary planes
    std::vector<long>   clamp_cells_bnd;    // ... which bcRadius extrapolates rather than integrates

    // ---- Restart checkpoint, ported from ATJUP (see FileIO_Sat.cpp) ----
    std::vector<Array*> restart_arrays();   // the prognostic 3D fields a checkpoint serialises
    void save_state(int iter);              // dump them to output_path/sat_restart_<iter>.bin
    bool load_state(int iter);              // restore them; false (run from scratch) if absent
    bool restart_state_is_clean();          // true when every serialised field is finite

    double layer_thickness_m(int i){
        if(i < 0 || i > im-2) return 0.0;
        return (double)(m_layer_heights[i+1] - m_layer_heights[i]) * 1.0e3;
    }   



private:

    static cSaturnModel* m_model;

    // Cached friend-class instances — allocated on first use, deleted in destructor
    ChemistrySat*    m_chem     = nullptr;
    PressureSolverSat* m_pressure = nullptr;

    ChemistrySat&     getChemistry();
    PressureSolverSat& getPressureSolver();

    PythonStream ps;
    std::streambuf *backup;

    const double c43 = 4.0/3.0, c13 = 1.0/3.0, c32 = 3.0/2.0, c42 = 4.0/2.0, c12 = 1.0/2.0;
    static const int im = 41, jm = 181, km = 361;

    int i_max = 40;  // corresponds to about 125 km above 10e6 Pa pressure level, maximum hight of the tropopause at equator
    int i_beg = 30;  // corresponds to about 100 km above 10e6 Pa pressure level, maximum hight of the tropopause at poles

    double mue_mix, k_mix, cp_mix, rg_mix, r_mix, R_mix, c_mix, x_mix;

// at 230K NH3 and H2S condense via a heterogenious reaction: NH3 + H2S -> NH4SH ( Planetary Scienses, de Pater, Lissauer)


// temperatures at triple point and ice formation
    double t_0_ch4 = 90.69;  // in K == -182.456°C, triple point
    // WAS 190.56 K, which is methane's CRITICAL temperature, not an ice-cloud bound — it sat 100 K
    // ABOVE t_0_ch4 and so inverted the ordering every other species has. The consequence was exact
    // and total: PrecipitationJup/Sat gate ice autoconversion on (T < t_frz && T >= t_low), which
    // for CH4 read (T < 90.69 && T >= 190.56) — an EMPTY band, so methane ice could never convert
    // to snow at any temperature. ATSAT was carrying a 45.6 g/m3 methane ice deck with zero methane
    // snow, and its clamp budget showed ch4_ice clipping 110 % of its own mass per 12 iterations
    // because the field had no sink at all.
    //
    // 67.36 K is t_0_ch4 scaled by the H2O/NH3 proportion, on Roger's instruction: H2O sits at
    // 210.15/273.15 = 0.7694 of its triple point and NH3 at 140.0/195.5 = 0.7161, mean 0.7427, and
    // 0.7427 * 90.69 = 67.36. That makes the CH4 ice band 67.36 .. 90.69 K. It is a proportion
    // carried across from two other substances, not a measured property of methane ice.
    double t_00_ch4 = 67.36;   // in K == -205.79 degC, ch4-ice cloud formation

    double t_0_h2o = 273.15;  // in K == 0°C, triple point
    double t_00_h2o = 210.15;  // in K == -67°C (Planetary Siences)
    double t_000 = 235.15;  // in K == -20°C (precipitation module)

    double t_0_h2s = 187.66;  // in K == -53.15°C, triple point
    double t_00_h2s = 140.0;  // in K == -133.15°C, nh3-ice cloud formation (Planetary Sciences, p. 96)

    double t_0_nh3 = 195.5;  // in K == -77.65°C, triple point, gas and liquid pase
    double t_00_nh3 = 140.0;  // in K == -133.15°C, nh3-ice cloud formation (Planetary Sciences, p. 96)

    double t_0_nh4sh = 230.0;  // in K == -43.15°C, nh4sh formation onset (Planetary Sciences, p. 96)
    double t_00_nh4sh = 140.0;  // in K == -133.15°C, nh3-ice cloud formation (Planetary Sciences, p. 96)


// pressures take from  Planetary Sciences, p. 96
    double p_0_ch4 = 0.1;  // in bar
    double p_00_ch4 = 1.1;  // in bar

    double p_0_h2o = 2.0;  // in bar
    double p_00_h2o = 5.1;  // in bar

    double p_0_h2s = 1.3;  // in bar
    double p_00_h2s = 2.2;  // in bar

    double p_0_nh3 = 0.32;  // in bar
    double p_00_nh3 = 2.2;  // in bar

    double p_0_nh4sh = 1.3;  // in bar
    double p_00_nh4sh = 2.2;  // in bar


// constants for Clausius-Clapeyron law
    double coeff_h2_A = -2000.0; // invented
    double coeff_h2_B = 8.0; // invented

    double coeff_he_A = -3000.0; // invented
    double coeff_he_B = 10.0; // invented


    double coeff_ch4_A = -1033.3;  // from triple and critical point values for methane
    double coeff_ch4_B = 6.3910;
    double coeff_ch4_A_i = -1033.3; // invented
    double coeff_ch4_B_i = 6.3910;  // invented


    double coeff_h2o_A = -4961.04;  // from triple and critical point values for water
    double coeff_h2o_B = 13.0662;  // from triple and critical point values for water

    double coeff_h2o_A_i = -4961.04; // invented
    double coeff_h2o_B_i = 13.0662; // invented


    double coeff_h2s_A = -2251.66;  // from triple and critical point values for hydrogen sulfide
    double coeff_h2s_B = 10.5253;  // from triple and critical point values for hydrogen sulfide

    double coeff_h2s_A_i = -2251.66; // invented
    double coeff_h2s_B_i = 10.5253; // invented


    double coeff_nh3_A = -2836.56;  // from triple and critical point values for ammonia
    double coeff_nh3_B = 11.7271;  // from triple and critical point values for ammonia

    double coeff_nh3_A_i = -2836.56; // invented
    double coeff_nh3_B_i = 11.7271; // invented


    double coeff_nh4sh_A = -2836.56;  // invented
    double coeff_nh4sh_B = 11.7271; // invented


// mass density (concentration) of species in g/cm³ == 10e6 g/m³ == 10e3 kg/m³                            Planetary sciences p. 90 2010
    double rg_h2 = 0.408; // mass density of vapour
    double rg_he = 0.1752; // mass density of vapour
    double rg_ch4 = 0.657; // mass density of vapour
    double rg_nh3 = 0.7623; // mass density of vapour
    double rg_nh4sh = 1170.0; // mass density of vapour
    double rg_h2s = 1.5357; // mass density of vapour
    double rg_h2o = 0.005; // mass density of vapour

// molecular weights
    double m_h2 = 2.016;  // molecular weight of hydrogen in kg/Kmol (molar mass)
    double m_he = 4.02602;  // molecular weight of helium in kg/Kmol
    double m_ch4 = 16.042;  // molecular weight of methane in kg/Kmol
    double m_nh3 = 17.03052;  // molecular weight of ammonia in kg/Kmol

    // Reaction enthalpy of NH3(g) + H2S(g) -> NH4SH(s), per kg of NH4SH formed.
    //
    // From standard enthalpies of formation:
    //     dHf NH3(g)    -45.90 kJ/mol
    //     dHf H2S(g)    -20.60 kJ/mol
    //     dHf NH4SH(s) -156.90 kJ/mol
    //     dH_rxn = -156.90 - (-45.90 - 20.60) = -90.40 kJ/mol   (exothermic)
    // divided by m_nh4sh = 51.1114 kg/kmol:  90.40e3 / 0.0511114 = 1.7687e6 J/kg.
    //
    // Sign convention here is HEAT RELEASED, so the source term is + dh_nh4sh * w_nh4sh and is
    // positive where NH4SH is forming. w_nh4sh is kg/(m3 s), so the product is W/m3 — the same
    // dimension as the enthalpy-diffusion term it is added to, which is the whole point.
    double dh_nh4sh = 1.7687e6;   // reaction enthalpy of NH4SH formation in J/kg
    double m_nh4sh = 51.1114;  // molecular weight of ammonium hydrosulfide in kg/Kmol
    double m_h2s = 34.08088;  // molecular weight of hydrogen sulfide in kg/Kmol
    double m_h2o = 18.01588;  // molecular weight of water in kg/Kmol

// vapour mass densities of gases                Planetary sciences p. 90 2010
    double r_h2 = 0.864;  // density of hydrogen vapour in kg/m³
    double r_he = 0.136;  // density of helium vapour  in kg/m³
    double r_ch4 = 0.19;  // density of methane vapour in kg/m³
    double r_nh3 = 0.04;  // density of ammonia vapour in kg/m³
    double r_h2s = 0.055;  // density of hydrogen sulfide vapour in kg/m³
    double r_h2o = 0.09;  // density of water vapour in kg/m³
//    double r_h2o = 55.0;  // density of water vapour in kg/m³
    double r_nh4sh = 0.007;  // density of ammonium hydrosulfide vapour in kg/m³  assumption
    double r_nh3_add = 0.09;  // density of ammonia vapour in kg/m³
/*
// vapour mass densities of clouds and ices               Planetary sciences p. 90 2010
    double r_nh3_cloud = 0.004;  // density of ammonia cloud in kg/m³
    double r_h2o_cloud = 0.009;  // density of water cloud in kg/m³
    double r_nh3_ice = 0.0003;  // density of ammonia ice in kg/m³
    double r_h2o_ice = 0.002;  // density of water ice in kg/m³
*/
// vapour molar densities of gases
    double c_h2 = r_h2/m_h2;  // density of hydrogen vapour in kmol/m³
    double c_he = r_he/m_he;  // density of helium vapour  in kmol/m³
    double c_ch4 = r_ch4/m_ch4;  // density of methane vapour in kmol/m³
    double c_nh3 = r_nh3/m_nh3;  // density of ammonia vapour in kmol/m³
    double c_h2s = r_h2s/m_h2s;  // density of hydrogen sulfide vapour in kmol/m³
    double c_h2o = r_h2o/m_h2o;  // density of water vapour in kmol/m³
    double c_nh4sh = r_nh4sh/m_nh4sh;  // density of ammonium hydrosulfide vapour in kmol/m³  assumption

 // ratio of vapour molecular weight to mean molecular weight
    double X_h2 = 0.864;
    double X_he = 0.136;
    double X_ch4 = 3.6e-5;
//    double X_h2o = 5.0e-5;
    double X_h2o = 1.7e-3;
    double X_nh3 = 2.0e-4;
    double X_h2s = 7.7e-5;
    double X_nh4sh = 3.6e-5;

// gas constants
    double R_h2 = 4124.2; // gas constant of hydrogen in J/(kg*K)
    double R_he = 2077.1; // gas constant of helium in J/(kg*K)
    double R_ch4 = 518.28; // gas constant of methane in J/(kg*K)
    double R_nh3 = 488.21; // gas constant of ammoinia in J/(kg*K)
    double R_nh4sh = 261.0; // gas constant of ammoinium hydrosulfide in J/(kg*K)   invented for the initial distribution of nh4sh
    double R_h2s = 243.96; // gas constant of hydrogen sulfide inJ/(kg*K) 
    double R_h2o = 461.52; // gas constant of water inJ/(kg*K)

// dynamic viscosities
    double mue_h2 = 0.84e-5; // dynamic viscosity of hydrogen in Ns/m²
    double mue_he = 1.87e-5; // dynamic viscosity of helium in Ns/m²
    // ATJUP's f71ef91 applied here. Methane's gas viscosity is 1.107e-2 CENTIPOISE, and the
    // centipoise figure was entered as if it were Pa*s — a factor 1000. Water's 1.308e-3 is the
    // LIQUID value; the vapour is of order 1e-5 like every other component here.
    //
    // It matters because ThermalProperties forms mue_mix as a MASS-WEIGHTED mean, so a component
    // that is 1000x too large dominates the sum outright: mue_mix came out 1.3254e-3 Ns/m2 where
    // every honest component is of order 1e-5 — liquid-water viscosity in a hydrogen atmosphere.
    //
    // ATJUP'S COMMIT SAYS "the only consumer is the NH4SH Stokes settling velocity". That is true
    // OF ATJUP, whose rhs_t carries no thermalmassflux term. It is NOT true here: mue_mix also
    // sets nu_mix, hence D_nh3 and D_h2s, hence the diffusive fluxes j_*, hence thermalmassflux —
    // which on the ice giants dominates rhs_t by four to six orders of magnitude. The scope note
    // was right where it was written and wrong for the models it was never propagated to.
    double mue_ch4 = 1.107e-5; // dynamic viscosity of methane in Ns/m²
    double mue_nh3 = 0.92e-5; // dynamic viscosity of ammonia in Ns/m²
    double mue_nh4sh = 0.99e-5; // dynamic viscosity of ammonium sulfide in Ns/m²
    double mue_h2s = 1.3e-5; // dynamic viscosity of hydrogen sulfid in Ns/m²
    double mue_h2o = 0.9e-5;   // dynamic viscosity of water VAPOUR in Ns/m² (was 1.308e-3, the liquid)

// thermal conductivities
    double k_h2 = 0.1317; // thermal conductivity of hydrogen in W/(m*K)
    double k_he = 0.1193; // thermal conductivity of helium in W/(m*K)
    double k_ch4 = 0.0; // thermal conductivity of methane in W/(m*K)
    double k_nh3 = 0.02102; // thermal conductivity of ammonia in W/(m*K)
    double k_nh4sh = 0.02102; // thermal conductivity of ammonium hydrosulfide in W/(m*K)
    double k_h2s = 0.013; // thermal conductivity of hydrogen sulfide in W/(m*K)
    double k_h2o = 0.0187; // thermal conductivity of water in W/(m*K)

// specific heat capacities
    double cp_h2 = 14.32e3;  // specific heat capacity of hydrogen in J/(kg*K)
    double cp_he = 5.19e3;  // specific heat capacity of helium in J/(kg*K)
    double cp_ch4 = 2.232e3;  // specific heat capacity of methane in J/(kg*K)
    double cp_nh3 = 2.19e3;  // specific heat capacity of ammonia in J/(kg*K)
    double cp_nh4sh = 2.00e3;  // specific heat capacity of ammonium hydrosulfide in J/(kg*K)
    double cp_h2s = 2.24e3;  // specific heat capacity of hydrogen sulfid in J/(kg*K)
    double cp_h2o = 1.93e3;  // specific heat capacity of water in J/(kg*K)

// ratios of gas constants of dry gas to vapour or vapour molecular weight to mean atmospheric molecular weight or m/m_mix
    double ep_h2 = 0.8572;  // ratio of the gas constants of dry air to h2 non-dimensional or m/m_mix
    double ep_he = 1.7152;  // ratio of the gas constants of dry he to h2 non-dimensional
    double ep_ch4 = 7.6752;  // ratio of the gas constants of dry methane to h2 non-dimensional
    double ep_h2o = 8.1253;  // ratio of the gas constants of dry air to h2 non-dimensional
    double ep_h2s = 14.5192;  // ratio of the gas constants of dry hydrogen sulfide to h2 non-dimensional
    double ep_nh3 = 7.6752;  // ratio of the gas constants of dry ammonia to h2 non-dimensional
    double ep_nh4sh = 21.7745;  // ratio of the gas constants of dry ammonia hydrosufide to h2 non-dimensional       invented
/*
// ratios of gas constants of dry gas to vapour or vapour molecular weight to mean atmospheric molecular weight or m/m_mix
    double ep_h2 = 0.8572;  // ratio of the gas constants of dry air to h2 non-dimensional or m/m_mix
    double ep_he = 1.7152;  // ratio of the gas constants of dry he to h2 non-dimensional
    double ep_h2o = 0.8715;  // ratio of the gas constants of dry air to h2 non-dimensional
    double ep_h2s = 1.6410;  // ratio of the gas constants of dry hydrogen sulfide to h2 non-dimensional
    double ep_nh3 = 0.8236;  // ratio of the gas constants of dry ammonia to h2 non-dimensional
    double ep_nh4sh = 2.4180;  // ratio of the gas constants of dry ammonia hydrosufide to h2 non-dimensional       invented
*/
// latent heat of evaporation
    double lv_ch4 = 5.11e5;  // latent heat of ch4 evaporation in J/kg
    double lv_h2o = 2.5009e6;  // latent heat of h2o evaporation at 0°C in J/kg
    double lv_h2s = 3.5340e6;  // latent heat of h2s evaporation at -73°C in J/kg
    double lv_nh3 = 1.3720e6;  // latent heat of nh3 evaporation at -33.33 in J/kg

// latent heat of sublimation
    double ls_ch4 = 5.11e5;  // latent heat of ch4 sublimation in J/kg
    double ls_h2o = 2.8339e6;  // latent heat of h2o sublimation at 0°C in J/kg
    double ls_h2s = 7.4500e5;  // latent heat of h2s sublimation at -98°C in J/kg
    double ls_nh3 = 1.8320e6;  // latent heat of nh3 sublimation at -93.15 in J/kg

// Schmidt number
    double sc_h2 = 0.20;  // Schmidt numbert of h2o, Sc = nue/D
    double sc_he = 0.22;  // Schmidt numbert of h2o, Sc = nue/D
    double sc_ch4 = 0.99;  // Schmidt number of ch4, Sc = nue/D
    double sc_h2o = 0.61;  // Schmidt numbert of h2o, Sc = nue/D
    double sc_h2s = 0.94;  // Schmidt number of h2s, Sc = nue/D 
    double sc_nh3 = 0.61;  // Schmidt number of nh3, Sc = nue/D 
    double sc_nh4sh = 0.7;  // Schmidt number of nh4sh, Sc = nue/D 

// Prantl numbers
    double Pr = 0.72;  // Prandtl-number 

// diffusion coefficients
    double D_nh3 = 1.5e-9; // ordinary diffusion coefficient of ammonia in m*m/s 
    double D_h2s = 1.36e-9; // ordinary diffusion coefficient of hydrogen sulfid in m*m/s 
    double D_nh4sh = 1.45e-9; // ordinary diffusion coefficient of ammonium hydrosulfide in m*m/s

// thermal diffusion coefficients
    double DT_nh3 = 1.54e-9; // thermal diffusion coefficient of ammonia in kg/(s*m)                           unklar
    double DT_h2s = 1.36e-9; // thermal diffusion coefficient of hydrogen sulfid in kg/(s*m)
    double DT_nh4sh = 1.45e-9; // thermal diffusion coefficient of ammonium hydrosulfide in kg/(s*m)

// constants for saturation vapour pressure and latent heat from the original paper by Sanchez-Lavega, Perez-Hoyos and Huesco, p. 770
    double C_ch4 = 1.627;   //  in bar
    double C_h2o = 25.096;  //  in bar
    double C_nh3 = 27.863;  //  in bar
    double C_h2s = 17.064;  //  in bar
    double C_nh4sh = 75.678;  //  in bar

    double L0_ch4 = 553.1;   //  in J/g
    double L0_h2o = 3148.2;  //  in J/g
    double L0_nh3 = 2016.0;  //  in J/g
    double L0_h2s = 747.0;  //  in J/g
    double L0_nh4sh = 2915.7;  //  in J/g
/*
    double L0_h2o = 2.8339e3;  //  in J/g
    double L0_nh3 = 1.3720e3;  //  in J/g
    double L0_h2s = 3.5340e3;  //  in J/g
    double L0_nh4sh = 2.915e3;  //  in J/g
*/
// alf and bet are empirical constants for each phase
    double del_alf_ch4 = 1.002;
    double del_bet_ch4 = -4.1e-3;
    double del_alf_ch4_ice = 1.002;
    double del_bet_ch4_ice = -4.1e-3;
    double C_ch4_ice = 1.627;
    double L0_ch4_ice = 553.1;

    double del_alf_h2o = 0.0;
    double del_bet_h2o = - 8.7e-3;

//    double del_alf_nh3 = - 0.888;  // original paper by Sanchez-Lavega, Perez-Hoyos and Huesco, p. 770
    // ---- Ice-phase SVP coefficients: PLACEHOLDERS, and named ones ----
    //
    // ATSAT's parameter set has no MEASURED ice pair for H2O or NH3 — only CH4 has one, and even
    // that is a copy of its liquid pair. The shared precipitation and saturation-adjustment code
    // reads one set of names on every planet, so rather than let ATSAT's copy of the microphysics
    // quietly wire the liquid coefficients into the ice slots (which is what it did while there
    // were two copies), the substitution is written down here, once, where the parameters live.
    //
    // CONSEQUENCE, and it is not small: at Saturn temperatures the ice branch is the one that
    // matters, so every ice saturation in this model is really a liquid saturation. This is
    // limitation 2 in the README and the reason ATSAT_SATADJ is still off. Replacing these six
    // numbers with measured ones is what settles it.
    double C_h2o_ice       = C_h2o;
    double L0_h2o_ice      = L0_h2o;
    double del_alf_h2o_ice = del_alf_h2o;
    double del_bet_h2o_ice = del_bet_h2o;

    double del_alf_nh3 = - 1.2;  // approximated
    double del_bet_nh3 = 0.0;

    double C_nh3_ice       = C_nh3;        // placeholder, see the note above
    double L0_nh3_ice      = L0_nh3;
    double del_alf_nh3_ice = del_alf_nh3;
    double del_bet_nh3_ice = del_bet_nh3;

    double del_alf_h2s = 0.0;
    double del_bet_h2s = - 2.9e-3;

    double del_alf_nh4sh = - 1.760;
    double del_bet_nh4sh = 7.8e-4;
 
    std::vector<std::vector<int> > j_ellipse;
    bool has_welcome_msg_printed;

    void init_layer_heights(){
        float h = L_atm/(im-1);
        for(int i=0; i<im; i++){
            m_layer_heights.push_back(i * h);
        } 
        return;
    }

    // Everything RHSSat needs to know about where a cell is, computed once per (i,j) column by
    // RungeKuttaSat instead of four times per cell by RHSSat. Reciprocals are stored, not the
    // quantities themselves, so the right-hand sides multiply where they used to divide — that
    // is what makes this worth passing around, and it is also why the change is not
    // bit-identical. Mirrors cJupiterModel::CellGeometry.
    struct CellGeometry {
        double rm, rm2, exp_rm, exp_2_rm;
        double sinthe, sinthe2, costhe, cotanthe;
        double inv_rm, inv_rm2;
        double inv_rmsinthe, inv_rm2sinthe, inv_rm2sinthe2;
        double costhe_inv_rm2sinthe;
        double inv_2dr, inv_2dthe, inv_2dphi;
        double inv_dr2, inv_dthe2, inv_dphi2;
    };

    void SetDefaultConfig();
    void RHSSat(int i, int j, int k, const CellGeometry& geo);
    void RungeKuttaSat();
    void SaturnPlotData();
    void paraview_vtk_longal(int n, int j_longal);
    void paraview_vtk_radial(int n, int i_radial);
    void paraview_vtk_zonal(int n, int k_zonal);
    void paraview_panorama_vts(int n);
    void paraview_sphere_vts(int n);

    void searchMinMax_3D(string, string, 
        string, Array &, double coeff=1., 
        std::function< double(double) > lambda = default_lambda,
        bool print_heading=false);

    void searchMinMax_2D(string, string, 
        string, Array_2D &, double coeff=1.0);

    void print_welcome_msg();
    void print_final_msg();
    void printMinMax();
    void initMsg();
    void writeResults();
    void writeData();

//    void BC_seamount();
//    void BC_solidground();
    void resetArrays();
    void TropopauseLocation();
    void SaturnCellStructure();
    void SaturnCellStructure_new();
    void Saturn_PlotData();

    void init_tropopause_layers();
    void init_u(Array &u, int j);
    void form_diagonals(Array &a, int start, int end);
    void init_v_or_w(Array &v_or_w, int j, double coeff_trop, double coeff_sl);
    void init_v_or_w_above_tropopause(Array &v_or_w, int j, double coeff);

    void init_temperature();
    void init_PressureStatic();
    void init_PressureDynamic();
//    void init_Density();

    void init_vapour(std::string gas, double &c_tropopause,
        double &coeff_A, double &coeff_B, double &coeff_A_i, double &coeff_B_i, 
        double &t_0, double &t_00,
        double &ep, double &r, double &m,
        double &C, double &L0, double &R, 
        double &del_alf, double &del_bet, double &X,
        Array &c, Array &cloud, Array &ice, Array &cloudiness);

    void init_vapour_cloud_ice(std::string gas, double &c_tropopause,
        double &coeff_A, double &coeff_B, double &coeff_A_i, double &coeff_B_i, 
        double &t_0, double &t_00,
        double &ep, double &r, double &m,
        double &C, double &L0, double &R, 
        double &del_alf, double &del_bet,
        Array &c, Array &cloud, Array &ice, Array &cloudiness);

    void init_h2s(std::string gas, double &c_tropopause,
        double &coeff_A, double &coeff_B, double &coeff_A_i, double &coeff_B_i, 
        double &t_0, double &t_00,
        double &ep, double &r, double &m,
        double &C, double &L0, double &R, 
        double &del_alf, double &del_bet, Array &c);

    void init_nh4sh(std::string gas, double &c_tropopause,
        double &coeff_A, double &coeff_B, 
        double &ep, double &r, double &m,
        double &C, double &L0, double &R, 
        double &del_alf, double &del_bet, Array &c);

    void Saturation_Adjustment(std::string gas, 
        double &coeff_A, double &coeff_B, double &coeff_A_i, double &coeff_B_i, 
        double &t_0, double &t_00, 
        double &ep, double &lv, double &ls, double &cp, double &r,
        double &C, double &L0, double &R, 
        double &del_alf, double &del_bet, double &m,
        Array &c, Array &cloud, Array &ice);

    void OneCategoryIceScheme();

    void ChemMassRateSat();
    void DiffMassFluxSat();
    void ThermalPropertiesSat();
    void Latent_Heat();
    void Forces();

    void steadyQuery();
    void restoreVar(double coeff);


    double Clausius_Clapeyron(double &T_K, double &A, double &B);  // temperature in °K
    double Humility_critical(double &x, double Hu_cr_max, double Hu_cr_mid);



    std::vector<int> im_tropopause; // keep the tropopause layer index

    // Snapshot of the lid temperature t.x[im-1][j][k], taken from the initial condition on the
    // first BC_radius() call. Empty until then, and filled only when ATSAT_BC_T_LID_PIN is set;
    // see the discussion at its call site in BC_Sat.cpp.
    std::vector<std::vector<double> > t_top_init;

    std::vector<float> m_layer_heights;
    std::vector<double> cloud_loc; // lateral cloudwater distribution
    std::vector<double> r_max; // lateral r_max distribution
    std::vector<double> r_max_add; // lateral r_max distribution
    std::vector<double> t_add; // lateral r_max distribution

    Array_1D rad;
    Array_1D the;
    Array_1D phi;

    Array_2D Topography; // topography
    Array_2D LatentHeat;        // areas of higher latent heat
    Array_2D precip_srf_h2o;    // surface precipitation, H2O   [kg/m2/s] (PrecipitationSat)
    Array_2D precip_srf_nh3;    // surface precipitation, NH3   [kg/m2/s]
    Array_2D precip_srf_ch4;    // surface precipitation, CH4   [kg/m2/s]
    Array_2D precip_srf_nh4sh;  // surface precipitation, NH4SH [kg/m2/s]
    Array_2D precip_srf_total;  // surface precipitation, all species [kg/m2/s]
    Array_2D vel_star;          // friction velocity u_tau at the first fluid layer [m/s]

    // Closure parameters. turb_model is now a configuration parameter like ATJUP's — declared in
    // param.py, so it appears in config_atsat.xml and is generated into SaturnParams.h.inc; it is
    // NOT declared here. ATSAT_TURB_MODEL still overrides it at runtime. coord_stretching stays a
    // member because ATSAT's radial coordinate is not stretched (init_layer_heights is linear) and
    // there is nothing to choose.
    double re_turb = 1.0;                       // = vel_star_ref*z_0/nue, set by TurbulenceSat
    double abl_height = 20000.0;                // boundary-layer height [m]
    bool   coord_stretching = false;            // ATSAT does not stretch the radial coordinate

    // THE turbulence gate. ATSAT_TURB and turb_model were two independent switches and only the
    // first of them decided anything: turb_model = "none" fell through parse_model's default and
    // silently ran k-omega SST. Both are folded into this one flag, set once in run_3D_loop_atm()
    // after the configuration and the ATSAT_TURB_MODEL override are both final, and read by
    // RHSSat, RungeKuttaSat, the three BC routines and the TurbulenceSat call sites. It is a plain
    // bool and not a predicate because RHSSat consults it once per cell per RK4 stage, where a
    // string comparison has no business being.
    bool   turb_active = false;
    Array_2D Precipitation;        // areas of higher precipitation
    Array_2D precipitable_water;// areas of precipitable water in the air
    Array_2D nh3_total;            // areas of higher nh3 concentration
    Array_2D nh3_cloud_total;    // areas of higher nh3_cloud concentration
    Array_2D nh3_ice_total;        // areas of higher nh3_ice concentration
    Array_2D aux_2D_v;            // auxilliar field v
    Array_2D aux_2D_w;            // auxilliar field w
    Array_2D tropopause_height; // local height of the tropopause

    Array t;                    // temperature
    Array u;                    // u-component velocity component in r-direction
    Array v;                    // v-component velocity component in theta-direction
    Array w;                    // w-component velocity component in phi-direction

    Array ch4;                    // methane vapour
    Array ch4_cloud;                // methane cloud
    Array ch4_ice;                    // methane ice
    Array h2o;                    // water vapour
    Array h2o_cloud;                // cloud water
    Array h2o_ice;                    // cloud ice
    Array h2s;                    // water vapour
    Array nh3;                    // nh3-vapour
    Array nh3_cloud;            // nh3-cloud
    Array nh3_ice;                // nh3-ice
    Array nh4sh;                    // nh4sh-vapour

    Array tn;                    // temperature new
    Array un;                    // u-velocity component in r-direction new
    Array vn;                    // v-velocity component in theta-direction new
    Array wn;                    // w-velocity component in phi-direction new
    Array ch4n;                    // ch4 new
    Array ch4_cloudn;            // ch4_cloud new
    Array ch4_icen;                // ch4_ice new
    Array h2on;                    // water vapour new
    Array h2o_cloudn;                    // water vapour new
    Array h2o_icen;                // cloud water new
    Array h2sn;                    // water vapour new
    Array nh3n;                    // nh3 new
    Array nh3_cloudn;            // nh3_cloud new
    Array nh3_icen;                 // nh3_ice new
    Array nh4shn;                    // nh4sh new

    Array massflux_h2s;   // mass flux h2s
    Array massflux_nh3;   // mass flux nh3
    Array massflux_nh4sh;   // mass flux nh4sh

    Array difflux_h2s;   // diffusive flux h2s
    Array difflux_nh3;   // diffusive flux nh3
    Array difflux_nh4sh;   // diffusive flux nh4sh

    Array fluxlim_nh4sh;   // TVD flux-limiter correction for nh4sh advection

    Array thermalmassflux;   // thermal massflux_h2s

    Array cloudiness_h2o; // cloudiness, N in literature
    Array cloudiness_nh3; // cloudiness, N in literature

    Array p_dyn;                // dynamic pressure
    Array p_dynn;                // dynamic pressure
    Array p_stat;                // static pressure
    Array rho_mix;              // mixture density [kg/m3], filled by computeMixtureDensity()
//    Array rho;                // density

    Array rhs_t;                // auxilliar field RHS temperature
    Array rhs_u;                // auxilliar field RHS u-velocity component
    Array rhs_v;                // auxilliar field RHS v-velocity component
    Array rhs_w;                // auxilliar field RHS w-velocity component

    Array rhs_ch4;                // auxilliar field RHS ch4
    Array rhs_ch4_cloud;        // auxilliar field RHS ch4_cloud
    Array rhs_ch4_ice;            // auxilliar field RHS ch4_ice
    Array rhs_h2o;                // auxilliar field RHS water vapour
    Array rhs_h2o_cloud;            // auxilliar field RHS cloud water
    Array rhs_h2o_ice;                // auxilliar field RHS cloud ice
    Array rhs_h2s;                // auxilliar field RHS water vapour
    Array rhs_nh3;                // auxilliar field RHS nh3
    Array rhs_nh3_cloud;        // auxilliar field RHS nh3_cloud
    Array rhs_nh3_ice;            // auxilliar field RHS nh3_ice
    Array rhs_nh4sh;                // auxilliar field RHS nh4sh

    Array aux;                // auxilliar field u-velocity component
    Array aux_u;                // auxilliar field u-velocity component
    Array aux_v;                // auxilliar field v-velocity component
    Array aux_w;                // auxilliar field w-velocity component

    Array Q_Latent;                // latent heat
    Array radiation;               // net thermal radiative flux [W/m2] (RadiationSat, diagnostic)
    Array epsilon;                 // layer emissivity 1 - exp(-tau) (RadiationSat, diagnostic)
    Array Q_rad;                   // radiative heating rate [W/m3] (RadiationSat, diagnostic)
    // Precipitation fluxes [kg/m2/s] and the latent heat they release (PrecipitationSat).
    Array P_rain;                  // H2O rain
    Array P_snow;                  // H2O snow
    Array P_graupel;               // H2O graupel
    Array P_nh3_rain;              // NH3 rain
    Array P_nh3_snow;              // NH3 snow
    Array P_nh3_graupel;           // NH3 graupel
    Array P_ch4_rain;              // CH4 rain
    Array P_ch4_snow;              // CH4 snow
    Array P_ch4_graupel;           // CH4 graupel
    Array P_nh4sh;                 // NH4SH settling flux
    Array Q_precip;                // latent heat released by precipitation [W/m3]

    // Precipitation SOURCE TERMS for the moisture equations, as NONDIMENSIONAL tendencies
    // (physical rate [kg/m3/s] already multiplied by L_atm[m]/u_0, the model's time unit), so
    // RHSSat adds them straight into rhs_*. Filled by PrecipitationSat only when
    // ATSAT_PRECIP_COUPLING is on; zero otherwise. See the note at the writeback in
    // PrecipitationSat.h for why the alternative — writing the depleted field in place — does
    // not survive the Runge-Kutta.
    // NH4SH needs no entry: its Stokes settling is already a term in rhs_nh4sh.
    Array S_precip_h2o;            // vapour source (rain evaporating back)  [nondim tendency]
    Array S_precip_h2o_cloud;      // cloud sink (autoconversion, accretion, riming)
    Array S_precip_h2o_ice;        // ice sink (ice autoconversion)
    Array S_precip_nh3;
    Array S_precip_nh3_cloud;
    Array S_precip_nh3_ice;
    Array S_precip_ch4;
    Array S_precip_ch4_cloud;
    Array S_precip_ch4_ice;
    // Turbulence closure (TurbulenceSat). k* and dis* are PROGNOSTIC when ATSAT_TURB is set:
    // RHS_Sat_Turb.cpp assembles rhs_tke/rhs_dis and RungeKutta_Sat_Turb.cpp integrates them.
    // With the closure off all of these stay identically zero and every run is bit-identical.
    Array tke;                     // turbulent kinetic energy k*      [dimensionless]
    Array dis;                     // dissipation eps* or omega*       [dimensionless]
    Array tken;                    // k* at the start of the RK4 step
    Array disn;                    // dis* at the start of the RK4 step
    Array rhs_tke;                 // tendency of k*
    Array rhs_dis;                 // tendency of dis*

    /*
     * RK4 STAGE ACCUMULATORS. y_{n+1} = y_n + dt/6 (k1 + 2k2 + 2k3 + k4), and the four k's are
     * evaluated in four separate passes over the grid, so the running sum needs somewhere to live
     * that is neither the start-of-step state (the *n arrays) nor the current stage input (the
     * live fields). One array per integrated field; see RungeKuttaSat for why the stages had to be
     * separated at all.
     */
    Array acc_t;
    Array acc_u;
    Array acc_v;
    Array acc_w;
    Array acc_h2o;
    Array acc_h2o_cloud;
    Array acc_h2o_ice;
    Array acc_ch4;
    Array acc_ch4_cloud;
    Array acc_ch4_ice;
    Array acc_h2s;
    Array acc_nh3;
    Array acc_nh3_cloud;
    Array acc_nh3_ice;
    Array acc_nh4sh;
    Array acc_tke;
    Array acc_dis;
    Array nue;                     // eddy viscosity nue* (the closure's own name)
    Array nue_t;                   // eddy viscosity as the RHS reads it [dimensionless]
    Array prod;                    // shear production P_k
    Array tke_source;              // P_k - Y_k
    Array dis_source;              // P_w - Y_w + D_w
    Array wall_nue;                // wall-adjacent eddy viscosity, filled once from the geometry
    Array Q_Sensible;            // sensible heat
    Array CoriolisForce;        // Coriolis force
    Array CentrifugalForce;             // centrifugal force
    Array BuoyancyForce;        // buoyancy force, Boussinesque approximation
    Array PresGradForce;// pressure gradient force
    Array SeaMount;             // sea mount contour

    Array w_nh3;                // reaction rate nh3
    Array w_h2s;                // reaction rate h2s
    Array w_nh4sh;                // reaction rate nh4sh

    Array w_nh3_cloud;                // reaction rate nh3_cloud
    Array w_nh3_ice;                // reaction rate nh3_ice

    Array j_nh3;                // ordinary-diffusion mass flux of nh3
    Array j_h2s;                // ordinary-diffusion mass flux of h2s
    Array j_nh4sh;              // ordinary-diffusion mass flux of nh4sh

    Array jT_nh3;                // thermo-diffusion mass flux of nh3
    Array jT_h2s;                // thermo-diffusion mass flux of h2s
    Array jT_nh4sh;              // thermo-diffusion mass flux of nh4sh
};
#endif
