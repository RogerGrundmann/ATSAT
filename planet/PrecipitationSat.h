/*
 * Saturn Atmosphere Circulation Model (ATSAT)
 * Ice / precipitation microphysics — PHASE 3 (H2O, NH3, CH4, NH4SH).
 *
 * Friend class of cSaturnModel (idiom of SaturationAdjustmentSat / RadiationSat).
 * PROVENANCE, AND WHAT OF IT IS MEASURED HERE. Ported from ATJUP's PrecipitationJup, which was
 * itself built in phases from ATOM_Precipitation's ThreeCatIceScheme (COSMO three-category
 * rain/snow/graupel), generalised there to H2O and NH3 with NH4SH settling by Stokes velocity.
 * The phase history in that file is ATJUP's, not ATSAT's, and is not repeated: here the scheme
 * arrives complete and UNVALIDATED. Nothing about it has been measured on Saturn.
 *
 * The microphysical constants below are ATOM's and ATJUP's. They are terrestrial-cloud numbers
 * carried through two ports, and they are the least transferable part of this file: an
 * autoconversion threshold in kg/m3, a fall speed in m/s and an evaporation coefficient all
 * belong to a particular atmosphere's density and drop-size distribution. They are left unchanged
 * so that what they do here is visible rather than hidden behind numbers invented for the
 * occasion — the same choice made for C_cia in RadiationSat, and for the same reason.
 *
 * ATSAT HAS NO OBSTACLE (BC_seamount is never called), so the column base is always i = 0; the
 * SeaMount scan is kept because it costs nothing and is correct if one is ever switched on.
 * i_topography must NOT be used: ATSAT declares it and never sizes it.
 */

#pragma once

#include <cmath>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

class cSaturnModel;
class Array;

class PrecipitationSat {
public:
    explicit PrecipitationSat(cSaturnModel& model) : m(model) {}

    // COSMO-style rate coefficients, CALIBRATED TO JUPITER'S ENERGY BUDGET.
    //
    // Why they are far below Earth/COSMO values: the condensate fields here are diagnosed
    // equilibrium concentrations from SaturationAdjustmentSat and reach q_c ~ 3e-2 kg/m3, some
    // 30x Earth's cloud water, sustained over a ~10 km deep cloud. Standard Kessler timescales
    // (c_c_au ~ 4e-4 /s) acting on that give P = c_c_au*q_c*H ~ 4e-2 kg/m2/s, i.e. a latent heat
    // flux Lv*P ~ 1e5 W/m2 against Saturn's total emitted 13.9 W/m2 (internal 5.4 + absorbed
    // solar ~8.5) — four orders of magnitude too much energy transport.
    //
    // The hard observable is that bound: Lv*P <= ~14 W/m2 gives P <= 14/2.5e6 ~ 5.6e-6 kg/m2/s
    // for H2O (~0.5 mm/day). c_c_au is set so P = c_c_au*q_c*H lands there, and the remaining
    // coefficients are scaled by the SAME factor (5e-5) so the relative balance between
    // autoconversion, accretion, riming and snow->graupel conversion is preserved. These are
    // therefore effective column-scale conversion efficiencies, not microphysical rates: they
    // absorb the gap between the model's diagnosed condensate loading and real cloud loading.
    static constexpr double c_c_au   = 2.0e-8;    // cloud-condensate autoconversion [1/s]
    static constexpr double c_i_au   = 5.0e-8;    // cloud-ice autoconversion [1/s]
    static constexpr double q_c_crit = 5.0e-4;    // autoconversion threshold [kg/m3]
    static constexpr double c_ac     = 1.0e-7;    // rain accretion of cloud
    static constexpr double c_rim    = 1.0e-8;    // riming of cloud by snow/graupel
    static constexpr double z_csg    = 2.5e-8;    // snow -> graupel conversion
    // NOT rescaled: rain evaporation is already bounded twice over (by the rain flux present and
    // by the saturation deficit), so it self-limits rather than setting the precipitation rate.
    static constexpr double a_ev     = 2.76e-3;   // rain evaporation coefficient

    // Phase-change handoff at the freezing level (phase 2b). Melting and freezing act on the
    // frozen/liquid flux ARRIVING from the layer above and are expressed as the fraction of
    // that flux which changes phase, so they can never convert more mass — or release more
    // latent heat — than is actually falling. f = min(1, c * |T - t_frz|): melting completes
    // over ~1/c_melt_K K below the freezing level, freezing over ~1/c_frz_K K above it.
    static constexpr double c_melt_K = 0.2;       // melt fraction per K of superheat  [1/K]
    static constexpr double c_frz_K  = 0.2;       // freeze fraction per K of supercooling [1/K]

    // Reference flux for the COSMO power-law collection terms. The normalised Rain/Snow/
    // Graupel ratios used to be P/P(bottom cell) with a 1e-6 floor, which pinned them at the
    // clamp whenever any flux existed and made them depend on the previous iteration's array
    // contents. A fixed physical scale (1e-4 kg/m2/s ~ 8.6 mm/day) removes both problems.
    static constexpr double P_ref    = 1.0e-4;    // [kg/m2/s]
    static constexpr double ratio_max = 5.0;      // clamp (geometric riming runaway otherwise)
    static constexpr double v_fall   = 2.0;       // condensate fall speed, sets residence time [m/s]

    // Per-species parameter bundle for the generalised column solver.
    struct Species {
        Array *vapour, *cloud, *ice;              // condensate fields (in/out)
        Array *P_r, *P_s, *P_g;                   // output fluxes (per species)
        Array *S_v, *S_c, *S_i;                   // source terms for the moisture RHS
        double C, L0, R, del_alf, del_bet;        // liquid/gas SVP coefficients
        double C_i, L0_i, del_alf_i, del_bet_i;   // ice SVP coefficients
        double ep, Lv, Ls, t_frz, t_low;          // gas-const ratio, latent heats, temps
    };

    // The mass half of phase 2c, mirrored from ATJUP. When on, column() stops writing the
    // depleted condensate into the fields and instead reports RATES into S_v/S_c/S_i, which
    // RHSSat adds to the moisture equations. Same knob as the latent-heat half in
    // RHS_Sat_Turb.cpp, because the two are one physical statement: the heat a conversion
    // releases and the mass it moved must enter the model together.
    static bool coupling_on(){
        static const bool v = [](){
            const char* e = getenv("ATSAT_PRECIP_COUPLING"); return e ? (atof(e) != 0.0) : true; }();
        return v;
    }

    void run();

private:
    cSaturnModel& m;
    void column(const Species& s);        // one species, three-category
    void sedimentNH4SH();                 // NH4SH crystal Stokes settling
    void surfaceMap();                    // all-species 2D surface precipitation map
};

// ---------------------------------------------------------------------------
#include "cSaturnModel.h"
#include "SaturationAdjustmentSat.h"

// Every species must have t_low < t_frz, or its ice band (T < t_frz && T >= t_low) is EMPTY and
// the ice can never convert to snow. That is not hypothetical: t_00_ch4 was methane's critical
// temperature, 100 K above its triple point, and methane ice was consequently sinkless in both
// models until 2026-07-31. Checked once, at the first call, and reported rather than asserted —
// a bad pair should be visible in the log of the run it spoiled, not abort a queued job.
static void check_phase_order(const char* gas, double t_frz, double t_low){
    if(!(t_low < t_frz))
        printf("      ATSAT: WARNING - %s has t_00 = %.2f K >= t_0 = %.2f K. The ice band is empty,"
               " so %s ice cannot convert to snow. See the note at t_00_ch4.\n",
               gas, t_low, t_frz, gas);
}

inline void PrecipitationSat::run(){
    std::cout << std::endl << "      ATSAT: PrecipitationSat (H2O+NH3+CH4 3-cat + NH4SH settling)" << std::endl;
    auto begin = std::chrono::high_resolution_clock::now();

    // Zero the shared latent-heat diagnostic; each species adds into it.
    // Q_precip and the nine moisture source terms: column() writes them per cell, so a cell it
    // does not reach this iteration must not keep last iteration's rate.
    Array* zero_me[] = { &m.Q_precip,
        &m.S_precip_h2o, &m.S_precip_h2o_cloud, &m.S_precip_h2o_ice,
        &m.S_precip_nh3, &m.S_precip_nh3_cloud, &m.S_precip_nh3_ice,
        &m.S_precip_ch4, &m.S_precip_ch4_cloud, &m.S_precip_ch4_ice };
    for(Array* a : zero_me){
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 0; i < m.im; i++)
            for(int j = 0; j < m.jm; j++)
                for(int k = 0; k < m.km; k++)
                    a->x[i][j][k] = 0.0;
    }

    // --- H2O three-category ---
    Species h2o;
    h2o.vapour = &m.h2o; h2o.cloud = &m.h2o_cloud; h2o.ice = &m.h2o_ice;
    h2o.P_r = &m.P_rain; h2o.P_s = &m.P_snow; h2o.P_g = &m.P_graupel;
    h2o.C = m.C_h2o; h2o.L0 = m.L0_h2o; h2o.R = m.R_h2o; h2o.del_alf = m.del_alf_h2o; h2o.del_bet = m.del_bet_h2o;
    // ATSAT HAS NO SEPARATE ICE-PHASE SATURATION COEFFICIENTS for H2O and NH3: its parameter set
    // carries one pair per species (C, L0, del_alf, del_bet) where ATJUP has both. The liquid
    // pair is used for the ice branch here — which is not a workaround but what ATSAT ALREADY
    // DOES: Weather_Sat.cpp:122-125 computes E_Rain and E_Ice from the same arguments. The
    // microphysics therefore saturates at exactly the value the adjustment does, which is what
    // the note further down requires of it.
    //
    // The physics that costs: with one pair, deposition and sublimation see the saturation vapour
    // pressure over LIQUID where they should see the one over ice, so the vapour-pressure
    // difference that drives the Bergeron process is understated and ice grows at the expense of
    // droplets more slowly than it should. Flagged rather than patched with invented
    // coefficients; the fix is to give ATSAT the ice pair, a parameter decision and not a port.
    h2o.C_i = m.C_h2o; h2o.L0_i = m.L0_h2o; h2o.del_alf_i = m.del_alf_h2o; h2o.del_bet_i = m.del_bet_h2o;
    h2o.ep = m.ep_h2o; h2o.Lv = m.lv_h2o; h2o.Ls = m.ls_h2o; h2o.t_frz = m.t_0_h2o; h2o.t_low = m.t_00_h2o;
    h2o.S_v = &m.S_precip_h2o; h2o.S_c = &m.S_precip_h2o_cloud; h2o.S_i = &m.S_precip_h2o_ice;
    check_phase_order("H2O", h2o.t_frz, h2o.t_low);
    column(h2o);

    // --- NH3 three-category (Saturn's main visible cloud deck) ---
    Species nh3;
    nh3.vapour = &m.nh3; nh3.cloud = &m.nh3_cloud; nh3.ice = &m.nh3_ice;
    nh3.P_r = &m.P_nh3_rain; nh3.P_s = &m.P_nh3_snow; nh3.P_g = &m.P_nh3_graupel;
    nh3.C = m.C_nh3; nh3.L0 = m.L0_nh3; nh3.R = m.R_nh3; nh3.del_alf = m.del_alf_nh3; nh3.del_bet = m.del_bet_nh3;
    nh3.C_i = m.C_nh3; nh3.L0_i = m.L0_nh3; nh3.del_alf_i = m.del_alf_nh3; nh3.del_bet_i = m.del_bet_nh3;   // see the note above
    nh3.ep = m.ep_nh3; nh3.Lv = m.lv_nh3; nh3.Ls = m.ls_nh3; nh3.t_frz = m.t_0_nh3; nh3.t_low = m.t_00_nh3;
    nh3.S_v = &m.S_precip_nh3; nh3.S_c = &m.S_precip_nh3_cloud; nh3.S_i = &m.S_precip_nh3_ice;
    check_phase_order("NH3", nh3.t_frz, nh3.t_low);
    column(nh3);

    // --- CH4 three-category ---
    //
    // The third condensable, and until now the one the scheme ignored — in ATJUP too, from which
    // this file was mirrored; both were fixed together. CH4 has vapour, cloud and ice fields, its
    // own saturation-vapour-pressure pair and its own latent heats, and SaturationAdjustmentSat
    // has been called for it on every iteration: it condensed, and the condensate then had
    // nowhere to go, so the CH4 budget could only grow.
    //
    // This matters more on Saturn than on Jupiter. Saturn's upper troposphere is cold enough for
    // methane to reach saturation, which is why the planet has a methane haze at all; Jupiter's
    // is not, so ATJUP's CH4 fluxes are expected to stay at zero and ATSAT's are not.
    Species ch4;
    ch4.vapour = &m.ch4; ch4.cloud = &m.ch4_cloud; ch4.ice = &m.ch4_ice;
    ch4.P_r = &m.P_ch4_rain; ch4.P_s = &m.P_ch4_snow; ch4.P_g = &m.P_ch4_graupel;
    ch4.C = m.C_ch4; ch4.L0 = m.L0_ch4; ch4.R = m.R_ch4; ch4.del_alf = m.del_alf_ch4; ch4.del_bet = m.del_bet_ch4;
    ch4.C_i = m.C_ch4_ice; ch4.L0_i = m.L0_ch4_ice; ch4.del_alf_i = m.del_alf_ch4_ice; ch4.del_bet_i = m.del_bet_ch4_ice;
    ch4.ep = m.ep_ch4; ch4.Lv = m.lv_ch4; ch4.Ls = m.ls_ch4; ch4.t_frz = m.t_0_ch4; ch4.t_low = m.t_00_ch4;
    ch4.S_v = &m.S_precip_ch4; ch4.S_c = &m.S_precip_ch4_cloud; ch4.S_i = &m.S_precip_ch4_ice;
    check_phase_order("CH4", ch4.t_frz, ch4.t_low);
    column(ch4);

    // --- NH4SH crystal sedimentation ---
    sedimentNH4SH();

    // --- all-species surface precipitation map ---
    surfaceMap();

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for PrecipitationSat\n", elapsed.count() * 1e-9);
    std::cout << "      ATSAT: PrecipitationSat ended" << std::endl;
}

// One condensing species, three-category (rain/snow/graupel).
inline void PrecipitationSat::column(const Species& s){
    const int im = m.im, jm = m.jm, km = m.km;
    const double t_ref = m.t_ref;
    const double Lf    = s.Ls - s.Lv;                // latent heat of fusion

    // Safety net well above the energy-budget scale (~5.6e-6 kg/m2/s), so the cap guards against
    // a numerical runaway instead of silently becoming the operative limiter as it used to.
    const double P_max_flux  = 1.0e-4;
    // Full depletion: precipitation is now allowed to remove all of the condensate it converts,
    // which is what closes the loop. The water/ammonia budget is then genuinely drawn down —
    // SaturationAdjustmentSat regenerates cloud from VAPOUR, so the vapour reservoir (resupplied
    // only by advection from the deep atmosphere) sets the equilibrium precipitation rate. With
    // the old dep_frac=0.5 the condensate was replenished faster than it could be consumed and
    // the scheme had no self-limiting mechanism at all.
    const double dep_frac    = 1.0;
    // Local copy: std::min binds by reference, and C++11 has no inline variables, so odr-using
    // the constexpr member directly would need an out-of-line definition in a header.
    const double ratio_cap   = ratio_max;

    #pragma omp parallel for collapse(2) schedule(static)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){

            int i_base = 0;
            while(i_base < m.im && m.SeaMount.x[i_base][j][k] == 1.0) i_base++;
            const int i_top  = im - 1;
            if(i_base >= i_top) continue;

            s.P_r->x[i_top][j][k] = 0.0;
            s.P_s->x[i_top][j][k] = 0.0;
            s.P_g->x[i_top][j][k] = 0.0;

            for(int i = i_top - 1; i >= i_base; i--){
                const double T   = m.t.x[i][j][k] * t_ref;
                const double p_u = m.p_stat.x[i][j][k];
                double q_v = std::max(0.0, s.vapour->x[i][j][k]);
                double q_c = std::max(0.0, s.cloud->x[i][j][k]);
                double q_i = std::max(0.0, s.ice->x[i][j][k]);

                // NOTE ON UNITS: the vapour/cloud/ice fields are mass CONCENTRATIONS [kg/m3],
                // not mixing ratios — SaturationAdjustmentSat builds them as r_mix*ep*E/p (a
                // mixing ratio scaled by the reference density) and converts their latent heat
                // with /(cp_mix*r_mix). So a conversion rate c*q is already [kg/(m3 s)]: the
                // flux increment is S*dz and the heating is Lf*S, with NO extra density factor.
                // Multiplying by a density here would double-count it.
                // cSaturnModel::saturation_vapour_pressure takes non-const references (it is
                // used with the model's own parameter members elsewhere), so the arguments are
                // copied into locals here rather than changing that signature.
                double T_ = T, C_ = s.C, L0_ = s.L0, R_ = s.R, da_ = s.del_alf, db_ = s.del_bet;
                const double E_r = m.saturation_vapour_pressure(T_, C_, L0_, R_, da_, db_);
                // Same density choice as SaturationAdjustmentSat — the two must agree or
                // the microphysics and the adjustment would saturate at different values.
                const double q_sat = (p_u > E_r) ? m.rho_at(i, j, k) * s.ep * E_r / p_u : q_v;

                // METRES — get_layer_height() is in km, and dz sets both the flux increment
                // S*dz and the volumetric heating F/dz, so the unit must be physical.
                const double dz = std::max(0.0, m.layer_thickness_m(i));
                if(dz <= 0.0){
                    s.P_r->x[i][j][k] = s.P_r->x[i + 1][j][k];
                    s.P_s->x[i][j][k] = s.P_s->x[i + 1][j][k];
                    s.P_g->x[i][j][k] = s.P_g->x[i + 1][j][k];
                    continue;
                }
                // Fluxes arriving from the layer above (already computed in this top-down sweep,
                // so the collection terms are driven by causally correct values rather than by
                // the previous iteration's contents of the same cell).
                const double F_r_in = s.P_r->x[i + 1][j][k];
                const double F_s_in = s.P_s->x[i + 1][j][k];
                const double F_g_in = s.P_g->x[i + 1][j][k];

                const double Rain    = std::min(ratio_cap, F_r_in / P_ref);
                const double Snow    = std::min(ratio_cap, F_s_in / P_ref);
                const double Graupel = std::min(ratio_cap, F_g_in / P_ref);
                const double Rain_79 = (Rain > 0.0) ? std::pow(Rain, 7.0/9.0) : 0.0;
                const double Rain_49 = (Rain > 0.0) ? std::pow(Rain, 4.0/9.0) : 0.0;

                const bool warm = (T >= s.t_frz);

                // --- condensate -> precipitation conversions ([kg/(m3 s)], see units note) ---
                const double S_c_au = (warm && q_c > q_c_crit) ? c_c_au * (q_c - q_c_crit) : 0.0;
                const double S_i_au = (!warm && T >= s.t_low)  ? c_i_au * q_i             : 0.0;
                const double S_ac   = (warm) ? c_ac * q_c * Rain_79 : 0.0;
                const double S_s_rim= (!warm) ? c_rim * q_c * Snow    : 0.0;
                const double S_g_rim= (!warm) ? c_rim * q_c * Graupel : 0.0;
                const double S_csg  = (!warm && q_c > q_c_crit) ? z_csg * q_c * Snow : 0.0;

                // --- freezing-level handoff, limited by the arriving flux ---
                // Melting turns falling snow/graupel into rain; freezing turns falling rain into
                // graupel. Both are fractions of the incoming flux, so the transfer is exact and
                // no phase change can occur where nothing is falling. This replaces the old
                // c_melt*(T-t_frz)*Snow rate, which was unbounded (it melted mass out of nothing
                // once the normalised Snow ratio saturated) and left the frozen flux to simply
                // vanish when the column crossed t_frz.
                const double f_melt = warm ? std::min(1.0, c_melt_K * (T - s.t_frz)) : 0.0;
                const double f_frz  = warm ? 0.0 : std::min(1.0, c_frz_K * (s.t_frz - T));
                const double F_s_melt = f_melt * F_s_in;     // snow    -> rain    [kg/m2/s]
                const double F_g_melt = f_melt * F_g_in;     // graupel -> rain
                const double F_r_frz  = f_frz  * F_r_in;     // rain    -> graupel

                // --- rain evaporation, limited by the rain actually present ---
                // q_sat = r_mix*ep*E_r/p_u grows without bound as p_u falls, so the raw COSMO
                // rate could evaporate far more than is there — creating vapour and a huge latent
                // cooling out of nothing (the -46000 W/m3 polar spike). Cap it at the rain flux
                // available in this cell and derive the vapour source from the capped flux.
                const double F_r_avail = std::max(0.0, F_r_in - F_r_frz + F_s_melt + F_g_melt
                                                 + (S_c_au + S_ac) * dz);
                // A second bound keeps the vapour it produces from oversaturating the cell:
                // over the residence time tau = dz/v_fall the added vapour S_ev*tau must not
                // exceed the saturation deficit, i.e. F_ev <= (q_sat-q_v)*v_fall.
                const double S_ev_raw  = (warm && Rain > 0.0 && q_v < q_sat)
                                       ? a_ev * (q_sat - q_v) * Rain_49 : 0.0;
                const double F_ev_sat  = std::max(0.0, q_sat - q_v) * v_fall;
                const double F_ev      = std::min(std::min(S_ev_raw * dz, F_r_avail), F_ev_sat);
                const double S_ev      = F_ev / dz;          // effective rate [kg/(m3 s)]

                // --- latent heat [W/m3]: riming/freezing release +Lf, melting -Lf, evaporation
                // -Lv. Flux terms become volumetric by dividing by the layer thickness. ---
                m.Q_precip.x[i][j][k] += Lf * (S_s_rim + S_g_rim + S_csg)
                                       + Lf * (F_r_frz - F_s_melt - F_g_melt) / dz
                                       - s.Lv * F_ev / dz;

                // --- bounded depletion of the condensate that fed the conversions ---
                // TWO WAYS TO APPLY IT, and only one of them survives. Mirrored from ATJUP.
                //
                // Writing the depleted value straight into the field (the `else` below) is what
                // this scheme has always done, and the Runge-Kutta throws it away: RHSSat runs
                // before it in the same iteration, but the integration is X = Xn + dt*(...) with
                // Xn the copy restoreVar() made at the END OF THE PREVIOUS iteration, which
                // predates this call. So the depleted value is overwritten a few lines later and
                // only the stage-1 right-hand side ever sees it. Measured on ATSAT: max
                // h2o_cloud over 50 iterations is 113.254038 with precipitation off and
                // 113.253984 with it on. The rain fell and the cloud never noticed.
                //
                // With the coupling on the depletion is reported as a RATE instead, and RHSSat
                // adds it to the moisture equations where the Runge-Kutta integrates it like any
                // other tendency. tau = dz/v_fall is the residence time of the falling
                // hydrometeors, NOT the model timestep, so d_c/tau is the conversion rate
                // itself, capped so it cannot remove more than the condensate present within one
                // residence time. Nondimensionalisation is L_atm[m]/u_0, the model's time unit,
                // the same factor the latent-heat half uses.
                const double tau = dz / v_fall;              // residence time of falling hydrometeors
                const double d_c = std::min(dep_frac * q_c, (S_c_au + S_ac + S_s_rim + S_g_rim + S_csg) * tau);
                const double d_i = std::min(dep_frac * q_i,  S_i_au * tau);
                if(coupling_on()){
                    const double nd = (m.L_atm * 1.0e3) / m.u_0;   // physical rate -> nondimensional
                    s.S_c->x[i][j][k] = - std::max(0.0, d_c) / tau * nd;
                    s.S_i->x[i][j][k] = - std::max(0.0, d_i) / tau * nd;
                    s.S_v->x[i][j][k] = + S_ev * nd;
                }else{
                    s.cloud->x[i][j][k]  = q_c - std::max(0.0, d_c);
                    s.ice->x[i][j][k]    = q_i - std::max(0.0, d_i);
                    s.vapour->x[i][j][k] = q_v + S_ev * tau;
                }

                // --- flux integration: incoming + phase handoff + local sources ---
                const double F_r = F_r_in - F_r_frz + F_s_melt + F_g_melt - F_ev
                                 + (S_c_au + S_ac) * dz;
                const double F_s = F_s_in - F_s_melt + (S_i_au + S_s_rim - S_csg) * dz;
                const double F_g = F_g_in - F_g_melt + F_r_frz + (S_g_rim + S_csg) * dz;

                s.P_r->x[i][j][k] = std::min(P_max_flux, std::max(0.0, F_r));
                s.P_s->x[i][j][k] = std::min(P_max_flux, std::max(0.0, F_s));
                s.P_g->x[i][j][k] = std::min(P_max_flux, std::max(0.0, F_g));
            }

            for(int i = 0; i < i_base; i++){
                s.P_r->x[i][j][k] = s.P_r->x[i_base][j][k];
                s.P_s->x[i][j][k] = s.P_s->x[i_base][j][k];
                s.P_g->x[i][j][k] = s.P_g->x[i_base][j][k];
            }
        }
    }
}

// All-species surface precipitation map [kg/m2/s], evaluated at the BASE of each column.
//
// The base is i_topography[j][k], not i=0, so any solid block is respected — column() does
// fill the solid cells below the base with the base value, but reading i_base states the intent
// and stays correct if that fill ever changes.
//
// Summing H2O, NH3 and NH4SH mass fluxes into one total mixes different condensates, which is
// why the per-species components are kept alongside it: the total answers "how much condensate
// mass reaches the bottom", and the components say which species is responsible.
//
// Interpreting the result: the solid boundary here is the SeaMount (BC_Sat.cpp BC_seamount) — a
// deliberately installed obstacle that makes a Great-Red-Spot-like storm possible in fluid-
// dynamical terms, i.e. a controlled experiment rather than an attempt at real Jovian geography.
// So "surface" means the bottom of the fluid domain, and which species lands there depends on
// where that boundary cuts the vertical structure. Over most of the globe i_base = 0 and the
// base is hot (~325 K), so only H2O arrives; over the SeaMount the base is lifted into the
// 200-230 K NH4SH stability band, and there the settling crystals dominate the total (they
// cannot reach any deeper because NH4SH decomposes back to NH3 + H2S above ~230 K). The NH4SH
// signal is therefore localised to the obstacle by construction, and is a property of the
// experiment's geometry — read the components, not just the total.
inline void PrecipitationSat::surfaceMap(){
    const int jm = m.jm, km = m.km;

    #pragma omp parallel for collapse(2) schedule(static)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_base = 0;
            while(i_base < m.im && m.SeaMount.x[i_base][j][k] == 1.0) i_base++;

            const double h2o_srf = m.P_rain.x[i_base][j][k]
                                 + m.P_snow.x[i_base][j][k]
                                 + m.P_graupel.x[i_base][j][k];
            const double nh3_srf = m.P_nh3_rain.x[i_base][j][k]
                                 + m.P_nh3_snow.x[i_base][j][k]
                                 + m.P_nh3_graupel.x[i_base][j][k];
            const double ch4_srf = m.P_ch4_rain.x[i_base][j][k]
                                 + m.P_ch4_snow.x[i_base][j][k]
                                 + m.P_ch4_graupel.x[i_base][j][k];
            const double nh4sh_srf = m.P_nh4sh.x[i_base][j][k];

            m.precip_srf_h2o.y[j][k]   = h2o_srf;
            m.precip_srf_nh3.y[j][k]   = nh3_srf;
            m.precip_srf_ch4.y[j][k]   = ch4_srf;
            m.precip_srf_nh4sh.y[j][k] = nh4sh_srf;
            m.precip_srf_total.y[j][k] = h2o_srf + nh3_srf + ch4_srf + nh4sh_srf;
        }
    }
}

// NH4SH crystals (solid, from the NH3+H2S reaction) settle by Stokes velocity.
// Local downward mass flux P_nh4sh = q_nh4sh * v_settle (q_nh4sh is already a mass
// concentration [kg/m3], see the units note in column()), computed only in the
// crystal-stability band [t_00_nh4sh, t_0_nh4sh].
inline void PrecipitationSat::sedimentNH4SH(){
    const int im = m.im, jm = m.jm, km = m.km;
    const double t_ref = m.t_ref;
    const double rho_ref = m.r_mix;   // reference density, only for the Stokes buoyancy term
    const double t_hi  = m.t_0_nh4sh;    // formation onset (upper T)
    const double t_lo  = m.t_00_nh4sh;   // formation end   (lower T)

    // Stokes settling velocity of the crystals: v = (2/9) r_p^2 (rho_c - rho_air) g / mu.
    constexpr double r_p   = 1.0e-5;     // NH4SH crystal radius [m] (r_p_nh4sh)
    constexpr double rho_c = 1200.0;     // NH4SH solid density [kg/m3]
    constexpr double mu    = 9.0e-6;     // H2/He dynamic viscosity [Pa s]
    const double v_settle  = (2.0/9.0) * r_p * r_p * (rho_c - rho_ref) * m.g / mu;   // ~0.07 m/s

    #pragma omp parallel for collapse(2) schedule(static)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            int i_base = 0;
            while(i_base < m.im && m.SeaMount.x[i_base][j][k] == 1.0) i_base++;
            // Fluid cells first, then fill the solid ones below. Doing it in one pass filled
            // i < i_base by reading P_nh4sh[i_base] before this pass had written it, i.e. from
            // the previous iteration.
            for(int i = i_base; i < im; i++){
                const double T = m.t.x[i][j][k] * t_ref;
                // q is already a crystal mass concentration [kg/m3], so the downward mass flux
                // is simply q*v_settle — no extra density factor (see units note in column()).
                const double q = std::max(0.0, m.nh4sh.x[i][j][k]);
                m.P_nh4sh.x[i][j][k] = (T >= t_lo && T <= t_hi) ? q * v_settle : 0.0;
            }
            const double base_flux = (i_base < im) ? m.P_nh4sh.x[i_base][j][k] : 0.0;
            for(int i = 0; i < i_base; i++) m.P_nh4sh.x[i][j][k] = base_flux;
        }
    }
}
