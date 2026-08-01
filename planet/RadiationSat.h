/*
 * Saturn Atmosphere Circulation Model (ATSAT)
 * Multi-layer grey-body radiation — diagnostic scaffold.
 *
 * PROVENANCE. Ported from RadiationJup (ATJUP), which in turn took only the ARCHITECTURE from
 * ATOM's MultiLayerRadiation: per-column grey layers, Stefan-Boltzmann emission sigma*T^4, layer
 * emissivity eps = 1 - exp(-tau), a two-stream up/down net-flux solve, and the flux divergence as
 * a heating rate. The Earth physics of the original (surface albedo feedback, ocean/land split,
 * Bignami and Atwater-Ball laws) applies to neither giant planet and is not here.
 *
 * WHAT CHANGED FROM THE JUPITER VERSION, AND WHY IT HAD TO. Every planetary number in the
 * radiation is Jupiter-specific, so a mechanical copy would have produced Jupiter's radiation
 * budget on Saturn's grid:
 *
 *                              ATJUP        ATSAT      source
 *     solar constant           50.5         14.8       1361 / a^2, a = 5.20 vs 9.58 AU
 *     Bond albedo              0.343        0.342      measured
 *     intrinsic flux F_int     5.4          2.01       measured
 *     x_H2 / x_He              0.863/0.134  0.96/0.032 Saturn's upper atmosphere is He-poor
 *     g                        25.9         from config (10.0)
 *
 * WHAT IS NOT CALIBRATED, AND MUST NOT BE READ AS IF IT WERE. C_cia was TUNED in ATJUP so that
 * Jupiter's thermal photosphere (tau = 1 from the top) lands near 0.5 bar. That tuning does not
 * transfer: the CIA layer optical depth goes as comp * P^2 / (T * g), and Saturn has 2.6x weaker
 * gravity and a larger H2 fraction, both of which deepen tau for the same coefficient. The value
 * below is therefore ATJUP's, carried over UNCHANGED and explicitly uncalibrated, so that the
 * first run measures where Saturn's photosphere actually falls instead of hiding the question
 * behind a number invented here. ATSAT_CIA_STRENGTH is the lever. The same applies to opac_cal
 * and the kappa values for the gas bands and clouds.
 *
 * SCOPE. The result goes to DIAGNOSTIC arrays only (radiation, epsilon, Q_rad); it does not touch
 * t or any rhs. Wiring the heating into the temperature equation is a separate step, and in ATJUP
 * it is still behind its own second knob. In radiative equilibrium the interface net flux is
 * uniform (= F_int) and Q_rad -> 0; a column out of equilibrium shows the tendency toward it, and
 * that is the intended self-test.
 *
 * THE DENSITY. Both models divide the species fields by the stored mixture density rho_mix,
 * filled by computeMixtureDensity() before the physics block. When this file was first written
 * ATSAT had no such array and formed its own density here; it has one now, so the departure is
 * gone and the two models again read the same quantity from the same place.
 */

#pragma once

#include <cmath>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif

class cSaturnModel;

class RadiationSat {
public:

    explicit RadiationSat(cSaturnModel& model) : m(model) {}

    // --- Physical constants, Saturn ---
    static constexpr double sigma  = 5.670374419e-8;  // Stefan-Boltzmann [W/m2/K4]
    static constexpr double F_int  = 2.01;            // Saturn intrinsic heat flux [W/m2]

    // --- Absorbed shortwave (solar) ---
    // Diurnally averaged insolation for a fast rotator: S*(1-A)*cos(lat)/pi, whose area-weighted
    // mean is exactly S*(1-A)/4 = 2.44 W/m2. Saturn's obliquity is 26.7 deg, far from Jupiter's
    // 3.1, so a seasonal cycle is real here and this annual mean is an approximation the Jupiter
    // version did not have to make. Equator gets 3.1 W/m2, poles 0.
    static constexpr double S_sat       = 14.83;   // solar constant at 9.58 AU [W/m2]
    static constexpr double albedo_bond = 0.342;   // Saturn Bond albedo
    static constexpr double k_sw        = 1.0;     // shortwave optical depth per bar, from the top

    // --- H2/He collision-induced absorption (CIA), grey parametrization ---
    // CIA needs two collision partners, each with number density n ~ P/T, so the volume
    // absorption goes as (P/T)^2; over a hydrostatic mass path dm = dP/g this gives a LAYER
    // optical depth dtau = C_cia * comp * (P/T) * dP / g. comp weights the H2-H2 and H2-He
    // collision-pair probabilities. C_cia IS NOT CALIBRATED FOR SATURN — see the header note.
    static constexpr double P0_phys  = 1.0e5;     // physical pressure [Pa] at the p_stat = 1 level
    static constexpr double x_H2     = 0.96;      // H2 mole fraction (Saturn)
    static constexpr double x_He     = 0.032;     // He mole fraction (Saturn, depleted vs Jupiter)
    static constexpr double he_ratio = 0.6;       // H2-He / H2-H2 grey binary-coefficient ratio
    static constexpr double C_cia    = 3.5e-6;    // carried from ATJUP, uncalibrated here

    // --- CH4/NH3 gas bands + cloud/ice continuum, additive to the CIA tau ---
    // dtau = kappa * q * dP/g with q the mass mixing ratio [kg/kg]. All of these are ATJUP's
    // values and share its calibration status: what was tuned there is the grey OPTICAL DEPTH
    // kappa*q, not kappa alone, so they are only as good as the species fields they multiply.
    static constexpr double kappa_ch4     = 0.003;
    static constexpr double kappa_nh3     = 1.5;
    static constexpr double kappa_cloud   = 25.0;
    static constexpr double kappa_ice     = 12.0;
    static constexpr double tau_cloud_cap = 0.4;
    static constexpr double opac_cal      = 0.25;

    void run();

private:
    cSaturnModel& m;
};

// ---------------------------------------------------------------------------
// Implementation (header-only, as in the ATJUP and ATOM originals).
// ---------------------------------------------------------------------------
#include "cSaturnModel.h"

inline void RadiationSat::run(){
    std::cout << std::endl << "      ATSAT: RadiationSat (grey CIA + CH4/NH3 + clouds)" << std::endl;
    auto begin = std::chrono::high_resolution_clock::now();

    const int im = m.im, jm = m.jm, km = m.km;
    const double t_ref = m.t_ref;
    const double g     = m.g;

    const double comp = x_H2 * x_H2 + he_ratio * x_H2 * x_He;
    static const double cia_mult  = [](){ const char* e = getenv("ATSAT_CIA_STRENGTH");     return e ? atof(e) : 1.0; }();
    static const double opac_mult = [](){ const char* e = getenv("ATSAT_OPACITY_STRENGTH"); return e ? atof(e) : 1.0; }();
    const double cia_coeff = C_cia * cia_mult * comp / g;

    static const int    solar_on   = [](){ const char* e = getenv("ATSAT_SOLAR");          return e ? atoi(e) : 1;   }();
    static const double solar_mult = [](){ const char* e = getenv("ATSAT_SOLAR_STRENGTH"); return e ? atof(e) : 1.0; }();
    static const double sw_tau_bar = [](){ const char* e = getenv("ATSAT_SW_TAU_PER_BAR"); return e ? atof(e) : k_sw; }();
    const double pi_ = 3.14159265358979323846;

    // Column diagnostics, so the first run can answer the calibration question it is here to ask.
    double olr_sum = 0.0, wsum = 0.0, tau1_p_sum = 0.0;
    long long n_tau1 = 0;

    #pragma omp parallel for collapse(2) schedule(static) \
        reduction(+:olr_sum, wsum, tau1_p_sum, n_tau1)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){

            // First fluid cell of the column.
            //
            // ATSAT HAS NO OBSTACLE. cSaturnModel::BC_seamount() exists (BC_Sat.cpp:856) and would
            // build one, but nothing ever calls it, so SeaMount keeps the 0.0 that initArray gave
            // it and i_base is always 0. The scan is therefore inert today; it is kept because it
            // costs nothing and is correct the moment a topography is switched on.
            //
            // What it must NOT do is read i_topography, the way ATJUP does. ATSAT declares that
            // vector (cSaturnModel.h:94) but never sizes it, so indexing it segfaults — which is
            // how this was found. Two other places read the same empty vector: the getter at
            // cSaturnModel.h:121 and Weather_Sat.cpp:635. Both are latent, and neither is repaired
            // from here.
            int i_base = 0;
            while(i_base < im && m.SeaMount.x[i_base][j][k] == 1.0) i_base++;
            const int i_top = im - 1;
            if(i_base >= i_top) continue;

            std::vector<double> B(im, 0.0);
            std::vector<double> eps(im, 0.0);
            std::vector<double> tau_l(im, 0.0);
            std::vector<double> Fd(im + 1, 0.0);
            std::vector<double> Fu(im + 1, 0.0);

            for(int i = i_base; i <= i_top; i++){
                const double T = m.t.x[i][j][k] * t_ref;
                B[i] = sigma * T * T * T * T;

                const double P_lo = m.p_stat.x[i][j][k] * P0_phys;
                const double P_hi = (i < i_top) ? m.p_stat.x[i + 1][j][k] * P0_phys : 0.0;
                double dP = P_lo - P_hi;
                if(dP < 0.0) dP = 0.0;
                const double P_mean = 0.5 * (P_lo + P_hi);
                const double dm = dP / g;

                const double tau_cia = (T > 0.0) ? cia_coeff * (P_mean / T) * dP : 0.0;

                // The density that turns the species DENSITIES into the mass mixing ratios the
                // kappa values want. It must be the CELL-CENTRED density, because the species
                // fields it divides are cell-centred: rho_mix, formed once per physics block by
                // computeMixtureDensity() from p_stat[i] and the same ideal gas law. Read
                // directly, not through rho_at(), so this is always the local field — as in
                // ATJUP, and for the same reason.
                //
                // This line used to form its own density from P_mean, the mean of the pressures
                // at the two INTERFACES of the layer. P_mean < p_stat[i] wherever pressure falls
                // upward, i.e. everywhere, so that density was systematically too small and every
                // mixing ratio derived from it too large. That was a leftover from before ATSAT
                // had a rho_mix array at all, not a deliberate departure.
                const double rho_c = m.rho_mix.x[i][j][k];
                const double inv_rho = (rho_c > 0.0 && std::isfinite(rho_c)) ? 1.0 / rho_c : 0.0;

                const double q_ch4 = std::max(0.0, m.ch4.x[i][j][k]) * inv_rho;
                const double q_nh3 = std::max(0.0, m.nh3.x[i][j][k]) * inv_rho;
                const double tau_gas = opac_mult * opac_cal * (kappa_ch4 * q_ch4 + kappa_nh3 * q_nh3) * dm;

                const double q_liq = (std::max(0.0, m.h2o_cloud.x[i][j][k])
                                    + std::max(0.0, m.nh3_cloud.x[i][j][k])
                                    + std::max(0.0, m.ch4_cloud.x[i][j][k])) * inv_rho;
                const double q_ice = (std::max(0.0, m.h2o_ice.x[i][j][k])
                                    + std::max(0.0, m.nh3_ice.x[i][j][k])
                                    + std::max(0.0, m.ch4_ice.x[i][j][k])) * inv_rho;
                double tau_cloud = opac_mult * opac_cal * (kappa_cloud * q_liq + kappa_ice * q_ice) * dm;
                if(tau_cloud > tau_cloud_cap) tau_cloud = tau_cloud_cap;

                tau_l[i] = tau_cia + tau_gas + tau_cloud;
                eps[i]   = 1.0 - std::exp(-tau_l[i]);
            }

            // Where does tau = 1 from the top fall? That is the photosphere, and the number this
            // scaffold exists to measure before anything is calibrated.
            {
                double cum = 0.0;
                for(int i = i_top; i >= i_base; i--){
                    cum += tau_l[i];
                    if(cum >= 1.0){
                        tau1_p_sum += m.p_stat.x[i][j][k];
                        n_tau1++;
                        break;
                    }
                }
            }

            std::vector<double> Fsw(im + 2, 0.0);
            double F_sun_toa = 0.0;
            if(solar_on != 0){
                const double colat   = pi_ * (double)j / (double)(jm - 1);
                const double cos_lat = std::max(0.0, std::sin(colat));
                F_sun_toa = solar_mult * S_sat * (1.0 - albedo_bond) * cos_lat / pi_;
                for(int i = i_base; i <= i_top + 1; i++){
                    const double p_bar = (i <= i_top) ? m.p_stat.x[i][j][k] : 0.0;
                    Fsw[i] = F_sun_toa * std::exp(-sw_tau_bar * std::max(0.0, p_bar));
                }
            }

            Fd[i_top + 1] = 0.0;
            for(int i = i_top; i >= i_base; i--)
                Fd[i] = Fd[i + 1] * (1.0 - eps[i]) + eps[i] * B[i];

            Fu[i_base] = Fd[i_base] + F_int + Fsw[i_base];

            for(int i = i_base; i <= i_top; i++)
                Fu[i + 1] = Fu[i] * (1.0 - eps[i]) + eps[i] * B[i];

            // Outgoing longwave at the top, area-weighted, for the balance check.
            {
                const double colat = pi_ * (double)j / (double)(jm - 1);
                const double wgt   = std::sin(colat);
                olr_sum += wgt * (Fu[i_top + 1] - Fd[i_top + 1]);
                wsum    += wgt;
            }

            for(int i = i_base; i <= i_top; i++){
                const double net_bot = Fu[i]     - Fd[i];
                const double net_top = Fu[i + 1] - Fd[i + 1];

                double dz = (i < i_top) ? m.layer_thickness_m(i) : m.layer_thickness_m(i - 1);
                if(dz <= 0.0) dz = 1.0;

                const double sw_absorbed = Fsw[i + 1] - Fsw[i];
                const double Q_sw        = sw_absorbed / dz;

                m.epsilon.x[i][j][k]   = eps[i];
                m.radiation.x[i][j][k] = 0.5 * (net_bot + net_top);
                m.Q_rad.x[i][j][k]     = (net_bot - net_top) / dz + Q_sw;
            }

            for(int i = 0; i < i_base; i++){
                m.epsilon.x[i][j][k]   = m.epsilon.x[i_base][j][k];
                m.radiation.x[i][j][k] = m.radiation.x[i_base][j][k];
                m.Q_rad.x[i][j][k]     = 0.0;
            }
        }
    }

    const double olr  = (wsum > 0.0) ? olr_sum / wsum : 0.0;
    const double tau1 = (n_tau1 > 0) ? tau1_p_sum / double(n_tau1) : 0.0;
    const double in_  = F_int + solar_mult * S_sat * (1.0 - albedo_bond) / 4.0;
    printf("      ATSAT: radiation — mean OLR %.3f W/m2 against %.3f in (F_int %.2f + solar %.2f);"
           " thermal photosphere (tau=1) at %.4f bar in %lld of %d columns\n",
           olr, in_, F_int, solar_mult * S_sat * (1.0 - albedo_bond) / 4.0,
           tau1, n_tau1, jm * km);

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for RadiationSat\n", elapsed.count() * 1e-9);
    std::cout << "      ATSAT: RadiationSat ended" << std::endl;
}
