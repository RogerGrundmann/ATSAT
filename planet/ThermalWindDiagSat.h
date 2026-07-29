/*
 * Saturn Atmosphere Circulation Model (ATSAT)
 * Thermal-wind residual — MEASUREMENT ONLY, nothing is modified.
 *
 * WHY ONLY THE MEASUREMENT. ATJUP has ThermalWindJup, which both measures this residual and
 * corrects the wind. That class cannot be ported mechanically: its coefficient 1e5/u_0 is the
 * ratio of two ATJUP-specific nondimensionalisation factors (buoyancy 1e5*L/u_0^2 against Coriolis
 * L/u_0), and ATSAT does not share that pair — its buoyancy term (RHS_Sat.cpp:294) is
 * `buoyancy * g * (p_stat + p_dyn)`, a pressure times g rather than a density anomaly, so the
 * density whose horizontal gradient the relation is about is never formed. Correcting the wind
 * here would mean rewriting ATSAT's buoyancy first, which is a modelling change and not a port.
 *
 * Measuring costs none of that. Everything below is in PHYSICAL units, formed locally, so no
 * nondimensionalisation factor from either model is involved.
 *
 * THE RELATION. Cross-differentiating the meridional geostrophic balance against the hydrostatic
 * one gives, for the zonal wind w and colatitude theta,
 *
 *     dw/dz = (g / (f * rho)) * (1/R) * drho/dtheta,      f = 2 * Omega * cos(theta)
 *
 * with R the planetary radius and rho = P/(R_mix*T) from the model's own mixture constant. A
 * planet whose jets are in thermal-wind balance satisfies it; one whose jets are barotropic while
 * its temperature field is baroclinic does not, and the size of the gap is what is printed.
 *
 * THE EQUATOR IS EXCLUDED, and not for numerical convenience: f vanishes there, so the relation is
 * singular, and geostrophy genuinely does not hold at the equator — which is why Saturn's
 * equatorial jet is not a thermal-wind feature. Columns within ATSAT_TW_LAT_TAPER degrees of the
 * equator are left out of the average entirely rather than tapered, because nothing is being
 * imposed here and a weighted average of a quantity that is being discarded would only blur it.
 *
 * A CAVEAT THAT BELONGS WITH THE NUMBER. The demanded shear uses Saturn's true radius, 58232 km.
 * The model's own metric does not: rad.z runs 1..2 and its 1/r factors are of order one, so the
 * curvature terms in RHS_Sat.cpp work on a sphere of about one length unit rather than on Saturn.
 * ATJUP had exactly this and it was fixed there (ATJUP_METRIC_RADIUS); ATSAT still has it. The
 * residual below is therefore the honest physical question — how far is this wind from balance —
 * and not a statement about what the model's own equations would call balanced.
 */

#pragma once

#include <cmath>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <iostream>

class cSaturnModel;

class ThermalWindDiagSat {
public:
    explicit ThermalWindDiagSat(cSaturnModel& model) : m(model) {}

    static constexpr double R_saturn = 5.8232e7;   // mean planetary radius [m]

    void run();

private:
    cSaturnModel& m;
};

#include "cSaturnModel.h"

inline void ThermalWindDiagSat::run(){
    std::cout << std::endl << "      ATSAT: ThermalWindDiagSat (measurement only)" << std::endl;
    auto begin = std::chrono::high_resolution_clock::now();

    const int im = m.im, jm = m.jm, km = m.km;

    static const double lat_taper_deg = [](){
        const char* e = getenv("ATSAT_TW_LAT_TAPER"); return e ? atof(e) : 10.0; }();
    const double pi_       = 3.14159265358979323846;
    const double sin_taper = std::sin(lat_taper_deg * pi_ / 180.0);

    // Metres per unit of rad.z, and metres per radial grid step.
    const double L_m  = m.L_atm * 1.0e3;
    const double dz_m = 2.0 * m.dr * L_m;          // the span of a centred difference

    double sum_act2 = 0.0, sum_dem2 = 0.0, sum_res2 = 0.0;
    double max_res = 0.0;
    long long n = 0;

    #pragma omp parallel for collapse(2) schedule(static) \
        reduction(+:sum_act2, sum_dem2, sum_res2, n) reduction(max:max_res)
    for(int j = 1; j < jm - 1; j++){
        for(int k = 1; k < km - 1; k++){

            const double colat  = m.the.z[j];
            const double costhe = std::cos(colat);          // = sin(latitude)
            if(std::fabs(costhe) < sin_taper) continue;      // equatorial band, dropped

            const double f = 2.0 * m.omega * costhe;

            for(int i = 1; i < im - 1; i++){
                // Local density from the ideal gas law with the model's own mixture constant.
                // R_mix is in J/(g K) in this model family, hence the 1e3; p_stat is in bars.
                auto rho_at = [&](int jj) -> double {
                    const double T = m.t.x[i][jj][k] * m.t_ref;
                    if(!(T > 0.0)) return 0.0;
                    return (m.p_stat.x[i][jj][k] * 1.0e5) / (m.R_mix * 1.0e3 * T);
                };
                const double rho = rho_at(j);
                if(!(rho > 0.0)) continue;

                const double drho_dthe = (rho_at(j + 1) - rho_at(j - 1)) / (2.0 * m.dthe);

                // Demanded shear [1/s] and the model's own shear [1/s].
                const double demanded = (m.g / (f * rho)) * drho_dthe / R_saturn;
                const double actual   = (m.w.x[i + 1][j][k] - m.w.x[i - 1][j][k]) * m.u_0 / dz_m;
                if(!std::isfinite(demanded) || !std::isfinite(actual)) continue;

                const double res = actual - demanded;
                sum_act2 += actual * actual;
                sum_dem2 += demanded * demanded;
                sum_res2 += res * res;
                if(std::fabs(res) > max_res) max_res = std::fabs(res);
                n++;
            }
        }
    }

    // Shears are per metre; multiplied by the shell thickness they are m/s across the whole
    // shell, which is the form a jet speed can be compared against.
    const double s = (n > 0) ? L_m / std::sqrt(double(n)) : 0.0;
    printf("      ATSAT: thermal wind — %lld cells outside the %.0f deg equatorial band;"
           " model shear %.4f, demanded %.4f, residual %.4f m/s per shell (max %.4f)\n",
           n, lat_taper_deg,
           std::sqrt(sum_act2) * s, std::sqrt(sum_dem2) * s, std::sqrt(sum_res2) * s,
           max_res * L_m);

    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for ThermalWindDiagSat\n", elapsed.count() * 1e-9);
    std::cout << "      ATSAT: ThermalWindDiagSat ended" << std::endl;
}
