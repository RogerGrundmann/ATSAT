/*
 * ATSAT_CHROM=1 (2026-10-30, default 0) -- the CHROMOPHORE TRACER, ported from ATJUP's ChromJup.h (ATJUP_CHROM,
 * 2026-10-24) on the user's word ("ATSAT + CHROM tracer"). Same equation, same defaults, same wiring.
 *
 * Every cloud this model condenses (NH3 ice, NH4SH, H2O, CH4) is white. The tints of the photographs come
 * from trace amounts of coloured compounds, the chromophores, made where solar UV breaks up NH3 (with PH3,
 * H2S, C2H2) high in the atmosphere and the fragments recombine. On Saturn the colours are muted -- a
 * thicker photochemical haze sits over the NH3 deck -- but where they form is the same argument as on
 * Jupiter: in the UV-lit top, from whatever NH3 the dynamics brings up there.
 *
 * So this is a passive species chrom [kg/m3], carried and diffused like NH3 (its Schmidt number), with
 *
 *     d chrom/dt |_chem = J * c_NH3 * F_sun(theta) * A_uv  -  chrom / tau
 *
 *     c_NH3    NH3 vapour + cloud + ice (ATSAT_CHROM_SRC=1, default) or vapour alone (0)
 *     F_sun    sin(colatitude): the diurnal-mean insolation shape, no season (1 at the equator). Saturn's
 *              26.7 deg obliquity and its ring shadow are NOT in it.
 *     A_uv     exp(-N_above / N_uv): the UV left after the NH3 VAPOUR column N_above [kg/m2] above the
 *              cell (from the cell's upper face to the lid) has absorbed it -- NH3 shields itself.
 *              N_uv = ATSAT_CHROM_NUV [kg/m2], default 3e-4 (ATJUP's: tau = 1 at ~1e22 molecules/m2,
 *              the NH3 cross-section ~1e-18 cm2 at 190-220 nm).
 *              ATSAT_CHROM_HAZE_P1=<bar> (2026-11-21, default 0 = off, byte-identical): times exp(-tau_haze),
 *              the UV optical depth of a photochemical HAZE with a uniform mixing ratio, so tau_haze grows
 *              linearly in pressure: tau_haze = min(p_up, p_base) / p1, p_up the static pressure at the cell's
 *              upper face (p_stat of the cell above; the lid's own), p1 the pressure where tau_haze = 1.
 *              ATSAT_CHROM_HAZE_PBASE=<bar> (default 0 = no base) ends the haze there (Saturn's upper-
 *              tropospheric haze sits over the NH3 cloud at ~0.1-0.3 bar, Sromovsky et al. arXiv 2108.00322).
 *              The haze above the lid only adds a column-uniform factor exp(-p_lid/p1), as J does, except
 *              that p_lid varies with latitude. Not modelled: haze that depends on latitude or season.
 *     J        ATSAT_CHROM_J, effective production rate [1/s], default 1e-6 (ATJUP's)
 *
 *     ATSAT_CHROM_AURORA=<A> (2026-11-25, row SAT-CHROM-REALISM step (3); default 0 = off, byte-identical): a second
 *     source that does NOT need NH3, the AURORAL one -- Saturn's poles are UV-dark from aerosol made by auroral
 *     particle precipitation (Pryor et al. 2019), which the NH3 source (max at +-45 deg / eq, poles 30-40x weaker)
 *     cannot give. It adds
 *         + A * exp(-((|lat| - lat_a) / w)^2 / 2) * exp(-p / p_a)        [kg/m3/s]
 *     lat_a = ATSAT_CHROM_AURORA_LAT [deg], default 75 (the main oval, ~70-80 deg in both hemispheres); w =
 *     ATSAT_CHROM_AURORA_W [deg], default 5; p_a = ATSAT_CHROM_AURORA_P [bar], default 0.1, p = p_stat of the cell:
 *     largest at the lid, standing in for aerosol made far above the model top and settling in. Uncalibrated like J.
 *     tau      ATSAT_CHROM_TAU [days of 86400 s], default 30 (ATJUP's)
 *
 * A_uv is built once per iteration (chromUV(), before RungeKuttaSat) into chrom_uv.
 *
 * UNCALIBRATED ON PURPOSE, as on ATJUP: read chrom as a relative colour intensity. It takes nothing from
 * NH3 and feeds nothing back (no radiation, no density). Not modelled: haze / cloud scattering in A_uv,
 * settling, scavenging by precipitation, a temperature-dependent destruction at depth.
 *
 * Off: chrom is never touched and stays 0, no output carries it and the restart file is unchanged, so a
 * default run is byte-identical. On: the restart file carries chrom and chromn after the other arrays; a
 * restart file without them starts chrom from 0 (load_state says so).
 */

#pragma once

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

namespace SatChrom {
    struct Par { double J, N_uv, tau_s; int src; double haze_p1, haze_pbase, aur, aur_lat, aur_w, aur_p; };
    inline const Par& par(){
        static const Par p = [](){
            auto d = [](const char* n, double v){ const char* e = getenv(n); return e ? atof(e) : v; };
            return Par{ d("ATSAT_CHROM_J", 1.0e-6), d("ATSAT_CHROM_NUV", 3.0e-4), d("ATSAT_CHROM_TAU", 30.0) * 86400.0,
                        (int)d("ATSAT_CHROM_SRC", 1.0), d("ATSAT_CHROM_HAZE_P1", 0.0), d("ATSAT_CHROM_HAZE_PBASE", 0.0),
                        d("ATSAT_CHROM_AURORA", 0.0), d("ATSAT_CHROM_AURORA_LAT", 75.0), d("ATSAT_CHROM_AURORA_W", 5.0),
                        d("ATSAT_CHROM_AURORA_P", 0.1) };
        }();
        return p;
    }
}

// chrom_uv = exp(-N_above / N_uv), N_above the NH3 vapour column from the upper face of the cell to the lid
// (the lid layer im-1 has no thickness, as in every column sum here); with ATSAT_CHROM_HAZE_P1 > 0 also
// times exp(-min(p_up, p_base) / p1), p_up = p_stat of the cell above (bar).
inline void cSaturnModel::chromUV(){
    if(!chrom_on()) return;
    const double n_uv = SatChrom::par().N_uv;
    const double p1 = SatChrom::par().haze_p1, pb = SatChrom::par().haze_pbase;
    #pragma omp parallel for collapse(2) schedule(static)
    for(int j = 0; j < jm; j++){
        for(int k = 0; k < km; k++){
            double n_above = 0.0;
            chrom_uv.x[im-1][j][k] = 1.0;
            for(int i = im - 2; i >= 0; i--){
                chrom_uv.x[i][j][k] = n_uv > 0.0 ? std::exp(-n_above / n_uv) : 1.0;
                if(p1 > 0.0){
                    double p_up = std::max(0.0, p_stat.x[i+1][j][k]);
                    if(pb > 0.0) p_up = std::min(p_up, pb);
                    chrom_uv.x[i][j][k] *= std::exp(-p_up / p1);
                }
                n_above += std::max(0.0, nh3.x[i][j][k]) * layer_thickness_m(i);
            }
        }
    }
}

inline double cSaturnModel::chromSource(int i, int j, int k){
    const SatChrom::Par& p = SatChrom::par();
    const double c = nh3.x[i][j][k] + (p.src ? nh3_cloud.x[i][j][k] + nh3_ice.x[i][j][k] : 0.0);
    const double f_sun = std::max(0.0, std::sin(the.z[j]));
    double rate = p.J * std::max(0.0, c) * f_sun * chrom_uv.x[i][j][k]
                      - (p.tau_s > 0.0 ? chrom.x[i][j][k] / p.tau_s : 0.0);     // kg/m3/s
    if(p.aur > 0.0){                                                            // ATSAT_CHROM_AURORA
        const double lat = 90.0 - the.z[j] * 180.0 / M_PI, x = (std::abs(lat) - p.aur_lat) / p.aur_w;
        const double pb = std::max(0.0, p_stat.x[i][j][k]);
        rate += p.aur * std::exp(-0.5 * x * x) * (p.aur_p > 0.0 ? std::exp(-pb / p.aur_p) : 1.0);
    }
    return rate * (L_atm * 1.0e3 / u_0);                                          // per nondimensional time
}

inline void cSaturnModel::chromInfo(){
    if(!chrom_on()) return;
    const SatChrom::Par& p = SatChrom::par();
    printf("      ATSAT: ATSAT_CHROM on -- chromophore tracer: J %.3e 1/s, N_uv %.3e kg/m2 (NH3 vapour self-shielding), tau %.1f d, source %s\n",
           p.J, p.N_uv, p.tau_s / 86400.0, p.src ? "NH3 vapour + cloud + ice" : "NH3 vapour");
    if(p.haze_p1 > 0.0)
        printf("      ATSAT: ATSAT_CHROM_HAZE_P1 %.3f bar -- UV haze tau = min(p, %s) / p1 in A_uv\n",
               p.haze_p1, p.haze_pbase > 0.0 ? (std::to_string(p.haze_pbase) + " bar").c_str() : "no base");
    if(p.aur > 0.0)
        printf("      ATSAT: ATSAT_CHROM_AURORA %.3e kg/m3/s -- auroral source at +-%.1f deg (w %.1f deg), exp(-p / %.3f bar)\n",
               p.aur, p.aur_lat, p.aur_w, p.aur_p);
}
