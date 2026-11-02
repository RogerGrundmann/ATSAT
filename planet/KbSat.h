/*
 * Saturn Atmosphere Circulation Model (ATSAT)
 * ATSAT_KB_DIAG=<n> -- the budget of the turbulent kinetic energy k by latitude band and level, PRINT-ONLY,
 * default off; printed every n iterations, right after the closure has run (cSaturnModel.cpp).
 *
 * Written 2026-11-01 for row GIANTS-TURB-POLE (user: "split the polar k budget"): in KB1000 the k maximum
 * sits at the poles (>= 87 deg) where the vertical shear of w is ~0, and the question is which term makes it.
 *
 * Each term of rhs_tke is re-evaluated from the CURRENT fields with the differences RHS_Sat_Turb.cpp uses
 * (central, no radial stretching -- ATSAT's coord_stretching is false), so the split is the tendency the
 * next Runge-Kutta step starts from, not an RK4-weighted mean:
 *     adv      -(u dk/dr + v dk/dthe / r + w dk/dphi / (r sin))
 *     diff_r   (1/re + nue/sig_k) (d2k/dr2 + 2/r dk/dr)
 *     diff_th  (1/re + nue/sig_k) (d2k/dthe2 / r^2 + cot / r^2 dk/dthe)
 *     diff_ph  (1/re + nue/sig_k) d2k/dphi2 / (r^2 sin^2)
 *     P_k      min(prod, 20 beta* k omega)      (SST, the limited production that enters tke_source)
 *     prod     the unlimited production tensor sum (Turbulence.h compute_sources)
 *     Y_k      beta* k omega
 * k in m2/s2, tendencies in m2/s3 (nondimensional x u_0^3 / L). Means over all longitudes and both
 * hemispheres (rows j and jm-1-j). Reads only; writes nothing the model reads.
 */
#pragma once

#include <cmath>
#include <cstdio>
#include <cstdlib>

inline void cSaturnModel::kbDiag(int iter){
    static const int every = [](){ const char* e = getenv("ATSAT_KB_DIAG"); return e ? atoi(e) : 0; }();
    if(every <= 0 || !turb_active || iter % every != 0) return;
    const double sig_k = 0.85, bet_star = 0.09;
    const double L_m = L_atm * 1.0e3, tend = u_0 * u_0 * u_0 / L_m, kd = u_0 * u_0;
    const int    rows[] = { 1, 2, 3, 4, 10, 30, 60, 80, 90 };          // 89 88 87 86 80 60 30 10 0 deg
    const int    levs[] = { 20, 30, 33, 36, 38, 39 };
    const bool   phi_per = ATPhys::phi_periodic<cSaturnModel>();
    printf("\n      ATSAT: [KB] k budget at iteration %d (k m2/s2, terms m2/s3; zonal mean, N and S averaged)%s\n",
           iter, coord_stretching ? "  ** coord_stretching ON: radial terms not stretched here **" : "");
    printf("      ATSAT: [KB] %4s %3s %10s %11s %11s %11s %11s %11s %11s %11s %11s\n",
           "lat", "i", "k", "adv", "diff_r", "diff_th", "diff_ph", "P_k", "prod", "Y_k", "net");
    for(int jr : rows){
        for(int i : levs){
            if(i < 1 || i > im - 2) continue;
            double s[9] = { 0 }; int n = 0;
            for(int h = 0; h < 2; h++){
                const int j = h ? jm - 1 - jr : jr;
                if(j < 1 || j > jm - 2 || (h && j == jr)) continue;
                const double rm = metricRadius(rad.z[i]), sn = std::sin(the.z[j]), cs = std::cos(the.z[j]);
                for(int k = 0; k < km; k++){
                    const int kw = ATPhys::phi_wrap(phi_per, k - 1, km), ke = ATPhys::phi_wrap(phi_per, k + 1, km);
                    if(kw < 0 || ke >= km) continue;
                    const double K = tke.x[i][j][k];
                    const double kr = (tke.x[i+1][j][k] - tke.x[i-1][j][k]) / (2.0 * dr);
                    const double krr = (tke.x[i+1][j][k] - 2.0 * K + tke.x[i-1][j][k]) / (dr * dr);
                    const double kt = (tke.x[i][j+1][k] - tke.x[i][j-1][k]) / (2.0 * dthe);
                    const double ktt = (tke.x[i][j+1][k] - 2.0 * K + tke.x[i][j-1][k]) / (dthe * dthe);
                    const double kp = (tke.x[i][j][ke] - tke.x[i][j][kw]) / (2.0 * dphi);
                    const double kpp = (tke.x[i][j][ke] - 2.0 * K + tke.x[i][j][kw]) / (dphi * dphi);
                    const double coef = 1.0 / re + nue.x[i][j][k] / sig_k;
                    const double om = std::max(dis.x[i][j][k], 1.0e-20), ks = std::max(K, 1.0e-20);
                    const double adv = -(u.x[i][j][k] * kr + v.x[i][j][k] * kt / rm + w.x[i][j][k] * kp / (rm * sn));
                    const double d_r = coef * (krr + 2.0 / rm * kr);
                    const double d_t = coef * (ktt / (rm * rm) + cs / (rm * rm * sn) * kt);
                    const double d_p = coef * kpp / (rm * rm * sn * sn);
                    const double Pk = std::min(prod.x[i][j][k], 20.0 * bet_star * ks * om);
                    const double Yk = bet_star * ks * om;
                    const double v9[9] = { K * kd, adv * tend, d_r * tend, d_t * tend, d_p * tend, Pk * tend,
                                           prod.x[i][j][k] * tend, Yk * tend, (adv + d_r + d_t + d_p + Pk - Yk) * tend };
                    for(int q = 0; q < 9; q++) s[q] += v9[q];
                    n++;
                }
            }
            if(n == 0) continue;
            printf("      ATSAT: [KB] %4d %3d %10.4f", 90 - jr, i, s[0] / n);
            for(int q = 1; q < 9; q++) printf(" %11.3e", s[q] / n);
            printf("\n");
        }
    }
}
