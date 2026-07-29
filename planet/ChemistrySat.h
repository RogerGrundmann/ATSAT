#pragma once

#include "cSaturnModel.h"

#include <cmath>
#include <chrono>
#include <cstdio>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

class ChemistrySat {
public:
    explicit ChemistrySat(cSaturnModel& model) : m(model) {}

    // -----------------------------------------------------------------------
    void ChemMassRateSat()
    {
        using namespace std;
        cout << endl << "      ATSAT: ChemMassRateSat" << endl;

        auto begin = chrono::high_resolution_clock::now();

        const int im = m.im, jm = m.jm, km = m.km;

        // reaction constants are the same for all cells
        const double keq = m.m_nh4sh / (m.m_nh3 * m.m_h2s);
        const double A   = 1500.0;
        const double B   = -0.3;
        const double T_d = 3020.0;

        #pragma omp parallel for collapse(3) schedule(static)
        for(int k = 1; k < km-1; k++){
            for(int j = 1; j < jm-1; j++){
                for(int i = 1; i < im-1; i++){
                    const double t_u = m.t.x[i][j][k] * m.t_ref;
                    double kf = react_rate_const(t_u, T_d, A, B);
                    double kb = kf / keq;

                    if((t_u <= m.t_0_nh4sh) && (t_u >= m.t_00_nh4sh))
                        react_rate_ordinary(i, j, k, kf, kb,
                            m.nh3, m.h2s, m.nh4sh, m.w_nh3, m.w_h2s, m.w_nh4sh);

                    m.massflux_h2s.x[i][j][k]   = m.w_h2s.x[i][j][k]   - m.difflux_h2s.x[i][j][k];
                    m.massflux_nh3.x[i][j][k]   = m.w_nh3.x[i][j][k]   - m.difflux_nh3.x[i][j][k];
                    m.massflux_nh4sh.x[i][j][k] = m.w_nh4sh.x[i][j][k] - m.difflux_nh4sh.x[i][j][k];
                }
            }
        }

        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for ChemMassRateSat\n", elapsed.count() * 1e-9);
        cout << "      ATSAT: ChemMassRateSat ended" << endl;
    }

    // -----------------------------------------------------------------------
    // TVD flux limiter for NH4SH advection (Superbee by default).
    // Total Variation Diminishing (TVD)
    // Computes the correction:
    //   fluxlim_nh4sh = transport_centered - transport_TVD
    //
    // Adding this to rhs_nh4sh in RHS_Sat replaces the centered-difference
    // advection with the Superbee-limited upwind scheme, preventing spurious
    // oscillations at the sharp NH4SH cloud-formation boundary.
    //
    // To switch limiter: replace superbee_phi() with van_leer_phi() below.
    // -----------------------------------------------------------------------
    void FluxLimiterNH4SH()
    {
        using namespace std;
        cout << endl << "      ATSAT: FluxLimiterNH4SH" << endl;

        auto begin = chrono::high_resolution_clock::now();

        const int    im = m.im, jm = m.jm, km = m.km;
        const double dr   = m.dr;
        const double dthe = m.dthe;
        const double dphi = m.dphi;
        constexpr double sinthe_min = 0.4;
        constexpr double eps = 1.0e-12;

        #pragma omp parallel for collapse(3) schedule(static)
        for(int k = 1; k < km-1; k++){
            for(int j = 1; j < jm-1; j++){
                for(int i = 1; i < im-1; i++){

                    const double q    = m.nh4sh.x[i][j][k];
                    const double q_rp = m.nh4sh.x[i+1][j][k];
                    const double q_rm = m.nh4sh.x[i-1][j][k];
                    const double q_tp = m.nh4sh.x[i][j+1][k];
                    const double q_tm = m.nh4sh.x[i][j-1][k];
                    const double q_pp = m.nh4sh.x[i][j][k+1];
                    const double q_pm = m.nh4sh.x[i][j][k-1];

                    const double u = m.u.x[i][j][k];
                    const double v = m.v.x[i][j][k];
                    const double w = m.w.x[i][j][k];

                    const double rm           = m.metricRadius(m.rad.z[i]);
                    const double sinthe       = max(sinthe_min, abs(sin(m.the.z[j])));
                    const double inv_rm       = 1.0 / rm;
                    const double inv_rmsinthe = 1.0 / (rm * sinthe);

                    double corr = 0.0;

                    // ---- r-direction ----
                    {
                        const double df    = q_rp - q;
                        const double db    = q    - q_rm;
                        const double denom = df + (df >= 0.0 ? eps : -eps);
                        const double r = (u >= 0.0)
                            ? db / denom
                            : ((i+2 < im ? m.nh4sh.x[i+2][j][k] : q_rp) - q_rp) / denom;
                        corr += (1.0 - superbee_phi(r)) * abs(u) * (df - db) / (2.0 * dr);
                    }

                    // ---- theta-direction ----
                    // Pole-symmetric handling: at j = 1 with v >= 0, q_tm = q[i][0][k]
                    // is the Neumann-extrapolated boundary value (c43*q - c13*q_tp),
                    // giving db = df/3 -> r = 1/3 (limiter partially active).
                    // The mirror case at j = jm-2 with v < 0 wants q[i][jm][k], which
                    // is off-grid; mirror the same Neumann extrapolation here so the
                    // limiter behaves symmetrically across the equator instead of
                    // collapsing to r = 0 (full antidiffusion) only at the south pole.
                    {
                        const double df    = q_tp - q;
                        const double db    = q    - q_tm;
                        const double denom = df + (df >= 0.0 ? eps : -eps);
                        const double q_far = (j+2 < jm)
                            ? m.nh4sh.x[i][j+2][k]
                            : (m.c43 * q_tp - m.c13 * q);   // Neumann extrap of q[jm]
                        const double r = (v >= 0.0)
                            ? db / denom
                            : (q_far - q_tp) / denom;
                        corr += (1.0 - superbee_phi(r)) * abs(v) * inv_rm * (df - db) / (2.0 * dthe);
                    }

                    // ---- phi-direction ----
                    {
                        const double df    = q_pp - q;
                        const double db    = q    - q_pm;
                        const double denom = df + (df >= 0.0 ? eps : -eps);
                        const double r = (w >= 0.0)
                            ? db / denom
                            : ((k+2 < km ? m.nh4sh.x[i][j][k+2] : q_pp) - q_pp) / denom;
                        corr += (1.0 - superbee_phi(r)) * abs(w) * inv_rmsinthe * (df - db) / (2.0 * dphi);
                    }

                    m.fluxlim_nh4sh.x[i][j][k] = corr;
                }
            }
        }

        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for FluxLimiterNH4SH\n", elapsed.count() * 1e-9);
        cout << "      ATSAT: FluxLimiterNH4SH ended" << endl;
    }

    // -----------------------------------------------------------------------
    void DiffMassFluxSat()
    {
        using namespace std;
        cout << endl << "      ATSAT: DiffMassFluxSat" << endl;

        auto begin = chrono::high_resolution_clock::now();

        const int im = m.im, jm = m.jm, km = m.km;

        // Diffusion coefficients are spatially constant — hoist out of the loop.
        m.D_nh3    = m.mue_nh3   / (m.rg_nh3   * m.sc_nh3);
        m.D_h2s    = m.mue_h2s   / (m.rg_h2s   * m.sc_h2s);
        m.D_nh4sh  = m.mue_nh4sh / (m.rg_nh4sh * m.sc_nh4sh);
        m.DT_nh3   = m.D_nh3;
        m.DT_h2s   = m.D_h2s;
        m.DT_nh4sh = m.D_nh4sh;

        // Pass 1: compute jT_* and j_* (ordinary + thermal diffusion fluxes).
        // Writes only at [i][j][k]; reads neighbouring cells of t/h2s/nh3/nh4sh.
        #pragma omp parallel for collapse(3) schedule(static)
        for(int k = 1; k < km-1; k++){
            for(int j = 1; j < jm-1; j++){
                for(int i = 1; i < im-1; i++){
                    const double rm       = m.metricRadius(m.rad.z[i]);
                    const double sinthe   = sin(m.the.z[j]);
                    const double rmsinthe = rm * sinthe;

                    double dtdr = 0.0,     dtdthe = 0.0,     dtdphi = 0.0;
                    double dh2sdr = 0.0,   dh2sdthe = 0.0,   dh2sdphi = 0.0;
                    double dnh3dr = 0.0,   dnh3dthe = 0.0,   dnh3dphi = 0.0;
                    double dnh4shdr = 0.0, dnh4shdthe = 0.0, dnh4shdphi = 0.0;

                    derivative_1_order(i, j, k, dtdr,     dtdthe,     dtdphi,     m.t);
                    derivative_1_order(i, j, k, dh2sdr,   dh2sdthe,   dh2sdphi,   m.h2s);
                    derivative_1_order(i, j, k, dnh3dr,   dnh3dthe,   dnh3dphi,   m.nh3);
                    derivative_1_order(i, j, k, dnh4shdr, dnh4shdthe, dnh4shdphi, m.nh4sh);

                    m.jT_nh3.x[i][j][k]   = m.DT_nh3   / m.t.x[i][j][k] * dtdr;
                    m.jT_h2s.x[i][j][k]   = m.DT_h2s   / m.t.x[i][j][k] * dtdr;
                    m.jT_nh4sh.x[i][j][k] = m.DT_nh4sh / m.t.x[i][j][k] * dtdr;

                    // Pole-symmetric sum of gradient components: dXdr and dXdphi are
                    // pole-symmetric for symmetric inputs, dXdthe is pole-antisymmetric.
                    // std::abs() on the theta term symmetrizes the resulting scalar.
                    const double dnh3   = dnh3dr   + std::abs(dnh3dthe)/rm   + dnh3dphi/rmsinthe;
                    const double dh2s   = dh2sdr   + std::abs(dh2sdthe)/rm   + dh2sdphi/rmsinthe;
                    const double dnh4sh = dnh4shdr + std::abs(dnh4shdthe)/rm + dnh4shdphi/rmsinthe;

                    m.j_nh3.x[i][j][k]   = m.c_mix/m.r_mix * (m.m_nh3   * m.D_nh3   * dnh3   - m.jT_nh3.x[i][j][k]);
                    m.j_h2s.x[i][j][k]   = m.c_mix/m.r_mix * (m.m_h2s   * m.D_h2s   * dh2s   - m.jT_h2s.x[i][j][k]);
                    m.j_nh4sh.x[i][j][k] = m.c_mix/m.r_mix * (m.m_nh4sh * m.D_nh4sh * dnh4sh - m.jT_nh4sh.x[i][j][k]);
                }
            }
        }

        // Pass 2: compute difflux_* (divergence of j_*) and thermalmassflux.
        // Reads j_* from neighbouring cells — all written by pass 1.
        #pragma omp parallel for collapse(3) schedule(static)
        for(int k = 1; k < km-1; k++){
            for(int j = 1; j < jm-1; j++){
                for(int i = 1; i < im-1; i++){
                    const double rm       = m.metricRadius(m.rad.z[i]);
                    const double sinthe   = sin(m.the.z[j]);
                    const double rmsinthe = rm * sinthe;

                    double dtdr = 0.0,       dtdthe = 0.0,       dtdphi = 0.0;
                    double dj_h2sdr = 0.0,   dj_h2sdthe = 0.0,   dj_h2sdphi = 0.0;
                    double dj_nh3dr = 0.0,   dj_nh3dthe = 0.0,   dj_nh3dphi = 0.0;
                    double dj_nh4shdr = 0.0, dj_nh4shdthe = 0.0, dj_nh4shdphi = 0.0;

                    derivative_1_order(i, j, k, dtdr,       dtdthe,       dtdphi,       m.t);
                    derivative_1_order(i, j, k, dj_h2sdr,   dj_h2sdthe,   dj_h2sdphi,   m.j_h2s);
                    derivative_1_order(i, j, k, dj_nh3dr,   dj_nh3dthe,   dj_nh3dphi,   m.j_nh3);
                    derivative_1_order(i, j, k, dj_nh4shdr, dj_nh4shdthe, dj_nh4shdphi, m.j_nh4sh);

                    m.difflux_h2s.x[i][j][k]   = dj_h2sdr   + std::abs(dj_h2sdthe)/rm   + dj_h2sdphi/rmsinthe;
                    m.difflux_nh3.x[i][j][k]   = dj_nh3dr   + std::abs(dj_nh3dthe)/rm   + dj_nh3dphi/rmsinthe;
                    m.difflux_nh4sh.x[i][j][k] = dj_nh4shdr + std::abs(dj_nh4shdthe)/rm + dj_nh4shdphi/rmsinthe;

                    m.thermalmassflux.x[i][j][k] =
                         (m.j_nh3.x[i][j][k]   * m.cp_nh3
                        + m.j_h2s.x[i][j][k]   * m.cp_h2s
                        + m.j_nh4sh.x[i][j][k] * m.cp_nh4sh)
                        * (dtdr + std::abs(dtdthe)/rm + dtdphi/rmsinthe)
                        + m.t.x[i][j][k] * (m.w_nh3.x[i][j][k]   * m.k_nh3
                                          + m.w_h2s.x[i][j][k]   * m.k_h2s
                                          + m.w_nh4sh.x[i][j][k] * m.k_nh4sh);
                }
            }
        }

        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for DiffMassFluxSat\n", elapsed.count() * 1e-9);
        cout << "      ATSAT: DiffMassFluxSat ended" << endl;
    }

    // -----------------------------------------------------------------------
    void ThermalPropertiesSat()
    {
        using namespace std;
        cout << endl << "      ATSAT: ThermalPropertiesSat" << endl;

        auto begin = chrono::high_resolution_clock::now();

        double r_mixture = 1.326;

        m.rg_mix = m.rg_h2 + m.rg_he + m.rg_h2s + m.rg_nh3 + m.rg_nh4sh + m.rg_h2o + m.rg_ch4;
        m.r_mix  = m.r_h2  + m.r_he  + m.r_h2s  + m.r_nh3  + m.r_nh4sh  + m.r_h2o  + m.r_ch4;
        m.c_mix  = m.c_h2  + m.c_he  + m.c_h2s  + m.c_nh3  + m.c_nh4sh  + m.c_h2o  + m.c_ch4;

        double M_mix = m.r_mix / m.c_mix;

        m.cp_mix  = m.r_h2/m.m_h2   * m.cp_h2   + m.r_h2s/m.m_h2s * m.cp_h2s
                  + m.r_he/m.m_he   * m.cp_he   + m.r_nh3/m.m_nh3 * m.cp_nh3
                  + m.r_nh4sh/m.m_nh4sh * m.cp_nh4sh + m.r_h2o/m.m_h2o * m.cp_h2o
                  + m.r_ch4/m.m_ch4 * m.cp_ch4;
        m.mue_mix = m.r_h2/m.m_h2   * m.mue_h2  + m.r_h2s/m.m_h2s * m.mue_h2s
                  + m.r_he/m.m_he   * m.mue_he  + m.r_nh3/m.m_nh3 * m.mue_nh3
                  + m.r_nh4sh/m.m_nh4sh * m.mue_nh4sh + m.r_h2o/m.m_h2o * m.mue_h2o
                  + m.r_ch4/m.m_ch4 * m.mue_ch4;
        m.k_mix   = m.r_h2/m.m_h2   * m.k_h2    + m.r_h2s/m.m_h2s * m.k_h2s
                  + m.r_he/m.m_he   * m.k_he    + m.r_nh3/m.m_nh3 * m.k_nh3
                  + m.r_nh4sh/m.m_nh4sh * m.k_nh4sh  + m.r_h2o/m.m_h2o * m.k_h2o
                  + m.r_ch4/m.m_ch4 * m.k_ch4;
        m.R_mix   = (m.r_h2 * m.R_h2 + m.r_he * m.R_he + m.r_h2o * m.R_h2o
                  + m.r_h2s * m.R_h2s + m.r_nh3 * m.R_nh3 + m.r_nh4sh * m.R_nh4sh
                  + m.r_ch4 * m.R_ch4)
                  / (m.r_h2 + m.r_he + m.r_h2o + m.r_h2s + m.r_nh3 + m.r_nh4sh + m.r_ch4);

        cout.precision(10);
        cout.setf(ios::fixed);
        cout << endl
            << "     r_mixture[kg/m³] = " << r_mixture << endl << endl
            << "     rg_mix[kg/m³] = " << m.rg_mix << endl << endl
            << "     r_mix[kg/m³] = " << m.r_mix << endl << endl
            << "     c_mix[kmol/m³] = " << m.c_mix << endl << endl
            << "     M_mix[kg/Kmol] = " << M_mix << endl << endl
            << "     cp_mix[J/(kg*K)] = " << m.cp_mix << endl
            << "     mue_mix[Ns/m²] = " << m.mue_mix << endl
            << "     k_mix[W/(m*K)] = " << m.k_mix << endl
            << "     R_mix[J/(kg*K)] = " << m.R_mix << endl << endl;

        auto end = chrono::high_resolution_clock::now();
        auto elapsed = chrono::duration_cast<chrono::nanoseconds>(end - begin);
        printf(" time measured: %.3f seconds for ThermalPropertiesSat\n", elapsed.count() * 1e-9);
        cout << "      ATSAT: ThermalPropertiesSat ended" << endl;
    }

private:
    cSaturnModel& m;

    // -----------------------------------------------------------------------
    void derivative_1_order(int i, int j, int k,
        double &dcdr, double &dcdthe, double &dcdphi, Array &c)
    {
        // radial gradients — SeaMount-aware
        if((m.SeaMount.x[i][j][k]==1.0) && (m.SeaMount.x[i+1][j][k]==0.0) && (m.SeaMount.x[i+2][j][k]==0.0)){
            c.x[i][j][k] = c.x[i+3][j][k] - 3.0*c.x[i+2][j][k] + 3.0*c.x[i+1][j][k];
            dcdr = (-3.0*c.x[i][j][k] + 4.0*c.x[i+1][j][k] - c.x[i+2][j][k]) / (2.0*m.dr);
        }
        if((m.SeaMount.x[i-1][j][k]==1.0) && (m.SeaMount.x[i][j][k]==0.0) && (m.SeaMount.x[i+1][j][k]==0.0))
            dcdr = (c.x[i+1][j][k] - c.x[i-1][j][k]) / (2.0*m.dr);
        if((m.SeaMount.x[i-1][j][k]==0.0) && (m.SeaMount.x[i][j][k]==0.0) && (m.SeaMount.x[i+1][j][k]==0.0))
            dcdr = (c.x[i+1][j][k] - c.x[i-1][j][k]) / (2.0*m.dr);

        // south (j+1) gradients
        if((m.SeaMount.x[i][j][k]==1.0) && (m.SeaMount.x[i][j+1][k]==0.0) && (m.SeaMount.x[i][j+2][k]==0.0)){
            c.x[i][j][k] = c43*c.x[i][j+1][k] - c13*c.x[i][j+2][k];
            dcdthe = (-3.0*c.x[i][j][k] + 4.0*c.x[i][j+1][k] - c.x[i][j+2][k]) / (2.0*m.dthe);
        }
        if((m.SeaMount.x[i][j-1][k]==1.0) && (m.SeaMount.x[i][j][k]==0.0) && (m.SeaMount.x[i][j+1][k]==0.0))
            dcdthe = (c.x[i][j+1][k] - c.x[i][j-1][k]) / (2.0*m.dthe);

        // north (j-1) gradients
        if((m.SeaMount.x[i][j][k]==1.0) && (m.SeaMount.x[i][j-1][k]==0.0) && (m.SeaMount.x[i][j-2][k]==0.0)){
            c.x[i][j][k] = c43*c.x[i][j-1][k] - c13*c.x[i][j-2][k];
            dcdthe = (-3.0*c.x[i][j][k] + 4.0*c.x[i][j-1][k] - c.x[i][j-2][k]) / (2.0*m.dthe);
        }
        if((m.SeaMount.x[i][j-1][k]==1.0) && (m.SeaMount.x[i][j][k]==0.0) && (m.SeaMount.x[i][j+1][k]==0.0))
            dcdthe = (c.x[i][j+1][k] - c.x[i][j-1][k]) / (2.0*m.dthe);
        if((m.SeaMount.x[i][j-1][k]==0.0) && (m.SeaMount.x[i][j][k]==0.0) && (m.SeaMount.x[i][j+1][k]==0.0))
            dcdthe = (c.x[i][j+1][k] - c.x[i][j-1][k]) / (2.0*m.dthe);

        // east (k+1) gradients
        if((m.SeaMount.x[i][j][k]==1.0) && (m.SeaMount.x[i][j][k+1]==0.0) && (m.SeaMount.x[i][j][k+2]==0.0)){
            c.x[i][j][k] = c43*c.x[i][j][k+1] - c13*c.x[i][j][k+2];
            dcdphi = (-3.0*c.x[i][j][k] + 4.0*c.x[i][j][k+1] - c.x[i][j][k+2]) / (2.0*m.dphi);
        }
        if((m.SeaMount.x[i][j][k-1]==1.0) && (m.SeaMount.x[i][j][k]==0.0) && (m.SeaMount.x[i][j][k+1]==0.0))
            dcdphi = (c.x[i][j][k+1] - c.x[i][j][k-1]) / (2.0*m.dphi);

        // west (k-1) gradients
        if((m.SeaMount.x[i][j][k]==1.0) && (m.SeaMount.x[i][j][k-1]==0.0) && (m.SeaMount.x[i][j][k-2]==0.0)){
            c.x[i][j][k] = c43*c.x[i][j][k-1] - c13*c.x[i][j][k-2];
            dcdphi = (-3.0*c.x[i][j][k] + 4.0*c.x[i][j][k-1] - c.x[i][j][k-2]) / (2.0*m.dphi);
        }
        if((m.SeaMount.x[i][j][k-1]==1.0) && (m.SeaMount.x[i][j][k]==0.0) && (m.SeaMount.x[i][j][k+1]==0.0))
            dcdphi = (c.x[i][j][k+1] - c.x[i][j][k-1]) / (2.0*m.dphi);
        if((m.SeaMount.x[i][j][k-1]==0.0) && (m.SeaMount.x[i][j][k]==0.0) && (m.SeaMount.x[i][j][k+1]==0.0))
            dcdphi = (c.x[i][j][k+1] - c.x[i][j][k-1]) / (2.0*m.dphi);
    }

    // -----------------------------------------------------------------------
    void derivative_1_order_boundary(int i, int j, int k,
        double &dcdr, double &dcdthe, double &dcdphi, Array &c)
    {
        if(i == 0)       dcdr   = (-3.0*c.x[0][j][k]    + 4.0*c.x[1][j][k]    - c.x[2][j][k])    / (2.0*m.dr);
        if(i == m.im-1)  dcdr   = (-3.0*c.x[m.im-1][j][k] + 4.0*c.x[m.im-2][j][k] - c.x[m.im-3][j][k]) / (2.0*m.dr);
        if(j == 0)       dcdthe = (-3.0*c.x[i][0][k]    + 4.0*c.x[i][1][k]    - c.x[i][2][k])    / (2.0*m.dthe);
        if(j == m.jm-1)  dcdthe = (-3.0*c.x[i][m.jm-1][k] + 4.0*c.x[i][m.jm-2][k] - c.x[i][m.jm-3][k]) / (2.0*m.dthe);
        if(k == 0)       dcdphi = (-3.0*c.x[i][j][0]    + 4.0*c.x[i][j][1]    - c.x[i][j][2])    / (2.0*m.dphi);
        if(k == m.km-1)  dcdphi = (-3.0*c.x[i][j][m.km-1] + 4.0*c.x[i][j][m.km-2] - c.x[i][j][m.km-3]) / (2.0*m.dphi);
    }

    // -----------------------------------------------------------------------
    void react_rate_ordinary(int i, int j, int k,
        double &f, double &b,
        Array &molecule_nh3, Array &molecule_h2s, Array &molecule_nh4sh,
        Array &mass_nh3,     Array &mass_h2s,     Array &mass_nh4sh)
    {
        double c_nh3   = molecule_nh3.x[i][j][k]   / m.m_nh3;
        double c_h2s   = molecule_h2s.x[i][j][k]   / m.m_h2s;
        double c_nh4sh = molecule_nh4sh.x[i][j][k] / m.m_nh4sh;

        double Rf = f * (c_nh3 * c_h2s) / m.dt;
        double Rb = b * c_nh4sh         / m.dt;
        double R_diff = Rf - Rb;

        mass_nh3.x[i][j][k]   = -m.m_nh3   * R_diff;
        mass_h2s.x[i][j][k]   = -m.m_h2s   * R_diff;
        mass_nh4sh.x[i][j][k] =  m.m_nh4sh * R_diff;
    }

    // -----------------------------------------------------------------------
    static double react_rate_const(double T_K, const double T_d,
        const double A, const double B)
    {
        return A * std::pow(T_K, B) * std::exp(-T_d / T_K);
    }

    // cSaturnModel::c43 / c13 shorthands
    static constexpr double c43 = 4.0/3.0;
    static constexpr double c13 = 1.0/3.0;

    // Superbee: most compressive TVD limiter — best for sharp cloud fronts.
    static double superbee_phi(double r){
        return std::max(0.0, std::max(std::min(2.0*r, 1.0), std::min(r, 2.0)));
    }

    // Van Leer: smooth, differentiable — good general-purpose alternative.
    static double van_leer_phi(double r){
        return (r + std::abs(r)) / (1.0 + std::abs(r));
    }
};
