/*
 * Saturn Atmosphere Circulation Model (ATSAT)
 * ATSAT_PROC_DIAG=<n> -- PROCESS fields, OUTPUT-ONLY, default 0 = off (2026-11-20, P2 of the visualisation
 * proposal; user: "start with P1, then P2").
 *
 * WHY. Every VTK file holds the STATE (vapour, cloud, ice, NH4SH, w, ...). Where NH4SH forms or decomposes,
 * where vapour condenses or evaporates, where precipitation takes condensate out and where the moist convective
 * adjustment fires could only be guessed from two snapshots. This instrument charges every change of t and
 * of the species to the stage of the time loop that made it -- as ATSAT_TB_DIAG / ATSAT_CWB_DIAG do for level
 * means and columns -- but keeps the field, zonally averaged:
 *
 *   stage          fields (zonal mean, time mean over the n iterations before the write)
 *   SatAdj         dT_satadj; cond_<g> (> 0) and evap_<g> (< 0) = cell-wise sign split of d(cloud + ice), g = H2O NH3 CH4
 *   NH4SH equil.   dT_chem; nh4sh_form (> 0), nh4sh_decomp (< 0)   (ChemistrySat::EquilibriumNH4SH, ATSAT_NH4SH_EQ)
 *   Precipitation  PrecipitationSat does not write the fields: it leaves tendencies that RHSSat applies over the
 *                  next two iterations (S_precip_*, nondimensional = kg/m3/s x L_atm[m] / u_0; Q_precip W/m3), all
 *                  times ATSAT_PRECIP_COUPLING. Sampled after each call and averaged over the calls:
 *                  conv_<g> = (S_c + S_i) u_0 / L (< 0: condensate taken out); revap_<g> = S_v u_0 / L (> 0:
 *                  re-evaporation of falling precipitation); dT_precip = Q_precip / (r_mix cp_mix), as rhs_t
 *                  forms it (without its 0.5 clamp); Pflux_<g> = rain + snow + graupel flux (Pflux_NH4SH = the
 *                  settling flux)
 *   MCA            dT_mca; dq_<g>_mca = d(vapour + cloud + ice); dNH4SH_mca; mca_freq = share of the MCA calls
 *                  that changed t in the cell (0..1)
 *
 * Rates are per PHYSICAL second: the summed change over the interval / (n * dt * L_atm * 1e3 / u_0); the
 * physics block runs on even iterations only, so its rates are already averaged over the odd ones. Units:
 * K/s; g/m3/s (species x 1e3, the scale of the slices and of printMinMax); P_* kg/m2/s as the arrays hold.
 * Means are over k = 0..km-2 (k = km-1 is the seam copy).
 *
 * Output: <output_path>/Saturn_proczm_0_<iter>.vtk, legacy ASCII structured grid with the geometry of the
 * zonal slice (DIMENSIONS jm im 1, x = 0.1 i, y = 0.05 j), plus 'height' (km) and 'PressureStat' (bar,
 * zonal mean at the write), so ParaView and satchk/tools/render.py read it as a zonal slice. ~2 MB per file.
 *
 * Cost when off: one branch per call. When on: 5 full-field snapshots (~110 MB at 41 x 181 x 361) and one
 * pass per stage. It only READS the model, so the prognostic fields are byte-identical with and without it.
 * The accumulators are not in the restart: a resumed run starts a fresh interval.
 */

#pragma once
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

namespace SatProc {

inline int stride(){
    static const int v = [](){ const char* e = getenv("ATSAT_PROC_DIAG"); return e ? atoi(e) : 0; }();
    return v;
}

enum Field {
    F_DT_SATADJ = 0, F_COND_H2O, F_EVAP_H2O, F_COND_NH3, F_EVAP_NH3, F_COND_CH4, F_EVAP_CH4,
    F_DT_CHEM, F_NH4SH_FORM, F_NH4SH_DECOMP,
    F_DT_PRECIP, F_CONV_H2O, F_CONV_NH3, F_CONV_CH4, F_REVAP_H2O, F_REVAP_NH3, F_REVAP_CH4,
    F_PFLUX_H2O, F_PFLUX_NH3, F_PFLUX_CH4, F_PFLUX_NH4SH,
    F_DT_MCA, F_DQ_H2O_MCA, F_DQ_NH3_MCA, F_DQ_CH4_MCA, F_DNH4SH_MCA, F_MCA_FREQ,
    NF
};
static const char* const NAME[NF] = {
    "dT_satadj", "cond_H2O", "evap_H2O", "cond_NH3", "evap_NH3", "cond_CH4", "evap_CH4",
    "dT_chem", "nh4sh_form", "nh4sh_decomp",
    "dT_precip", "conv_H2O", "conv_NH3", "conv_CH4", "revap_H2O", "revap_NH3", "revap_CH4",
    "Pflux_H2O", "Pflux_NH3", "Pflux_CH4", "Pflux_NH4SH",
    "dT_mca", "dq_H2O_mca", "dq_NH3_mca", "dq_CH4_mca", "dNH4SH_mca", "mca_freq" };

template<class M>
class Diag {
public:
    explicit Diag(M& m) : m(m), on(stride() > 0), n_it(0), n_precip(0), n_mca(0) {
        if(!on) return;
        n3 = (size_t)m.im * m.jm * m.km;
        acc.assign((size_t)NF * m.im * m.jm, 0.0);
        for(int s = 0; s < NS; s++) snap[s].assign(n3, 0.0);
        printf("\n      ATSAT: [PROC] ATSAT_PROC_DIAG=%d: process fields, zonal means every %d iterations -> %s/Saturn_proczm_0_<iter>.vtk\n",
               stride(), stride(), m.output_path.c_str());
    }
    bool active() const { return on; }

    // ---- SatAdj: before the three SaturationAdjustmentSat calls / after NH3 (before the NH4SH equilibrium)
    void satadj_begin(){ if(!on) return; take(S_T, T_); take(S_C1, CND_H2O); take(S_C2, CND_NH3); take(S_C3, CND_CH4); }
    void satadj_end(){
        if(!on) return;
        diff_T(S_T, F_DT_SATADJ);
        diff_split(S_C1, CND_H2O, F_COND_H2O, F_EVAP_H2O);
        diff_split(S_C2, CND_NH3, F_COND_NH3, F_EVAP_NH3);
        diff_split(S_C3, CND_CH4, F_COND_CH4, F_EVAP_CH4);
    }
    // ---- NH4SH equilibrium
    void chem_begin(){ if(!on) return; take(S_T, T_); take(S_N, NH4SH_); }
    void chem_end(){ if(!on) return; diff_T(S_T, F_DT_CHEM); diff_split(S_N, NH4SH_, F_NH4SH_FORM, F_NH4SH_DECOMP); }
    // ---- Precipitation
    void precip_begin(){}
    void precip_end(){
        if(!on) return;
        add_mean(QP_T, F_DT_PRECIP);
        add_mean(CV_H2O, F_CONV_H2O); add_mean(CV_NH3, F_CONV_NH3); add_mean(CV_CH4, F_CONV_CH4);
        add_mean(RV_H2O, F_REVAP_H2O); add_mean(RV_NH3, F_REVAP_NH3); add_mean(RV_CH4, F_REVAP_CH4);
        add_mean(PF_H2O, F_PFLUX_H2O); add_mean(PF_NH3, F_PFLUX_NH3); add_mean(PF_CH4, F_PFLUX_CH4); add_mean(PF_NH4SH, F_PFLUX_NH4SH);
        n_precip++;
    }
    // ---- MCA
    void mca_begin(){
        if(!on) return;
        take(S_T, T_); take(S_C1, TOT_H2O); take(S_C2, TOT_NH3); take(S_C3, TOT_CH4); take(S_N, NH4SH_);
    }
    void mca_end(){
        if(!on) return;
        diff_T(S_T, F_DT_MCA, F_MCA_FREQ);
        diff_net(S_C1, TOT_H2O, F_DQ_H2O_MCA); diff_net(S_C2, TOT_NH3, F_DQ_NH3_MCA); diff_net(S_C3, TOT_CH4, F_DQ_CH4_MCA);
        diff_net(S_N, NH4SH_, F_DNH4SH_MCA);
        n_mca++;
    }
    // ---- end of every iteration: count it; write + reset at the stride
    void iteration_done(int iter_n){
        if(!on) return;
        n_it++;
        if(iter_n % stride() == 0 && n_it > 0){ write(iter_n); reset(); }
    }

private:
    M& m;
    bool on;
    size_t n3;
    int n_it, n_precip, n_mca;
    std::vector<double> acc;                         // [field][i][j], sums over k and over the calls
    enum { S_T = 0, S_C1, S_C2, S_C3, S_N, NS };
    std::vector<double> snap[NS];
    enum Q { T_, CND_H2O, CND_NH3, CND_CH4, VAP_H2O, VAP_NH3, VAP_CH4, TOT_H2O, TOT_NH3, TOT_CH4, NH4SH_,
             PF_H2O, PF_NH3, PF_CH4, PF_NH4SH, QP_T, CV_H2O, CV_NH3, CV_CH4, RV_H2O, RV_NH3, RV_CH4 };
    static double coupling(){
        static const double v = [](){ const char* e = getenv("ATSAT_PRECIP_COUPLING"); return e ? atof(e) : 1.0; }();
        return v;
    }

    double q(Q w, int i, int j, int k) const {
        switch(w){
        case T_:       return m.t.x[i][j][k] * m.t_ref;                                       // K
        case CND_H2O:  return 1e3 * (m.h2o_cloud.x[i][j][k] + m.h2o_ice.x[i][j][k]);          // g/m3
        case CND_NH3:  return 1e3 * (m.nh3_cloud.x[i][j][k] + m.nh3_ice.x[i][j][k]);
        case CND_CH4:  return 1e3 * (m.ch4_cloud.x[i][j][k] + m.ch4_ice.x[i][j][k]);
        case VAP_H2O:  return 1e3 * m.h2o.x[i][j][k];
        case VAP_NH3:  return 1e3 * m.nh3.x[i][j][k];
        case VAP_CH4:  return 1e3 * m.ch4.x[i][j][k];
        case TOT_H2O:  return 1e3 * (m.h2o.x[i][j][k] + m.h2o_cloud.x[i][j][k] + m.h2o_ice.x[i][j][k]);
        case TOT_NH3:  return 1e3 * (m.nh3.x[i][j][k] + m.nh3_cloud.x[i][j][k] + m.nh3_ice.x[i][j][k]);
        case TOT_CH4:  return 1e3 * (m.ch4.x[i][j][k] + m.ch4_cloud.x[i][j][k] + m.ch4_ice.x[i][j][k]);
        case NH4SH_:   return 1e3 * m.nh4sh.x[i][j][k];
        case PF_H2O:   return m.P_rain.x[i][j][k] + m.P_snow.x[i][j][k] + m.P_graupel.x[i][j][k];    // kg/m2/s
        case PF_NH3:   return m.P_nh3_rain.x[i][j][k] + m.P_nh3_snow.x[i][j][k] + m.P_nh3_graupel.x[i][j][k];
        case PF_CH4:   return m.P_ch4_rain.x[i][j][k] + m.P_ch4_snow.x[i][j][k] + m.P_ch4_graupel.x[i][j][k];
        case PF_NH4SH: return m.P_nh4sh.x[i][j][k];
        case QP_T:     return (m.r_mix > 0.0 && m.cp_mix > 0.0) ? coupling() * m.Q_precip.x[i][j][k] / (m.r_mix * m.cp_mix) : 0.0;  // K/s
        case CV_H2O:   return 1e3 * coupling() * (m.S_precip_h2o_cloud.x[i][j][k] + m.S_precip_h2o_ice.x[i][j][k]) * m.u_0 / (m.L_atm * 1e3);  // g/m3/s
        case CV_NH3:   return 1e3 * coupling() * (m.S_precip_nh3_cloud.x[i][j][k] + m.S_precip_nh3_ice.x[i][j][k]) * m.u_0 / (m.L_atm * 1e3);
        case CV_CH4:   return 1e3 * coupling() * (m.S_precip_ch4_cloud.x[i][j][k] + m.S_precip_ch4_ice.x[i][j][k]) * m.u_0 / (m.L_atm * 1e3);
        case RV_H2O:   return 1e3 * coupling() * m.S_precip_h2o.x[i][j][k] * m.u_0 / (m.L_atm * 1e3);
        case RV_NH3:   return 1e3 * coupling() * m.S_precip_nh3.x[i][j][k] * m.u_0 / (m.L_atm * 1e3);
        case RV_CH4:   return 1e3 * coupling() * m.S_precip_ch4.x[i][j][k] * m.u_0 / (m.L_atm * 1e3);
        }
        return 0.0;
    }
    size_t id(int i, int j, int k) const { return ((size_t)i * m.jm + j) * m.km + k; }
    double& a(int f, int i, int j){ return acc[((size_t)f * m.im + i) * m.jm + j]; }

    void take(int s, Q w){
        std::vector<double>& b = snap[s];
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 0; i < m.im; i++) for(int j = 0; j < m.jm; j++) for(int k = 0; k < m.km; k++) b[id(i, j, k)] = q(w, i, j, k);
    }
    // each (i, j) is owned by one thread and sums its k in order: deterministic at any thread count
    void diff_T(int s, int f, int f_freq = -1){
        const std::vector<double>& b = snap[s];
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 0; i < m.im; i++) for(int j = 0; j < m.jm; j++){
            double d = 0.0, c = 0.0;
            for(int k = 0; k < m.km - 1; k++){ const double x = q(T_, i, j, k) - b[id(i, j, k)]; d += x; if(x != 0.0) c += 1.0; }
            a(f, i, j) += d; if(f_freq >= 0) a(f_freq, i, j) += c;
        }
    }
    void diff_net(int s, Q w, int f){
        const std::vector<double>& b = snap[s];
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 0; i < m.im; i++) for(int j = 0; j < m.jm; j++){
            double d = 0.0;
            for(int k = 0; k < m.km - 1; k++) d += q(w, i, j, k) - b[id(i, j, k)];
            a(f, i, j) += d;
        }
    }
    void diff_split(int s, Q w, int f_pos, int f_neg){
        const std::vector<double>& b = snap[s];
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 0; i < m.im; i++) for(int j = 0; j < m.jm; j++){
            double p = 0.0, n = 0.0;
            for(int k = 0; k < m.km - 1; k++){ const double x = q(w, i, j, k) - b[id(i, j, k)]; if(x > 0.0) p += x; else n += x; }
            a(f_pos, i, j) += p; a(f_neg, i, j) += n;
        }
    }
    void add_mean(Q w, int f){
        #pragma omp parallel for collapse(2) schedule(static)
        for(int i = 0; i < m.im; i++) for(int j = 0; j < m.jm; j++){
            double d = 0.0;
            for(int k = 0; k < m.km - 1; k++) d += q(w, i, j, k);
            a(f, i, j) += d;
        }
    }
    void reset(){ std::fill(acc.begin(), acc.end(), 0.0); n_it = n_precip = n_mca = 0; }

    void write(int iter_n){
        const double s_it = m.dt * m.L_atm * 1.0e3 / m.u_0;     // physical seconds per iteration
        const double t_int = n_it * s_it, nk = m.km - 1;
        const std::string name = m.output_path + "/Saturn_proczm_0_" + std::to_string(iter_n) + ".vtk";
        std::ofstream f(name);
        if(!f.is_open()){ fprintf(stderr, "ATSAT: [PROC] cannot open %s\n", name.c_str()); return; }
        char buf[64];
        f << "# vtk DataFile Version 3.0\n";
        f << "Process_zonal_mean_Sat iterations " << (iter_n - n_it + 1) << ".." << iter_n << " (" << t_int << " s)\n";
        f << "ASCII\nDATASET STRUCTURED_GRID\nDIMENSIONS " << m.jm << " " << m.im << " 1\nPOINTS " << m.im * m.jm << " float\n";
        for(int i = 0; i < m.im; i++) for(int j = 0; j < m.jm; j++){ snprintf(buf, sizeof buf, "%.4f %.4f 0.0\n", 0.1 * i, 0.05 * j); f << buf; }
        f << "POINT_DATA " << m.im * m.jm << "\n";
        auto field = [&](const char* nm, const std::function<double(int, int)>& val){
            f << "SCALARS " << nm << " float 1\nLOOKUP_TABLE default\n";
            for(int i = 0; i < m.im; i++) for(int j = 0; j < m.jm; j++){ snprintf(buf, sizeof buf, "%.6e\n", val(i, j)); f << buf; }
        };
        field("height", [&](int i, int) -> double { return m.get_layer_height(i); });
        field("PressureStat", [&](int i, int j) -> double { double s = 0.0; for(int k = 0; k < m.km - 1; k++) s += m.p_stat.x[i][j][k]; return s / nk; });
        for(int fi = 0; fi < NF; fi++){
            double norm;
            if(fi >= F_DT_PRECIP && fi <= F_PFLUX_NH4SH) norm = n_precip > 0 ? n_precip * nk : 0.0;   // mean over the calls
            else if(fi == F_MCA_FREQ)                    norm = n_mca > 0 ? n_mca * nk : 0.0;          // share of MCA calls
            else                                         norm = t_int * nk;                            // rate per second
            field(NAME[fi], [&](int i, int j) -> double { return norm > 0.0 ? acc[((size_t)fi * m.im + i) * m.jm + j] / norm : 0.0; });
        }
        f.close();
        printf("      ATSAT: [PROC] wrote %s (%d iterations = %.1f s, %d precipitation calls, %d MCA calls)\n",
               name.c_str(), n_it, t_int, n_precip, n_mca);
    }
};

} // namespace SatProc
