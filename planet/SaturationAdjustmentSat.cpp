/*
 * Mirrored saturation adjustment for ATSAT — the counterpart of
 * SaturationAdjustmentJup.cpp. See SaturationAdjustmentSat.h for what differs from the
 * inherited cSaturnModel::Saturation_Adjustment() in Weather_Sat.cpp, which remains the
 * default and which this file does not touch.
 */

#include "SaturationAdjustmentSat.h"
#include "cSaturnModel.h"

using namespace std;

// Dispatch. The legacy routine keeps every call site working unchanged.
void SaturationAdjustmentSat::run(std::string gas,
    double &coeff_A,   double &coeff_B,
    double &coeff_A_i, double &coeff_B_i,
    double &t_0,       double &t_00,
    double &ep,        double &lv,  double &ls,
    double &cp,        double &r,
    double &C,         double &L0,  double &R,
    double &del_alf,   double &del_bet,  double &m_mol,
    Array &c, Array &cloud, Array &ice)
{
    if(mirrored_enabled() == 0){
        m.Saturation_Adjustment(gas, coeff_A, coeff_B, coeff_A_i, coeff_B_i,
            t_0, t_00, ep, lv, ls, cp, r, C, L0, R, del_alf, del_bet, m_mol,
            c, cloud, ice);
        return;
    }

    // The ice-phase quadruple is the LIQUID one, because ATSAT's parameter set has no ice
    // pair. This is the single place to change when real values for Saturn's H2O, NH3 and CH4
    // ices exist; the algorithm below already treats the two pairs separately, so nothing else
    // has to move. Point 4 in the header note records what it costs meanwhile.
    run_mirrored(gas, t_0, t_00, ep, lv, ls,
                 C, L0, R, del_alf, del_bet,
                 C, L0,    del_alf, del_bet,
                 c, cloud, ice);
}


void SaturationAdjustmentSat::run_mirrored(const std::string& gas,
        double t_0,       double t_00,
        double ep,        double lv,  double ls,
        double C,         double L0,  double R,
        double del_alf,   double del_bet,
        double C_i,       double L0_i,
        double del_alf_i, double del_bet_i,
        Array& c,         Array& cloud,   Array& ice)
{
    cout << endl << "      SaturationAdjustment of " << gas
         << " begin (mirrored from ATJUP)" << endl;

    auto begin = std::chrono::high_resolution_clock::now();

    const double exp_pressure = m.g / (m.gam * m.R_ref);
    const double t_range_inv  = 1.0 / (t_0 - t_00);

    // Shared diagnostic state — written under omp critical, unlike the legacy routine, which
    // writes it from every thread at once (point 8 in the header note).
    bool   sat_found       = false;
    int    iter_prec_found = 0;
    int    i_sat = 0, j_sat = 0, k_sat = 0;
    double height_sat = 0.0;
    double t_latent = 0.0, p_latent = 0.0;
    double t_sat    = 0.0, p_sat    = 0.0;
    double t_u_sat  = 0.0, p_u_sat  = 0.0;
    double q_v_b_sat = 0.0, q_c_b_sat = 0.0, q_i_b_sat = 0.0;
    double saturation = 0.0;

    // -----------------------------------------------------------------------
    // Main saturation-adjustment loop — fully independent per cell
    // -----------------------------------------------------------------------
    #pragma omp parallel for collapse(2) schedule(static)
    for(int k = 0; k < m.km; k++){
        for(int j = 0; j < m.jm; j++){
            for(int i = 0; i < m.im; i++){

                const double t_u = m.t.x[i][j][k] * m.t_ref;   // [K]
                const double p_u = m.p_stat.x[i][j][k];        // [bar]

                // Cells with no usable thermodynamic state are skipped. Both formulas below
                // are singular there, and the entry test further down cannot catch a NaN
                // because every comparison against one is false. Written as !(x > 0.0) so a
                // NaN that arrived from elsewhere is skipped rather than propagated.
                // ATJUP additionally excludes its SeaMount interior; ATSAT has no solid body.
                if(!(t_u > 0.0) || !(p_u > 0.0)) continue;

                // Enforce physical bounds
                if(t_u > t_0)              ice.x[i][j][k]   = 0.0;
                if(c.x[i][j][k]     < 0.0) c.x[i][j][k]     = 0.0;
                if(cloud.x[i][j][k] < 0.0) cloud.x[i][j][k] = 0.0;
                if(ice.x[i][j][k]   < 0.0) ice.x[i][j][k]   = 0.0;

                // Density for both the saturation quantity and the latent-heat divisor.
                // rho_at() is now ATSAT's single gated accessor, exactly as ATJUP's is: r_mix
                // unless ATSAT_LOCAL_RHO is set. PrecipitationSat goes through the same gate, so
                // the microphysics and the adjustment saturate on the same density by
                // construction — which is what PrecipitationSat.h has always claimed.
                //
                // Why the gate defaults to the reference density: switching this to the local
                // one diverged the run, T to 4.65e5 degC in eight iterations. Both halves of the
                // feedback push the same way — the divisor cp_mix*rho is small in a thin cell so
                // the heating per unit condensate is large, and the target rho*ep*E/p is scaled
                // by the same small rho so more vapour reads as supersaturated.
                const double rho_c = m.rho_at(i, j, k);

                const double E_Rain_0 = saturation_vapour_pressure(t_u, C, L0, R,
                                                                   del_alf, del_bet);
                const double q_Rain_0 = rho_c * ep * E_Rain_0 / p_u;

                // Skip: subsaturated, already saturated, or the SVP underflowed to zero in a
                // very cold cell.
                if(c.x[i][j][k] <= q_Rain_0 || q_Rain_0 <= 0.0) continue;

                // ---- Mixed-phase iteration (Tao et al. 1988) ----
                double q_v_b   = c.x[i][j][k];
                double q_c_b   = cloud.x[i][j][k];
                double q_i_b   = ice.x[i][j][k];
                double T       = t_u;
                double q_v_hyp = q_v_b;

                bool cell_found = false;
                int  cell_iter  = 0;

                for(int itr = 1; itr <= iter_prec_end; itr++){
                    double CND = (T - t_00) * t_range_inv;
                    double DEP = (t_0 - T)  * t_range_inv;
                    if(T <= t_00){ CND = 0.0; DEP = 1.0; }
                    if(T >= t_0) { CND = 1.0; DEP = 0.0; }

                    const double d_q_v = q_v_hyp - q_v_b;
                    const double d_q_c = -d_q_v * CND;
                    const double d_q_i = -d_q_v * DEP;

                    T     += (lv * d_q_c + ls * d_q_i) / (m.cp_mix * rho_c);
                    q_v_b += d_q_v;
                    q_c_b += d_q_c;
                    q_i_b += d_q_i;

                    if(q_v_b < 0.0) q_v_b = 0.0;
                    if(q_c_b < 0.0) q_c_b = 0.0;
                    if(q_i_b < 0.0) q_i_b = 0.0;

                    // AT THE UPDATED TEMPERATURE T, not at the entry temperature. This is the
                    // feedback the scheme exists to resolve: latent heat raises T, which raises
                    // the saturation vapour pressure, which limits further condensation.
                    const double E_Rain = saturation_vapour_pressure(T, C,   L0,   R,
                                                                     del_alf,   del_bet);
                    const double E_Ice  = saturation_vapour_pressure(T, C_i, L0_i, R,
                                                                     del_alf_i, del_bet_i);
                    const double q_Rain = rho_c * ep * E_Rain / p_u;
                    const double q_Ice  = rho_c * ep * E_Ice  / p_u;

                    if(q_c_b > 0.0 && q_i_b > 0.0)
                        q_v_hyp = (q_c_b * q_Rain + q_i_b * q_Ice) / (q_c_b + q_i_b);
                    else if(q_i_b == 0.0) q_v_hyp = q_Rain;
                    else                  q_v_hyp = q_Ice;

                    if(T >= t_0) q_i_b = 0.0;

                    const double q_diff = std::fabs(q_v_b - q_v_hyp) / (q_v_hyp + 1e-20);
                    if(q_diff <= q_diff_min){
                        cell_found = true;
                        cell_iter  = itr;
                        break;
                    }
                    q_v_hyp = 0.5 * (q_v_hyp + q_v_b);   // has smoothing effect
                }

                // Write converged cell values back
                c.x[i][j][k]     = q_v_b;
                cloud.x[i][j][k] = q_c_b;
                ice.x[i][j][k]   = q_i_b;
                m.t.x[i][j][k]   = T / m.t_ref;

                // KEPT FROM ATSAT, NOT FROM ATJUP: the hydrostatic pressure is rewritten from
                // the adjusted temperature. ATJUP does not touch p_stat here. See the closing
                // note in the header — this is flagged, not settled.
                m.p_stat.x[i][j][k] = m.p_ref
                    * std::pow(m.t.x[i][j][k], exp_pressure);   // [bar]

                if(cell_found){
                    #pragma omp critical
                    {
                        sat_found       = true;
                        iter_prec_found = cell_iter;
                        i_sat = i; j_sat = j; k_sat = k;
                        height_sat = m.get_layer_height(i_sat);
                        t_u_sat    = t_u;
                        p_u_sat    = p_u;
                        t_sat      = T;
                        p_sat      = m.p_ref * std::pow(t_sat / m.t_ref, exp_pressure);
                        t_latent   = t_sat - t_u_sat;
                        p_latent   = p_sat - p_u_sat;
                        q_v_b_sat  = q_v_b;
                        q_c_b_sat  = q_c_b;
                        q_i_b_sat  = q_i_b;
                        saturation = q_v_b - q_Rain_0;
                    }
                }
            }
        }
    }

    // -----------------------------------------------------------------------
    // Clamp negative values left by the iteration
    // -----------------------------------------------------------------------
    #pragma omp parallel for collapse(2) schedule(static)
    for(int k = 0; k < m.km; k++){
        for(int j = 0; j < m.jm; j++){
            for(int i = 0; i < m.im; i++){
                if(c.x[i][j][k]     <= 0.0) c.x[i][j][k]     = 0.0;
                if(cloud.x[i][j][k] <= 0.0) cloud.x[i][j][k] = 0.0;
                if(ice.x[i][j][k]   <= 0.0) ice.x[i][j][k]   = 0.0;
            }
        }
    }

    auto end     = std::chrono::high_resolution_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin);
    printf(" time measured: %.3f seconds for SaturationAdjustment\n", elapsed.count() * 1e-9);

    if(!sat_found){
        cout << "      NO saturation in SaturationAdjustment of " << gas << " found" << endl
             << "      iter_prec_end = " << iter_prec_end << endl;
    } else {
        cout << "      saturation in SaturationAdjustment of " << gas << " found" << endl
             << "      iter_prec_found = " << iter_prec_found
             << "   iter_prec_end = "      << iter_prec_end   << endl
             << "      i_sat = "     << i_sat
             << "   j_sat = "        << j_sat
             << "   k_sat = "        << k_sat
             << "   height_sat[km] = " << height_sat  << endl
             << "      p_stat[bar] = "  << p_sat
             << "   p_u[bar] = "        << p_u_sat
             << "   p_latent[bar] = "   << p_latent   << endl
             << "      T[°C] = "        << t_sat    - m.t_ref
             << "   t_u[°C] = "         << t_u_sat  - m.t_ref
             << "   t_latent[°C] = "    << t_latent  << endl
             << "      saturation[g/m³] = " << saturation * 1e3 << endl
             << "      " << gas << " humid[g/m³] = "  << q_v_b_sat * 1e3
             << "   cloud[g/m³] = "  << q_c_b_sat * 1e3
             << "   ice[g/m³] = "    << q_i_b_sat * 1e3 << endl;
    }

    cout << "      SaturationAdjustment of " << gas << " ended" << endl;
}
