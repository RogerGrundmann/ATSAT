/*
 * ATSAT's entry point into the saturation adjustment.
 *
 * The ALGORITHM is not here — it is the shared SaturationAdjustment<Planet> template, which
 * ATJUP instantiates too. What is here is the two things that are ATSAT's alone: the choice
 * between the inherited routine and the shared one, and the ice coefficients ATSAT does not have.
 */

#include "SaturationAdjustmentSat.h"
#include "cSaturnModel.h"

using namespace std;

// Dispatch. The legacy routine keeps every call site working unchanged and remains the DEFAULT;
// ATSAT_SATADJ=1 selects the shared algorithm. See the header for what the two differ by and
// what the 16 % it costs in peak cloud water is evidence of.
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

    // The ice quadruple, by gas. Until 2026-10-09 this passed the liquid one in its place; the
    // model's ice coefficients are real now (cSaturnModel.h has the note). A gas without an ice
    // pair of its own gets its liquid pair, which is what CH4's is anyway.
    double C_i = C, L0_i = L0, del_alf_i = del_alf, del_bet_i = del_bet;
    if(gas == "H2O"){
        C_i = m.C_h2o_ice; L0_i = m.L0_h2o_ice; del_alf_i = m.del_alf_h2o_ice; del_bet_i = m.del_bet_h2o_ice;
    } else if(gas == "NH3"){
        C_i = m.C_nh3_ice; L0_i = m.L0_nh3_ice; del_alf_i = m.del_alf_nh3_ice; del_bet_i = m.del_bet_nh3_ice;
    } else if(gas == "CH4"){
        C_i = m.C_ch4_ice; L0_i = m.L0_ch4_ice; del_alf_i = m.del_alf_ch4_ice; del_bet_i = m.del_bet_ch4_ice;
    }
    SaturationAdjustment<cSaturnModel>(m).run(gas, t_0, t_00, ep, lv, ls,
                                              C,   L0,   R, del_alf,   del_bet,
                                              C_i, L0_i,    del_alf_i, del_bet_i,
                                              c, cloud, ice);
}
