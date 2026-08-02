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

    // THE ICE QUADRUPLE IS THE LIQUID ONE. ATSAT's parameter set has no ice pair for H2O, NH3 or
    // CH4, and inventing numbers for Saturn's ices is a physics decision rather than a port. This
    // is the single place to change when real values exist: the shared algorithm already treats
    // the liquid and ice pairs separately, so nothing else has to move. Point 4 in the header
    // note records what the substitution costs meanwhile.
    SaturationAdjustment<cSaturnModel>(m).run(gas, t_0, t_00, ep, lv, ls,
                                              C, L0, R, del_alf,   del_bet,
                                              C, L0,    del_alf,   del_bet,
                                              c, cloud, ice);
}
