#pragma once

#include <string>

class cSaturnModel;
class Array;

class SaturationAdjustmentSat {
public:
    explicit SaturationAdjustmentSat(cSaturnModel& model) : m(model) {}

    void run(std::string gas,
             double &coeff_A,   double &coeff_B,
             double &coeff_A_i, double &coeff_B_i,
             double &t_0,       double &t_00,
             double &ep,        double &lv,  double &ls,
             double &cp,        double &r,
             double &C,         double &L0,  double &R,
             double &del_alf,   double &del_bet,  double &m_mol,
             Array &c, Array &cloud, Array &ice);

private:
    cSaturnModel& m;
};

// Implementation is in Weather_Sat.cpp (as cSaturnModel::Saturation_Adjustment,
// delegated through the run() wrapper in cSaturnModel.cpp).
