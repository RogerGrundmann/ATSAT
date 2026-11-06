/*
 * Saturn Atmosphere Circulation Model (ATSAT)
 * Turbulence closure — ATSAT's binding of the SHARED implementation in Turbulence.h.
 */

#pragma once

#include "Turbulence.h"

class cSaturnModel;

typedef Turbulence<cSaturnModel> TurbulenceSat;

// <TAG>_TURB_BUOY=3: Saturn's condensable gases for the moist N^2 in the shared closure -- the three the
// moist convective adjustment carries (H2O, NH3, CH4; cSaturnModel.cpp builds the same list for it), n = 0..2.
// Declared here so argument-dependent lookup finds it when Turbulence<cSaturnModel> is instantiated; the
// definition is in cSaturnModel.cpp (it reads private members through cSaturnModel::moist_gas_for_turb).
bool turb_moist_gas(cSaturnModel& m, int n, Array*& v, Array*& c, Array*& ice,
                    double& ep, double& lv, double& ls, double& t_0);
