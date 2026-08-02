/*
 * SHARED PHYSICS — small helpers every planet model's physics needs. MUST BE BYTE-IDENTICAL IN
 * EVERY MODEL THAT USES IT; `make check-shared` verifies that against planet/SHARED.md5.
 *
 * Nothing here knows which planet it is running on. The environment helpers take the model's own
 * planet_tag() so one shared implementation reads each model's knobs under its own prefix, and
 * the saturation vapour pressure is the single expression that ATJUP, ATSAT and their saturation
 * adjustments were each computing separately — character for character the same formula in three
 * places, which is three places for it to drift.
 */

#pragma once

#include <cmath>
#include <cstdlib>
#include <string>

namespace ATPhys {


// Look up "<PLANET>_<name>", e.g. ATJUP_CONV_ADJ_LAPSE, so one shared implementation reads each
// model's own knobs under its own prefix. Not cached: the caller caches where it matters.
inline const char* env_for(const char* planet_tag, const char* name){
    std::string key(planet_tag);
    key += "_";
    key += name;
    return getenv(key.c_str());
}

inline double env_double(const char* planet_tag, const char* name, double fallback){
    const char* e = env_for(planet_tag, name);
    return e ? atof(e) : fallback;
}

inline int env_int(const char* planet_tag, const char* name, int fallback){
    const char* e = env_for(planet_tag, name);
    return e ? atoi(e) : fallback;
}


// Saturation vapour pressure [bar] from the Sanchez-Lavega / Clausius-Clapeyron coefficient set.
//
// This is the expression SaturationAdjustmentJup::saturation_vapour_pressure,
// SaturationAdjustmentSat::saturation_vapour_pressure and cSaturnModel::saturation_vapour_pressure
// all evaluated, in that order of discovery and with identical operand order. All four copies are
// gone now and every caller in both models comes here; R arrives in J/(kg K) and the 1e-3 takes it
// to the J/(g K) the coefficients are tabulated against.
inline double saturation_vapour_pressure(double T_K, double C, double L0, double R,
                                         double del_alf, double del_bet){
    return std::exp(C + (-L0 / T_K + del_alf * std::log(T_K) + del_bet * T_K) / (1e-3 * R));
}


// Two-coefficient Clausius-Clapeyron, E = exp(A/T + B). The short form used where a full
// Sanchez-Lavega coefficient set is not tabulated; Thermo_Jup.cpp is its only caller.
//
// It arrives here from SaturationAdjustmentJup and SaturationAdjustmentSat, which held it as
// character-for-character identical statics. ATSAT's copy had NO caller at all — it existed only
// because the class was written to mirror ATJUP's — so this is one live implementation replacing
// one live and one ornamental copy.
inline double clausius_clapeyron(double T_K, double A, double B){
    return std::exp(A / T_K + B);
}

}  // namespace ATPhys
