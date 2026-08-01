/*
 * Saturn Atmosphere Circulation Model (ATSAT)
 * Grey multi-layer radiation — ATSAT's binding of the SHARED implementation in Radiation.h.
 * Saturn's own radiative constants live in cSaturnModel.h, next to the rest of its parameters.
 */

#pragma once

#include "Radiation.h"

class cSaturnModel;

typedef Radiation<cSaturnModel> RadiationSat;
