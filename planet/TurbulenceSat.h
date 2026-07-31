/*
 * Saturn Atmosphere Circulation Model (ATSAT)
 * Turbulence closure — ATSAT's binding of the SHARED implementation in Turbulence.h.
 */

#pragma once

#include "Turbulence.h"

class cSaturnModel;

typedef Turbulence<cSaturnModel> TurbulenceSat;
