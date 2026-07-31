/*
 * Saturn Atmosphere Circulation Model (ATSAT)
 * Ice / precipitation microphysics — ATSAT's binding of the SHARED implementation.
 */

#pragma once

#include "Precipitation.h"

class cSaturnModel;

typedef Precipitation<cSaturnModel> PrecipitationSat;
