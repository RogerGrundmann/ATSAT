/*
 * Saturn Atmosphere Circulation Model (ATSAT)
 * Dry convective adjustment — ATSAT's binding of the SHARED implementation.
 *
 * The algorithm, the reasoning and the knobs live in ConvectiveAdjustment.h, which is
 * byte-identical to ATJUP's copy and knows nothing about which planet it runs on. This file
 * exists only so the call sites keep their familiar name.
 */

#pragma once

#include "ConvectiveAdjustment.h"

class cSaturnModel;

typedef ConvectiveAdjustment<cSaturnModel> ConvectiveAdjustmentSat;
