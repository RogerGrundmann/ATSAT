/*
 * Saturn Atmosphere Circulation Model (ATSAT)
 * Dynamic-pressure Poisson solve — ATSAT's binding of the SHARED implementation in
 * PressureSolver.h. A bare typedef, as ATJUP's is.
 *
 * It was briefly a wrapper: ATSAT carried an inherited second solver,
 * cSaturnModel::computePressure(), and ATSAT_PRESS_SOLVER chose between them. Both are gone. The
 * shared solver had been the default since it was measured better, computePressure() had no
 * caller left but the gate, and a second implementation kept alive only to be switched off is the
 * duplication the shared headers exist to remove. What it knew that the shared solver did not —
 * this model's boundary conditions on the aux and rhs fields — was already moved into
 * cSaturnModel::prepareProjectionBoundaries() when the solver was shared, and stays there.
 */

#pragma once

#include "PressureSolver.h"

class cSaturnModel;

typedef PressureSolver<cSaturnModel> PressureSolverSat;
