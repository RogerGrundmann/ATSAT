/*
 * Saturn Atmosphere Circulation Model (ATSAT)
 * Dynamic-pressure Poisson solve — ATSAT's binding of the SHARED implementation in
 * PressureSolver.h.
 *
 * This is NOT a bare typedef, because ATSAT has something ATJUP does not: an inherited solver,
 * cSaturnModel::computePressure() in Pressure_Sat.cpp, which was the default for this model's
 * whole history. ATSAT_PRESS_SOLVER selects between them — DEFAULT 1, the shared solver; 0 goes
 * back to the inherited one, which is kept precisely so that comparison stays one environment
 * variable away.
 *
 * Everything that made the two solvers differ now lives either in the shared file or in
 * cSaturnModel: has_obstacle(), press_rigid_lid() and prepareProjectionBoundaries() are the three
 * things the shared solver asks this planet, and the header of PressureSolver.h says why each one
 * is a model's business rather than the solver's.
 */

#pragma once

#include "PressureSolver.h"

class cSaturnModel;

class PressureSolverSat {
public:
    explicit PressureSolverSat(cSaturnModel& model) : m(model) {}

    // DEFAULT 1 = the shared solver. It became the default on a measured 43% lower Poisson
    // residual and a dynamic pressure that no longer goes negative; see the commit.
    static int shared_enabled(){
        static const int v = [](){
            const char* e = getenv("ATSAT_PRESS_SOLVER"); return e ? atoi(e) : 1; }();
        return v;
    }

    void run(){
        if(shared_enabled() == 0){ m.computePressure(); return; }
        PressureSolver<cSaturnModel>(m).run();
    }

private:
    cSaturnModel& m;
};
