#pragma once

class cSaturnModel;

class PressureSolverSat {
public:
    explicit PressureSolverSat(cSaturnModel& model) : m(model) {}

    void run();

private:
    cSaturnModel& m;
};

// Implementation delegates to cSaturnModel::computePressure() in Pressure_Sat.cpp.
