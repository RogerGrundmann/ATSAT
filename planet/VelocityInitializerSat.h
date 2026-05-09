#pragma once

class cSaturnModel;

class VelocityInitializerSat {
public:
    explicit VelocityInitializerSat(cSaturnModel& model) : m(model) {}

    void compute();

private:
    cSaturnModel& m;
};

// Implementation delegates to cSaturnModel::SaturnCellStructure() in InitVelocity_Sat.cpp.
