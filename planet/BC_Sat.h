#pragma once

class cSaturnModel;

class BC_Sat {
public:
    explicit BC_Sat(cSaturnModel& model) : m(model) {}

    void bcRadius();
    void bcTheta();
    void bcPhi();

private:
    cSaturnModel& m;
};

// Implementations delegate to cSaturnModel::BC_radius/theta/phi() in BC_Sat.cpp.
