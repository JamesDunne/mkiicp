#pragma once
#include "MNASolver.h"

class TubeStageV1B : public MNASolver<2> { // 2 nodes: Plate, Cathode
public:
    TubeStageV1B(double sampleRate);
    double process(double v_grid_N020) { return MNASolver::process(v_grid_N020); }

protected:
    void stampLinear(Matrix<2>& G, Matrix<2>& C) override;
    void stampInputs(std::vector<std::pair<Vector<2>, Vector<2>>>& i) override {} // N/A
    void stampNonlinear(Matrix<2>& J, Vector<2>& b, const Vector<2>& x, const std::vector<double>& u) override;
    int getOutputNodeIndex() const override { return 0; } // Node 0 = Plate (N009)
};
