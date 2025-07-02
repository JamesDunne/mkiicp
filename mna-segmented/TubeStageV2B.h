#pragma once
#include "MNASolver.h"

class TubeStageV2B : public MNASolver<2> {
public:
    TubeStageV2B(double sampleRate);
    double process(double v_grid_N002) { return MNASolver::process(v_grid_N002); }

protected:
    void stampLinear(Matrix<2>& G, Matrix<2>& C) override;
    void stampInputs(std::vector<std::pair<Vector<2>, Vector<2>>>& i) override {}
    void stampNonlinear(Matrix<2>& J, Vector<2>& b, const Vector<2>& x, const std::vector<double>& u) override;
    int getOutputNodeIndex() const override { return 0; } // Plate (N021)
};
