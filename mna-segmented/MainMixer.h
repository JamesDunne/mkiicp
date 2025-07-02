#pragma once
#include "MNASolver.h"

class MainMixer : public MNASolver<1> { // Only 1 unknown node: N002
public:
    MainMixer(double sampleRate);
    double process(double v_in_N001, double v_in_N026) {
        return MNASolver::process({v_in_N001, v_in_N026});
    }

protected:
    void stampLinear(Matrix<1>& G, Matrix<1>& C) override;
    void stampInputs(std::vector<std::pair<Vector<1>, Vector<1>>>& i) override;
    void stampNonlinear(Matrix<1>&, Vector<1>&, const Vector<1>&, const std::vector<double>&) override {}
    int getOutputNodeIndex() const override { return 0; }
};
