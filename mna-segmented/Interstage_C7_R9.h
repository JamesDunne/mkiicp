#pragma once
#include "MNASolver.h"

class Interstage_C7_R9 : public MNASolver<1> {
public:
    Interstage_C7_R9(double sampleRate);
    double process(double v_in_N009) { return MNASolver::process(v_in_N009); }

protected:
    void stampLinear(Matrix<1>& G, Matrix<1>& C) override;
    void stampInputs(std::vector<std::pair<Vector<1>, Vector<1>>>& i) override;
    void stampNonlinear(Matrix<1>&, Vector<1>&, const Vector<1>&, const std::vector<double>&) override {}
    int getOutputNodeIndex() const override { return 0; }
};
