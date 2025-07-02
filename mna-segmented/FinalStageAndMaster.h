#pragma once
#include "MNASolver.h"

class FinalStageAndMaster : public MNASolver<8> {
public:
    FinalStageAndMaster(double sampleRate);
    void setMaster(double master); // 0.0 to 1.0
    double process(double v_in_N021) { return MNASolver::process(v_in_N021); }

protected:
    void stampLinear(Matrix<8>& G, Matrix<8>& C) override;
    void stampInputs(std::vector<std::pair<Vector<8>, Vector<8>>>& i) override;
    void stampNonlinear(Matrix<8>& J, Vector<8>& b, const Vector<8>& x, const std::vector<double>& u) override;
    int getOutputNodeIndex() const override { return 7; } // N014

private:
    double p_master;
};
