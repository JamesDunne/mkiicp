#pragma once
#include "MNASolver.h"

class LeadDriveNetwork : public MNASolver<3> {
public:
    LeadDriveNetwork(double sampleRate);
    void setGain(double gain); // 0.0 to 1.0
    double process(double v_in_N001) { return MNASolver::process(v_in_N001); }

protected:
    void stampLinear(Matrix<3>& G, Matrix<3>& C) override;
    void stampInputs(std::vector<std::pair<Vector<3>, Vector<3>>>& i) override;
    void stampNonlinear(Matrix<3>&, Vector<3>&, const Vector<3>&, const std::vector<double>&) override {}
    int getOutputNodeIndex() const override { return 2; } // N029

private:
    double p_gain;
};
