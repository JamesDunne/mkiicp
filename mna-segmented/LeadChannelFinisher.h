#pragma once
#include "MNASolver.h"

// Node mapping for clarity
namespace LCF_Nodes {
    constexpr int N018 = 0, N030 = 1, PlateV4A = 2, CathodeV4A = 3, N026 = 4;
}

class LeadChannelFinisher : public MNASolver<5> {
public:
    LeadChannelFinisher(double sampleRate);
    double process(double v_in_N017) { return MNASolver::process(v_in_N017); }

protected:
    void stampLinear(Matrix<5>& G, Matrix<5>& C) override;
    void stampInputs(std::vector<std::pair<Vector<5>, Vector<5>>>& i) override;
    void stampNonlinear(Matrix<5>& J, Vector<5>& b, const Vector<5>& x, const std::vector<double>& u) override;
    int getOutputNodeIndex() const override { return 4; } // N026
};
