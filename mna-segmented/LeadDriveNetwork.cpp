#include "LeadDriveNetwork.h"

namespace LDN_Nodes { constexpr int N011 = 0, N024 = 1, N029 = 2; }

LeadDriveNetwork::LeadDriveNetwork(double sr) : MNASolver<3>(sr, false), p_gain(0.5) {
    flagForRecalculation();
}

void LeadDriveNetwork::setGain(double gain) {
    p_gain = gain;
    flagForRecalculation();
}

void LeadDriveNetwork::stampLinear(Matrix<3>& G, Matrix<3>& C) {
    using namespace LDN_Nodes;
    const double C21 = 0.02e-6, R21 = 680e3, R22 = 475e3;
    double gain_sq = p_gain * p_gain;
    const double RA_LD = 1e6 * (1.0 - gain_sq), RC_LD = 1e6 * gain_sq;

    // C21
    C[N011][N011] += C21;
    // R21
    G[N011][N011] += 1.0/R21; G[N024][N024] += 1.0/R21;
    G[N011][N024] -= 1.0/R21; G[N024][N011] -= 1.0/R21;
    // Lead Drive Pot
    G[N024][N024] += 1.0/RA_LD; G[N029][N029] += 1.0/RA_LD;
    G[N024][N029] -= 1.0/RA_LD; G[N029][N024] -= 1.0/RA_LD;
    G[N029][N029] += 1.0/RC_LD;
    // R22
    G[N029][N029] += 1.0/R22;
}

void LeadDriveNetwork::stampInputs(std::vector<std::pair<Vector<3>, Vector<3>>>& i) {
    using namespace LDN_Nodes;
    const double C21 = 0.02e-6;
    double g_s = 2.0 * 48000.0;
    Vector<3> C_vec, D_vec;
    C_vec[N011] = g_s * C21;
    D_vec[N011] = -g_s * C21;
    i.push_back({C_vec, D_vec});
}
