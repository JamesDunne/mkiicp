#include "Interstage_C7_R9.h"

Interstage_C7_R9::Interstage_C7_R9(double sampleRate) : MNASolver<1>(sampleRate, false) {
    flagForRecalculation();
}

void Interstage_C7_R9::stampLinear(Matrix<1>& G, Matrix<1>& C) {
    const double R9 = 100e3;
    const double C7 = 0.1e-6;
    G[0][0] += 1.0 / R9;
    C[0][0] += C7;
}

void Interstage_C7_R9::stampInputs(std::vector<std::pair<Vector<1>, Vector<1>>>& i) {
    const double C7 = 0.1e-6;
    double g_s = 2.0 * 48000.0; // Assume SR, better to get from base class
    Vector<1> C_vec, D_vec;
    C_vec[0] = g_s * C7;
    D_vec[0] = -g_s * C7;
    i.push_back({C_vec, D_vec});
}
