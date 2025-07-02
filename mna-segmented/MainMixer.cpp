#include "MainMixer.h"

MainMixer::MainMixer(double sr) : MNASolver<1>(sr, false) {
    flagForRecalculation();
}

void MainMixer::stampLinear(Matrix<1>& G, Matrix<1>& C) {
    // All components loading N002
    const double R10=3.3e6, C10=20e-12, R11=680e3;
    const double R31=220e3, C31=250e-12;
    const double C11=47e-12, C32=500e-12, R32=100e3;

    G[0][0] = 1.0/R10 + 1.0/R11 + 1.0/R31 + 1.0/R32;
    C[0][0] = C10 + C31 + C11 + C32;
}

void MainMixer::stampInputs(std::vector<std::pair<Vector<1>, Vector<1>>>& i) {
    const double R10=3.3e6, C10=20e-12;
    const double R31=220e3, C31=250e-12;
    double g_s = 2.0 * 48000.0;

    // Input 1 from N001
    Vector<1> C1, D1;
    C1[0] = 1.0/R10 + g_s*C10;
    D1[0] = 1.0/R10 - g_s*C10;
    i.push_back({C1, D1});

    // Input 2 from N026
    Vector<1> C2, D2;
    C2[0] = 1.0/R31 + g_s*C31;
    D2[0] = 1.0/R31 - g_s*C31;
    i.push_back({C2, D2});
}
