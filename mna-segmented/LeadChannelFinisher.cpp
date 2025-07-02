#include "LeadChannelFinisher.h"
#include "TubeModel.h"

namespace LCF_Nodes {
    constexpr int N018 = 0, N030 = 1, PlateV4A = 2, CathodeV4A = 3, N026 = 4;
}

LeadChannelFinisher::LeadChannelFinisher(double sr) : MNASolver<5>(sr, true) {
    flagForRecalculation();
}

void LeadChannelFinisher::stampLinear(Matrix<5>& G, Matrix<5>& C) {
    using namespace LCF_Nodes;
    auto g = [](double r){ return 1.0/r; };
    // Interstage network
    const double C25 = 0.022e-6, R25 = 270e3, R24 = 68e3, C24 = 1000e-12;
    // V4A stage
    const double R27 = 274e3, R30 = 3.3e3, C29 = 0.22e-6;
    // Output coupling
    const double C30 = 0.047e-6, R31 = 220e3, C31 = 250e-12; // R31/C31 are load

    // C25 (Input -> N018) handled in stampInputs
    // R25/R24/C24 network
    G[N018][N018] += g(R25); G[N030][N030] += g(R25); G[N018][N030] -= g(R25); G[N030][N018] -= g(R25);
    G[N030][N030] += g(R24);
    C[N030][N030] += C24;
    // V4A Stage
    G[PlateV4A][PlateV4A] += g(R27);
    G[CathodeV4A][CathodeV4A] += g(R30);
    C[CathodeV4A][CathodeV4A] += C29;
    // Output Coupling C30
    C[PlateV4A][PlateV4A] += C30; C[N026][N026] += C30;
    C[PlateV4A][N026] -= C30; C[N026][PlateV4A] -= C30;
    // Output Load (from MainMixer)
    G[N026][N026] += g(R31);
    C[N026][N026] += C31;
}

void LeadChannelFinisher::stampInputs(std::vector<std::pair<Vector<5>, Vector<5>>>& i) {
    using namespace LCF_Nodes;
    const double C25 = 0.022e-6;
    double g_s = 2.0 * 48000.0;
    Vector<5> C_vec, D_vec;
    C_vec[N018] = g_s * C25;
    D_vec[N018] = -g_s * C25;
    i.push_back({C_vec, D_vec});
}

void LeadChannelFinisher::stampNonlinear(Matrix<5>& J, Vector<5>& b, const Vector<5>& x, const std::vector<double>& u) {
    using namespace LCF_Nodes;
    // The input to the tube grid is the voltage at node N030, which is part of our system x.
    // The overall segment input from N017 is handled by the linear solver part.
    double v_g_k = x[N030];
    TubeModel::stamp(J, b, x, v_g_k, PlateV4A, CathodeV4A);
}
