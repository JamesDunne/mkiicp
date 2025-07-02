#include "FinalStageAndMaster.h"
#include "TubeModel.h"

namespace FS_Nodes {
    constexpr int P001 = 0, N022 = 1, N032 = 2, GridV2A = 3, PlateV2A = 4, CathodeV2A = 5, N013 = 6, N014 = 7;
}

FinalStageAndMaster::FinalStageAndMaster(double sr) : MNASolver<8>(sr, true), p_master(0.5) {
    flagForRecalculation();
}

void FinalStageAndMaster::setMaster(double master) {
    p_master = master;
    flagForRecalculation();
}

void FinalStageAndMaster::stampLinear(Matrix<8>& G, Matrix<8>& C) {
    using namespace FS_Nodes;
    auto g = [](double r){ return 1.0/r; };
    // Input network
    const double C9 = 0.047e-6, R105=47e3, R46=47e3, R102=150e3, R101=4.7e3;
    // V2A Stage
    const double R19=120e3, R103=47e3, C16=.47e-6, C15=15e-6, R104=1e3, R12=2.2e3;
    // Master Volume
    const double C12=0.047e-6, R106=15e3;
    double master_sq = p_master * p_master;
    const double RA_MASTER = 1e6 * master_sq, RC_MASTER = 1e6 * (1.0 - master_sq); // Assuming reverse log pot

    // C9 (Input -> P001) handled in stampInputs
    // R105/R46
    G[P001][P001] += g(R105); G[N022][N022] += g(R105) + g(R46); G[P001][N022] -= g(R105); G[N022][P001] -= g(R105);
    // R102
    G[N022][N022] += g(R102); G[N032][N032] += g(R102); G[N022][N032] -= g(R102); G[N032][N022] -= g(R102);
    // R101
    G[N032][N032] += g(R101);
    // R12 (feedback)
    G[PlateV2A][PlateV2A] += g(R12); G[N032][N032] += g(R12); G[PlateV2A][N032] -= g(R12); G[N032][PlateV2A] -= g(R12);
    // R103 (grid leak for V2A)
    G[GridV2A][GridV2A] += g(R103);
    // V2A Plate/Cathode resistors
    G[PlateV2A][PlateV2A] += g(R19);
    G[CathodeV2A][CathodeV2A] += g(R104);
    C[CathodeV2A][CathodeV2A] += C16 + C15;
    // C12
    C[N013][N013] += C12; C[PlateV2A][PlateV2A] += C12; C[N013][PlateV2A] -= C12; C[PlateV2A][N013] -= C12;
    // R106 & Master Pot
    G[N013][N013] += g(R106); G[N014][N014] += g(R106); G[N013][N014] -= g(R106); G[N014][N013] -= g(R106);
    G[N014][N014] += g(RA_MASTER) + g(RC_MASTER);
}

void FinalStageAndMaster::stampInputs(std::vector<std::pair<Vector<8>, Vector<8>>>& i) {
    using namespace FS_Nodes;
    const double C9 = 0.047e-6;
    double g_s = 2.0 * 48000.0;
    Vector<8> C_vec, D_vec;
    C_vec[P001] = g_s * C9;
    D_vec[P001] = -g_s * C9;
    i.push_back({C_vec, D_vec});
}

void FinalStageAndMaster::stampNonlinear(Matrix<8>& J, Vector<8>& b, const Vector<8>& x, const std::vector<double>& u) {
    using namespace FS_Nodes;
    double v_g_k = x[GridV2A];
    TubeModel::stamp(J, b, x, v_g_k, PlateV2A, CathodeV2A);
}
