#include "LeadChannelFinisher.h"
#include "TubeModel.h"

using namespace LCF_Nodes;

LeadChannelFinisher::LeadChannelFinisher(double sr)
    : MNASolver<5>(sr, true) { // true for non-linear
    flagForRecalculation();
}

void LeadChannelFinisher::stampLinear(Matrix<5>& G, Matrix<5>& C) {
    auto g = [](double r){ return (r > 1e-12) ? 1.0/r : 1e12; };

    // Component Values
    const double C25 = 0.022e-6, R25 = 270e3, R24 = 68e3, C24 = 1000e-12;
    const double R27 = 274e3, R30 = 3.3e3, C29 = 0.22e-6;
    const double C30 = 0.047e-6;
    const double R_load = 220e3, C_load = 250e-12; // Load from MainMixer

    // Stamp C25's contribution to its own node's admittance.
    // This models the -C*d/dt*V(N018) part of the capacitor current.
    C[N018][N018] += C25;

    // Stamp the rest of the linear network
    G[N018][N018] += g(R25); G[N030][N030] += g(R25); G[N018][N030] -= g(R25); G[N030][N018] -= g(R25);
    G[N030][N030] += g(R24);
    C[N030][N030] += C24;
    G[PlateV4A][PlateV4A] += g(R27);
    G[CathodeV4A][CathodeV4A] += g(R30);
    C[CathodeV4A][CathodeV4A] += C29;
    C[PlateV4A][PlateV4A] += C30; C[N026][N026] += C30;
    C[PlateV4A][N026] -= C30; C[N026][PlateV4A] -= C30;
    G[N026][N026] += g(R_load);
    C[N026][N026] += C_load;
}

/**
 * This function is required by the MNASolver base class but is not used
 * for non-linear circuits because all source contributions are handled
 * dynamically within stampNonlinear. We provide an empty implementation
 * to satisfy the pure virtual requirement.
 */
void LeadChannelFinisher::stampInputs(std::vector<std::pair<Vector<5>, Vector<5>>>& i) {
    // This function is intentionally left empty.
}

void LeadChannelFinisher::stampNonlinear(Matrix<5>& J, Vector<5>& b, const Vector<5>& x_k, const std::vector<double>& u_n) {
    // --- 1. Model the driving current from the input capacitor C25 ---
    const double C25 = 0.022e-6;

    double v_in_current = u_n[0];
    double v_in_prev = (prev_inputs.empty()) ? 0.0 : prev_inputs[0];

    // Calculate the driving current I = C * dV/dt, approximated by the Bilinear Transform.
    // This models the +C*d/dt*V(input) part of the capacitor current.
    double i_in = C25 * g_s * (v_in_current - v_in_prev);

    // Add this current source to the RHS vector `b` at the row for node N018.
    b[N018] += i_in;

    // --- 2. Model the internal non-linear tube V4A ---
    // The grid voltage for tube V4A is the voltage at internal node N030.
    double v_g_k = x_k[N030];

    // Delegate the tube calculations to the helper.
    TubeModel::stamp(J, b, x_k, v_g_k, PlateV4A, CathodeV4A);
}
