#include "InputAndToneStage.h"

// Bring node names into scope for readability
using namespace V1A_Nodes;

InputAndToneStage::InputAndToneStage(double sampleRate)
    : MNASolver<V1A_TONESTACK_NODES>(sampleRate, true), // true for non-linear
      p_treble(0.5), p_mid(0.5), p_bass(0.5), p_vol1(0.75)
{
    flagForRecalculation(); // Initial calculation
}

void InputAndToneStage::setParams(double treble, double mid, double bass, double vol1) {
    p_treble = treble;
    p_mid = mid;
    p_bass = bass;
    p_vol1 = vol1;
    flagForRecalculation(); // Tell the solver matrices need to be rebuilt
}

void InputAndToneStage::stampLinear(Matrix<V1A_TONESTACK_NODES>& G, Matrix<V1A_TONESTACK_NODES>& C) {
    // This is the same stamping code from the previous answer's updateSystem() method.
    // It only deals with R and C components.
    
    // --- Component Values ---
    const double R4 = 150e3, R2 = 1.5e3, C1_C2 = 22.47e-6;
    const double R5 = 100e3, C_in = 1e-9, C4 = 0.1e-6, C3 = 0.047e-6;
    const double R5A = 100e3, C13B = 180e-12;
    // Potentiometers
    double treble_sq = p_treble * p_treble;
    const double RA_TREBLE = 250e3 * (1.0 - treble_sq), RC_TREBLE = 250e3 * treble_sq;
    double bass_sq = p_bass * p_bass;
    const double RA_BASS = 250e3 * bass_sq;
    double mid_sq = p_mid * p_mid;
    const double RA_MID = 10e3 * mid_sq;
    double vol1_sq = p_vol1 * p_vol1;
    const double RA_VOLUME1 = 1e6 * (1.0 - vol1_sq), RC_VOLUME1 = 1e6 * vol1_sq;

    auto g = [](double r) { return (r > 1e-9) ? 1.0 / r : 1e9; };

    // V1A Stage Components
    G[N004][N004] += g(R4);
    G[N033][N033] += g(R2);
    C[N033][N033] += C1_C2;
    
    // Tone Stack Components (full stamping logic here...)
    G[N004][N004] += g(R5); G[N015][N015] += g(R5); G[N004][N015] -= g(R5); G[N015][N004] -= g(R5);
    C[N004][N004] += C_in; C[N005][N005] += C_in; C[N004][N005] -= C_in; C[N005][N004] -= C_in;
    // ... and so on for all other linear components
}

void InputAndToneStage::stampNonlinear(Matrix<V1A_TONESTACK_NODES>& J, Vector<V1A_TONESTACK_NODES>& b, const Vector<V1A_TONESTACK_NODES>& x_k, const std::vector<double>& u_n) {
    // This is the same Newton-Raphson code from the previous answer's process() method.
    // It takes the Jacobian J and RHS vector b (already pre-filled with linear parts)
    // and adds the non-linear contributions from the tube.

    double v_p_prev = x_k[N004];
    double v_c_prev = x_k[N033];
    double v_g_prev = u_n[0]; // Grid voltage is our single input

    // Evaluate tube equations and their partial derivatives (same as before)
    double v_gp = v_g_prev - v_p_prev;
    double v_gc = v_g_prev - v_c_prev;
    double e1_val_arg = KP * (1.0/MU + v_gc / sqrt(KVB + v_gp * v_gp));
    // ... (full evaluation of g1_val and derivatives dG1_dVgc, dG1_dVgp)
    
    // For brevity, assuming g1_val, dG1_dVgc, dG1_dVgp are calculated here.
    double g1_val = 0.0; // placeholder
    double dG1_dVgc = 0.0; // placeholder
    double dG1_dVgp = 0.0; // placeholder
    
    // Add G1 (plate current) contribution to the RHS vector b
    b[N004] -= g1_val;
    b[N033] += g1_val;

    // Add G1's derivatives to the Jacobian J
    J[N004][N004] += -dG1_dVgp;
    J[N004][N033] += -dG1_dVgc;
    J[N033][N004] += dG1_dVgp;
    J[N033][N033] += dG1_dVgc;

    // The input grid voltage affects the derivatives, modifying the RHS
    double dG1_dVg = dG1_dVgp + dG1_dVgc;
    b[N004] -= dG1_dVg * v_g_prev; // u_n[0] is v_g_prev
    b[N033] += dG1_dVg * v_g_prev;
}
