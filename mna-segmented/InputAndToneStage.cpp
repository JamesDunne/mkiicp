#include "InputAndToneStage.h"
#include "TubeModel.h" // Includes our helper for the non-linear tube calculations

// Use the namespace for node indices for improved readability in the implementation.
using namespace V1A_Nodes;

InputAndToneStage::InputAndToneStage(double sampleRate)
    : MNASolver<V1A_TONESTACK_NODES>(sampleRate, true), // true = this is a non-linear circuit
      p_treble(0.5), p_mid(0.5), p_bass(0.5), p_vol1(0.75) {
    // Flag for the initial calculation of the system matrices.
    flagForRecalculation();
}

void InputAndToneStage::setParams(double treble, double mid, double bass, double vol1) {
    p_treble = treble;
    p_mid = mid;
    p_bass = bass;
    p_vol1 = vol1;
    // Notify the base solver that parameters have changed and a recalculation is needed.
    flagForRecalculation();
}

void InputAndToneStage::stampLinear(Matrix<V1A_TONESTACK_NODES>& G, Matrix<V1A_TONESTACK_NODES>& C) {
    // Helper lambda for converting resistance to conductance.
    auto g = [](double r) { return (r > 1e-12) ? 1.0 / r : 1e12; };

    // --- Define Component Values ---
    // V1A Stage Components
    const double R4 = 150e3;        // Plate load V1A
    const double R2 = 1.5e3;        // Cathode resistor V1A
    const double C1_C2 = 22.47e-6;  // C1+C2 in parallel

    // Tone Stack Fixed Components
    const double R5 = 100e3;
    const double C_in = 1e-9;       // C5+C6
    const double C4 = 0.1e-6;
    const double C3 = 0.047e-6;
    const double R5A = 100e3;
    const double C13B = 180e-12;    // Bright cap

    // --- Calculate Potentiometer Resistances ---
    double treble_sq = p_treble * p_treble;
    const double RA_TREBLE = 250e3 * (1.0 - treble_sq);
    const double RC_TREBLE = 250e3 * treble_sq;

    double bass_sq = p_bass * p_bass;
    const double RA_BASS = 250e3 * bass_sq;

    double mid_sq = p_mid * p_mid;
    const double RA_MID = 10e3 * mid_sq;

    double vol1_sq = p_vol1 * p_vol1;
    const double RA_VOLUME1 = 1e6 * (1.0 - vol1_sq);
    const double RC_VOLUME1 = 1e6 * vol1_sq;

    // --- Stamp Linear Components into G (conductance) and C (capacitance) matrices ---

    // V1A Stage
    G[N004][N004] += g(R4); // Plate load to B+ (AC ground)
    G[N033][N033] += g(R2); // Cathode resistor to ground
    C[N033][N033] += C1_C2; // Cathode bypass caps to ground

    // Tone Stack Components
    // R5 (N004 <-> N015)
    G[N004][N004] += g(R5); G[N015][N015] += g(R5);
    G[N004][N015] -= g(R5); G[N015][N004] -= g(R5);

    // C_in (N004 <-> N005)
    C[N004][N004] += C_in; C[N005][N005] += C_in;
    C[N004][N005] -= C_in; C[N005][N004] -= C_in;

    // C4 (N015 <-> N016)
    C[N015][N015] += C4; C[N016][N016] += C4;
    C[N015][N016] -= C4; C[N016][N015] -= C4;

    // C3 (N015 <-> N028)
    C[N015][N015] += C3; C[N028][N028] += C3;
    C[N015][N028] -= C3; C[N028][N015] -= C3;

    // Treble Pot (RA_TREBLE: N005<->N007, RC_TREBLE: N007<->N016)
    G[N005][N005] += g(RA_TREBLE); G[N007][N007] += g(RA_TREBLE);
    G[N005][N007] -= g(RA_TREBLE); G[N007][N005] -= g(RA_TREBLE);
    G[N007][N007] += g(RC_TREBLE); G[N016][N016] += g(RC_TREBLE);
    G[N007][N016] -= g(RC_TREBLE); G[N016][N007] -= g(RC_TREBLE);

    // Bass Pot (RA_BASS: N016 <-> N028)
    G[N016][N016] += g(RA_BASS); G[N028][N028] += g(RA_BASS);
    G[N016][N028] -= g(RA_BASS); G[N028][N016] -= g(RA_BASS);

    // Mid Pot (RA_MID: N028 -> Gnd)
    G[N028][N028] += g(RA_MID);

    // R5A (N008 <-> N007)
    G[N008][N008] += g(R5A); G[N007][N007] += g(R5A);
    G[N008][N007] -= g(R5A); G[N007][N008] -= g(R5A);

    // Volume Pot (RA_VOLUME1: N020<->N008, RC_VOLUME1: N020->Gnd)
    G[N020][N020] += g(RC_VOLUME1);
    G[N020][N020] += g(RA_VOLUME1); G[N008][N008] += g(RA_VOLUME1);
    G[N020][N008] -= g(RA_VOLUME1); G[N008][N020] -= g(RA_VOLUME1);

    // Bright Cap C13B (N020 <-> N008)
    C[N020][N020] += C13B; C[N008][N008] += C13B;
    C[N020][N008] -= C13B; C[N008][N020] -= C13B;
}

void InputAndToneStage::stampNonlinear(Matrix<V1A_TONESTACK_NODES>& J, Vector<V1A_TONESTACK_NODES>& b, const Vector<V1A_TONESTACK_NODES>& x_k, const std::vector<double>& u_n) {
    // This function adds the contribution of the non-linear V1A tube.
    // The grid voltage for V1A is the external input to this entire block.
    double v_grid = u_n[0];

    // Delegate the complex calculations to the TubeModel helper function.
    // We provide the Jacobian (J), RHS vector (b), previous state (x_k),
    // current grid voltage (v_grid), and the node indices for the tube's plate and cathode.
    TubeModel::stamp(J, b, x_k, v_grid, N004, N033);
}
