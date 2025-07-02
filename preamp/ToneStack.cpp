
#include "ToneStack.h"

#ifndef USE_TONESTACK_MNA

#include <ranges>
#include "Preamp.h"

// --- Tone Stack (Final, Stable Biquad Model) ---

void make_hpf(IIRBiquad& filter, double sampleRate, double cutoff_freq);

ToneStack::ToneStack() { /* ... */ }

void ToneStack::prepare(double sampleRate, double v1a_Vp_dc) {
    this->sampleRate = sampleRate;
    calculateCoefficients();

    make_hpf(dcBlocker, sampleRate, 5.0); // 5Hz cutoff for fast settling
    dcBlocker.prime(v1a_Vp_dc); // Prime with the DC voltage it will receive
}

void ToneStack::reset() {
    filter.reset();
    dcBlocker.reset();
}

void ToneStack::setParams(double treble, double mid, double bass, double volume) {
    p_treble = treble; p_mid = mid; p_bass = bass; p_vol = volume;
    calculateCoefficients();
    vol_gain = p_vol * p_vol;
}

// Uses the correct, canonical transfer function for a TMB stack
void ToneStack::calculateCoefficients() {
    const double C1 = 1000e-12, C2 = 0.047e-6, R1 = 100e3;
    const double P1_R = 250e3, P2_R = 250e3, P3_R = 10e3;

    double R_t = P1_R * p_treble, R_b = P2_R * p_bass, R_m = P3_R * p_mid;

    double d2 = C1*C2*R1*R_t;
    double d1 = C1*R1 + C1*R_t + C2*R1 + C2*R_t + C2*R_m;
    double d0 = 1.0;

    double n2 = C1*C2*R_t*R_b;
    double n1 = C1*R_t + C2*R_m;
    double n0 = 1.0;

    double k = 2.0 * sampleRate;
    double k2 = k*k;

    double a0 = d2*k2 + d1*k + d0, a1 = 2*(d0 - d2*k2), a2 = d2*k2 - d1*k + d0;
    double b0 = n2*k2 + n1*k + n0, b1 = 2*(n0 - n2*k2), b2 = n2*k2 - n1*k + n0;

    filter.setCoefficients(b0, b1, b2, a0, a1, a2);
}

double ToneStack::process(double in) {
    double ac_signal = dcBlocker.process(in);
    double eq_signal = filter.process(ac_signal);
    return eq_signal * vol_gain;
}

#else

#include <iostream>

ToneStackFilter::ToneStackFilter() :
    sampleRate(48000.0), T(1.0 / sampleRate),
    p_treble(0.5), p_mid(0.5), p_bass(0.5), p_vol1(0.75),
    u_prev(0.0)
{
    // Initialize previous voltage vector to zeros
    for(int i = 0; i < NUM_NODES; ++i) v_prev[i] = 0.0;
    updateCoefficients();
}

void ToneStackFilter::setSampleRate(double newSampleRate) {
    if (newSampleRate > 0 && newSampleRate != sampleRate) {
        sampleRate = newSampleRate;
        T = 1.0 / sampleRate;
        updateCoefficients();
    }
}

void ToneStackFilter::setParams(double treble, double mid, double bass, double vol1) {
    p_treble = treble;
    p_mid = mid;
    p_bass = bass;
    p_vol1 = vol1;
    updateCoefficients();
}

void ToneStackFilter::updateCoefficients() {
    // --- Define Component Values ---
    const double R5 = 100e3;
    const double C_in = 1e-9; // C5+C6 = 250p + 750p
    const double C4 = 0.1e-6;
    const double C3 = 0.047e-6;
    const double R5A = 100e3;
    const double C13B = 180e-12;

    // --- Calculate Potentiometer Resistances ---
    // The x*x curve provides a more audio-like taper than a linear pot
    double treble_sq = p_treble * p_treble;
    const double RA_TREBLE = 250e3 * (1.0 - treble_sq);
    const double RC_TREBLE = 250e3 * treble_sq;

    double bass_sq = p_bass * p_bass;
    const double RA_BASS = 250e3 * bass_sq; // Bass pot is wired reverse

    double mid_sq = p_mid * p_mid;
    const double RA_MID = 10e3 * mid_sq;

    double vol1_sq = p_vol1 * p_vol1;
    const double RA_VOLUME1 = 1e6 * (1.0 - vol1_sq);
    const double RC_VOLUME1 = 1e6 * vol1_sq;

    // --- Build G (conductance) and C (capacitance) matrices ---
    // These represent the system: (G + sC)V = F
    Matrix<NUM_NODES> G, C;
    Vector<NUM_NODES> Fg, Fc; // Input vectors F = (Fg + s*Fc)*u_in

    // Conductance mapping for readability
    auto g = [](double r) { return (r > 1e-9) ? 1.0 / r : 1e9; };

    // Node 1: N020 (Output)
    G[0][0] = g(RA_VOLUME1) + g(RC_VOLUME1);
    G[0][4] = -g(RA_VOLUME1);
    C[0][0] = C13B;
    C[0][4] = -C13B;

    // Node 2: N015
    G[1][1] = g(R5);
    C[1][1] = C4 + C3;
    C[1][5] = -C4;
    C[1][6] = -C3;
    Fg[1] = g(R5); // Connection to input node N004 via R5

    // Node 3: N005
    G[2][2] = g(RA_TREBLE);
    G[2][3] = -g(RA_TREBLE);
    C[2][2] = C_in;
    Fc[2] = C_in; // Connection to input node N004 via C_in

    // Node 4: N007
    G[3][3] = g(R5A) + g(RA_TREBLE) + g(RC_TREBLE);
    G[3][2] = -g(RA_TREBLE);
    G[3][4] = -g(R5A);
    G[3][5] = -g(RC_TREBLE);

    // Node 5: N008
    G[4][4] = g(R5A) + g(RA_VOLUME1);
    G[4][0] = -g(RA_VOLUME1);
    G[4][3] = -g(R5A);
    C[4][4] = C13B;
    C[4][0] = -C13B;

    // Node 6: N016
    G[5][5] = g(RC_TREBLE) + g(RA_BASS);
    G[5][1] = 0; // Filled from other side
    G[5][3] = -g(RC_TREBLE);
    G[5][6] = -g(RA_BASS);
    C[5][5] = C4;
    C[5][1] = -C4;

    // Node 7: N028
    G[6][6] = g(RA_BASS) + g(RA_MID);
    G[6][1] = 0;
    G[6][5] = -g(RA_BASS);
    C[6][6] = C3;
    C[6][1] = -C3;

    // Symmetrize matrices (MNA produces symmetric G and C for R/C circuits)
    for (int i = 0; i < NUM_NODES; ++i) {
        for (int j = i + 1; j < NUM_NODES; ++j) {
            G[j][i] = G[i][j];
            C[j][i] = C[i][j];
        }
    }

    // --- Calculate discrete-time system matrices ---
    // Using the Bilinear Transform method
    double g_s = 2.0 / T;
    Matrix<NUM_NODES> A_mat = G + g_s * C;
    B_mat = G - g_s * C;
    C_vec = Fg + g_s * Fc;
    D_vec = Fg - g_s * Fc;

    // The most expensive step: inverting the A matrix.
    // This is only done when parameters change, not per-sample.
    try {
        A_inv = Matrix<NUM_NODES>::invert(A_mat);
    } catch (const std::runtime_error& e) {
        // Fallback or error logging
        std::cerr << "Error during matrix inversion: " << e.what() << std::endl;
        // You might want to load a safe "bypass" identity matrix here
    }
}

double ToneStackFilter::process(double inputSample) {
    // Current input
    double u_current = inputSample;

    // The core update equation:
    // A_inv * (C_vec*u[n] + D_vec*u[n-1] - B_mat*v[n-1])
    Vector<NUM_NODES> term1 = u_current * C_vec;
    Vector<NUM_NODES> term2 = u_prev * D_vec;
    Vector<NUM_NODES> term3 = B_mat * v_prev;

    Vector<NUM_NODES> rhs = term1 + term2 - term3;

    Vector<NUM_NODES> v_current = A_inv * rhs;

    // The output is the voltage at node N020 (our first unknown node)
    double outputSample = v_current[0];

    // Update state for the next iteration
    v_prev = v_current;
    u_prev = u_current;

    return outputSample;
}

#endif
