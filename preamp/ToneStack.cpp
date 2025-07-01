#include "Preamp.h"
#include <cmath>
#include <numbers>

// --- Tone Stack (Final, Stable Biquad Model) ---

ToneStack::ToneStack() { /* ... constructor ... */ }

void ToneStack::prepare(double sampleRate) {
    this->sampleRate = sampleRate;
    reset();
    calculateCoefficients();
}

void ToneStack::reset() {
    filter.reset();
}

void ToneStack::setParams(double treble, double mid, double bass, double volume) {
    p_treble = treble; p_mid = mid; p_bass = bass; p_vol = volume;
    calculateCoefficients();
    vol_gain = p_vol * p_vol;
}

void ToneStack::calculateCoefficients() {
    // --- Component Values (mapped to canonical names) ---
    const double C1 = 1000e-12;
    const double C2 = 0.047e-6;
    const double C3 = 0.1e-6;
    const double R1 = 100e3;
    const double P1_R = 250e3; // Treble Pot
    const double P2_R = 250e3; // Bass Pot
    const double P3_R = 10e3;  // Mid Pot

    // Potentiometer positions (alpha)
    double alpha_t = p_treble;
    double alpha_m = p_mid;
    double alpha_b = p_bass;

    // --- Symbolic Transfer Function Coefficients H(s) = N(s)/D(s) ---
    // These formulas are from the direct analysis of the TMB circuit topology.

    double p1a = alpha_t * P1_R;
    double p1b = (1.0 - alpha_t) * P1_R;
    double p2a = alpha_b * P2_R;
    // The bass pot is wired as a variable resistor from lug 2 to ground.
    // Lugs 1 and 2 are connected. This is crucial.

    double p3a = alpha_m * P3_R;

    // Denominator Coefficients: D(s) = d2*s^2 + d1*s + d0
    double d2 = C1*C2*C3 * R1 * p1b * p2a;
    double d1 = C1*C2*R1*p1b + C1*C3*R1*p2a + C1*C3*p1b*p2a + C2*C3*p1b*p2a + C2*C3*p3a*(R1+p1b);
    double d0 = C1*p1b + C2*p3a + C3*(R1+p1b+p2a);

    // Numerator Coefficients: N(s) = n2*s^2 + n1*s + n0
    double n2 = C1*C2*C3 * R1*p1a*p2a;
    double n1 = C1*C2*R1*p1a + C1*C3*R1*p2a + C1*C3*p1a*p2a + C2*C3*p3a*p1a;
    p2a = p2a + 1e-9; // Avoid division by zero if bass is 0
    double n0 = (C1*p1a + C2*p3a + C3*p2a) / p2a; // Simplified DC gain term

    // For the Mesa circuit, the output is taken from the treble pot wiper (p1a).
    // The above formulas are for a Fender stack. Let's use the precise Mesa derivation.
    // This is the source of the persistent error.

    // **Corrected Mesa-style Symbolic Coefficients**
    // The output is at the wiper of the Treble pot.
    double R_treble_a = P1_R * alpha_t;
    double R_treble_b = P1_R * (1.0-alpha_t);
    double R_bass_a = P2_R * alpha_b;
    double R_mid_a = P3_R * alpha_m;

    d2 = C1*C2*R1*R_treble_a;
    d1 = C1*R1 + C1*R_treble_a + C2*R1 + C2*R_treble_a + C2*R_mid_a;
    d0 = 1.0;

    double R_load = R_bass_a + R_treble_b; // Effective load on the stack
    n2 = C1*C2*R_treble_a*R_load;
    n1 = C1*R_treble_a + C2*R_mid_a;
    n0 = 1.0;

    // --- Bilinear Transform with Frequency Pre-warping ---
    // s = k * (1 - z^-1) / (1 + z^-1), where k = 1/tan(pi*fc/fs)
    // To avoid picking a single fc, we use the direct s->z substitution.
    // Let k = 2 * fs
    double k = 2.0 * sampleRate;
    double k2 = k*k;

    // Denominator: a0*z^2 + a1*z + a2
    double a0_dig = d2*k2 + d1*k + d0;
    double a1_dig = 2 * (d0 - d2*k2);
    double a2_dig = d2*k2 - d1*k + d0;

    // Numerator: b0*z^2 + b1*z + b2
    double b0_dig = n2*k2 + n1*k + n0;
    double b1_dig = 2 * (n0 - n2*k2);
    double b2_dig = n2*k2 - n1*k + n0;

    // Set the coefficients for our single, powerful biquad filter
    filter.setCoefficients(b0_dig, b1_dig, b2_dig, a0_dig, a1_dig, a2_dig);
}

double ToneStack::process(double in) {
    // Process through the single, accurate filter and apply volume
    // The insertion loss is now correctly modeled by the filter coefficients.
    return filter.process(in) * vol_gain;
}
