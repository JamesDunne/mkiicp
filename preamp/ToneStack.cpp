#include "Preamp.h"
#include <cmath>
#include <numbers>

// --- Tone Stack (Final, Stable Biquad Model) ---

void make_hpf(IIRBiquad& filter, double sampleRate, double cutoff_freq);

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
    make_hpf(dcBlocker, sampleRate, 1.0);
}

/**
 * @brief Calculates the digital filter coefficients from the knob positions.
 *
 * This method uses a direct symbolic transfer function derived from circuit analysis
 * of the Mesa-style TMB (Treble-Mid-Bass) tone stack. The output is taken from the
 * treble pot wiper, which correctly models this specific topology. The formulas
 * capture the complex interactions between all the controls.
 */
void ToneStack::calculateCoefficients() {
    // --- 1. Define physical component values from the SPICE schematic ---
    const double C_Treble = 1000e-12; // Combination of C5 and C6
    const double C_Mid = 0.047e-6;    // C3 in schematic
    const double R_In = 100e3;        // R5 in schematic
    const double R_Pot_Treble = 250e3;
    const double R_Pot_Bass = 250e3;
    const double R_Pot_Mid = 10e3;

    // --- 2. Calculate variable resistances from knob positions (alpha: 0.0 to 1.0) ---
    double trebleWiperRes = R_Pot_Treble * p_treble;
    double trebleRemRes = R_Pot_Treble * (1.0 - p_treble);
    double bassLoadRes = R_Pot_Bass * p_bass;
    double midScoopRes = R_Pot_Mid * p_mid;

    // --- 3. Calculate Analog Transfer Function Coefficients ---
    // H(s) = N(s) / D(s) = (n2*s^2 + n1*s + n0) / (d2*s^2 + d1*s + d0)

    // Denominator coefficients (define the filter's poles/resonances)
    double d2 = C_Treble * C_Mid * R_In * trebleWiperRes;
    double d1 = C_Treble * (R_In + trebleWiperRes) + C_Mid * (R_In + trebleWiperRes + midScoopRes);
    double d0 = 1.0;

    // Numerator coefficients (define the filter's zeros/notches)
    // The "effective load" on the stack is a key part of the sound.
    double effectiveLoad = bassLoadRes + trebleRemRes;
    double n2 = C_Treble * C_Mid * trebleWiperRes * effectiveLoad;
    double n1 = C_Treble * trebleWiperRes + C_Mid * midScoopRes;
    double n0 = 1.0;

    // --- 4. Convert Analog to Digital via Bilinear Transform ---
    // This maps the continuous-time H(s) to the discrete-time H(z).
    // The formula used is s = 2*fs * (1-z^-1)/(1+z^-1)
    double k = 2.0 * sampleRate;
    double k_squared = k * k;

    // Digital denominator coefficients (a_coeffs)
    double a0_digital = d2 * k_squared + d1 * k + d0;
    double a1_digital = 2.0 * (d0 - d2 * k_squared);
    double a2_digital = d2 * k_squared - d1 * k + d0;

    // Digital numerator coefficients (b_coeffs)
    double b0_digital = n2 * k_squared + n1 * k + n0;
    double b1_digital = 2.0 * (n0 - n2 * k_squared);
    double b2_digital = n2 * k_squared - n1 * k + n0;

    // --- 5. Set coefficients on the IIR filter ---
    filter.setCoefficients(b0_digital, b1_digital, b2_digital,
                           a0_digital, a1_digital, a2_digital);
}

double ToneStack::process(double in) {
    // First, process the signal through the EQ filter.
    double out = filter.process(in);

    // Then, apply the volume control.
    out *= vol_gain;

    // Finally, remove the DC offset.
    return dcBlocker.process(out);
}
