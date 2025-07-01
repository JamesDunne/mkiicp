
#include "ToneStack.h"
#include "IIRFilter.h"
#include "Preamp.h"
#include <numbers>

// --- Tone Stack (Corrected Biquad Cascade/Parallel Model) ---

ToneStack::ToneStack() {
    sampleRate = 48000.0;
    p_treble = 0.8; p_mid = 0.5; p_bass = 0.25; p_vol = 0.75;
    reset();
}

void ToneStack::prepare(double sampleRate) {
    this->sampleRate = sampleRate;
    reset();
    calculateCoefficients();
}

void ToneStack::reset() {
    bassMidFilter.reset();
    trebleFilter.reset();
}

void ToneStack::setParams(double treble, double mid, double bass, double volume) {
    p_treble = treble; p_mid = mid; p_bass = bass; p_vol = volume;
    calculateCoefficients();
    volumeGain = p_vol * p_vol;
}

// In Preamp.cpp

// ... other includes and methods ...

void ToneStack::calculateCoefficients() {
    // --- Component Values ---
    const double R5 = 100e3;
    const double C_bright = 1000e-12;
    const double C4 = 0.1e-6;
    const double C3 = 0.047e-6;

    // --- Potentiometer values (0-1 range) ---
    double treble_alpha = p_treble * p_treble;
    double bass_alpha = p_bass * p_bass;
    double mid_alpha = 1.0 - (p_mid * p_mid);

    // --- Treble Filter (High-shelf / Peaking) ---
    // This models the path through the bright cap and treble pot.
    double R_treble_eq = 250e3 * (1.0 - treble_alpha) + 1.0; // Resistance opposing the bright cap
    double freq_treble = 1.0 / (2.0 * std::numbers::pi * R_treble_eq * C_bright);
    double K_t = tan(std::numbers::pi * freq_treble / sampleRate);

    // The gain V_t is constant, so the filter is always a peaking type.
    const double V_t = pow(10.0, 12.0 / 20.0); // Constant boost of ~12dB

    // Direct calculation for the peaking filter coefficients
    // The Quality Factor (Q) is implicitly set to 1.0 here for a broad peak.
    double norm_t = 1.0 / (1.0 + K_t + K_t * K_t);
    double b0_t = (1.0 + V_t * K_t + K_t * K_t) * norm_t;
    double b1_t = 2.0 * (K_t * K_t - 1.0) * norm_t;
    double b2_t = (1.0 - V_t * K_t + K_t * K_t) * norm_t;
    double a1_t = b1_t; // For this peaking filter topology, a1 = b1
    double a2_t = (1.0 - K_t + K_t * K_t) * norm_t;

    trebleFilter.setCoefficients(b0_t, b1_t, b2_t, 1.0, a1_t, a2_t);

    // --- Bass/Mid Filter (Band-pass/scoop) ---
    // This models the main path through C4 and the bass/mid controls
    double R_bass_eq = 250e3 * (1.0 - bass_alpha) + R5;
    double R_mid_eq = 10e3 * (mid_alpha) + 1.0;
    double freq_bass = 1.0 / (2.0 * std::numbers::pi * R_bass_eq * C4);

    double freq_mid = 1.0 / (2.0 * std::numbers::pi * R_mid_eq * C3);

    // Here, freq_mid is being used to define the Q of the bass band-pass filter,
    // which is a clever way to model their interaction.
    double Q_bm = freq_mid > 0 ? freq_bass / freq_mid : 0.1;

    double K_bm = tan(std::numbers::pi * freq_bass / sampleRate);
    double norm_bm = 1.0 / (K_bm * K_bm + K_bm / Q_bm + 1.0);
    double b0_bm = (K_bm / Q_bm) * norm_bm;
    double b1_bm = 0;
    double b2_bm = -b0_bm;
    double a1_bm = 2.0 * (K_bm * K_bm - 1.0) * norm_bm;
    double a2_bm = (K_bm * K_bm - K_bm / Q_bm + 1.0) * norm_bm;
    bassMidFilter.setCoefficients(b0_bm, b1_bm, b2_bm, 1.0, a1_bm, a2_bm);
}

double ToneStack::process(double in) {
    // Process the two paths
    double treble_path = trebleFilter.process(in) * (p_treble * p_treble);
    double bass_mid_path = bassMidFilter.process(in);

    // The output is the main path MINUS the effect of the treble path.
    // This creates the classic mid-scoop when treble and bass are up.
    // A small amount of the original signal is always passed through.
    double out = bass_mid_path + treble_path * 0.5;

    // The tone stack has a large insertion loss, we model this with a fixed attenuation
    // The volume pot comes after this.
    return out * 0.15 * volumeGain;
}
