#include "Preamp.h"
#include <numbers>

constexpr double PI = std::numbers::pi;

// --- Utility function to create a 1-pole High Pass Filter ---
void make_hpf(IIRBiquad& filter, double sampleRate, double cutoff_freq) {
    if (cutoff_freq <= 0.0) { // Avoid issues with zero or negative frequencies
        filter.setCoefficients(1.0, 0.0, 0.0, 1.0, 0.0, 0.0); // Pass-through
        return;
    }
    double wc = 2.0 * PI * cutoff_freq;
    double K = tan(wc / (2.0 * sampleRate));
    double norm = 1.0 / (K + 1.0);

    filter.setCoefficients(
        1.0 * norm,      // b0
        -1.0 * norm,     // b1
        0.0,             // b2
        1.0,             // a0 (normalized)
        (K - 1.0) * norm,// a1
        0.0              // a2
    );
}

// --- V1A ---
void V1AStage::prepare(double sampleRate) {
    make_hpf(cathodeBypassFilter, sampleRate, 4.7);
    make_hpf(outputCouplingFilter, sampleRate, 31.0);
    reset();
}

void V1AStage::reset() {
    cathodeBypassFilter.reset();
    outputCouplingFilter.reset();
}

double V1AStage::process(double in) {
    double out = cathodeBypassFilter.process(in);
    out = tubeSaturate(out, gain);
    out = outputCouplingFilter.process(out);
    return out;
}

// --- Tone Stack ---
void ToneStack::prepare(double sampleRate) {
    this->sampleRate = sampleRate;
    calculateCoefficients();
    reset();
}

void ToneStack::reset() {
    biquad.reset();
    onePole.reset();
}

void ToneStack::setParams(double treble, double mid, double bass, double volume) {
    treble_p = treble;
    mid_p = mid;
    bass_p = bass;
    volume_p = volume;
    calculateCoefficients();
}

void ToneStack::calculateCoefficients() {
    // This is a highly simplified model of the Mesa tone stack's transfer function.
    // Full symbolic derivation is complex. This captures the general behavior.
    double R5 = 100e3;
    double C4 = 0.1e-6;
    double C3 = 0.047e-6;

    double treble_pot = 250e3 * treble_p * treble_p;
    double bass_pot = 250e3 * bass_p * bass_p;
    double mid_pot = 10e3 * mid_p * mid_p;

    double v = pow(10, volume_p - 1.0);

    double c0 = C3 * C4 * mid_pot * R5;
    double c1 = (C3 + C4) * mid_pot + C4 * R5 + C3 * bass_pot;
    double c2 = mid_pot + bass_pot;
    double c3 = C3 * C4 * R5 * treble_pot;
    double c4 = C4 * R5 + C3 * treble_pot;

    double b0_a = c0;
    double b1_a = c1 * 0.5;
    double b2_a = c2 * 0.1;
    double a0_a = c3 + 1e-7;
    double a1_a = c4 * 0.5;
    double a2_a = c2 * 0.2 + treble_pot * 0.1;

    double T = 1.0 / sampleRate;
    double K = tan(PI * 1000 / sampleRate);
    double K2 = K * K;

    double b0 = b0_a*K2 + b1_a*K + b2_a;
    double b1 = 2 * (b2_a - b0_a*K2);
    double b2 = b0_a*K2 - b1_a*K + b2_a;
    double a0 = a0_a*K2 + a1_a*K + a2_a;
    double a1 = 2 * (a2_a - a0_a*K2);
    double a2 = a0_a*K2 - a1_a*K + a2_a;

    biquad.setCoefficients(b0, b1, b2, a0, a1, a2);

    double R_treb_vol = treble_pot + 1e6 * v;
    make_hpf(onePole, sampleRate, 1.0 / (2 * PI * C4 * R_treb_vol + 1e-9));
}

double ToneStack::process(double in) {
    double out = biquad.process(in);
    out = onePole.process(out);
    return out * volume_p * volume_p;
}

// --- V1B ---
void V1BStage::prepare(double sampleRate) {
    make_hpf(cathodeBypassFilter, sampleRate, 4.8);
    make_hpf(outputCouplingFilter, sampleRate, 0.4);
    reset();
}

void V1BStage::reset() {
    cathodeBypassFilter.reset();
    outputCouplingFilter.reset();
}

double V1BStage::process(double in) {
    double out = cathodeBypassFilter.process(in);
    out = tubeSaturate(out, gain);
    out = outputCouplingFilter.process(out);
    return out;
}

// --- V3B and V4A (Lead Gain) ---
void V3BV4AStage::prepare(double sampleRate) {
    make_hpf(v3b_inputFilter, sampleRate, 5.0);
    make_hpf(v3b_cathodeBypass, sampleRate, 48.2);
    make_hpf(v3b_outputCoupling, sampleRate, 26.0);
    make_hpf(v4a_cathodeBypass, sampleRate, 220.0);
    make_hpf(v4a_outputCoupling, sampleRate, 15.0);
    reset();
}

void V3BV4AStage::reset() {
    v3b_inputFilter.reset();
    v3b_cathodeBypass.reset();
    v3b_outputCoupling.reset();
    v4a_cathodeBypass.reset();
    v4a_outputCoupling.reset();
}

void V3BV4AStage::setGain(double g) {
    drive = g * g;
}

double V3BV4AStage::process(double in) {
    double out = v3b_inputFilter.process(in * drive);
    out = v3b_cathodeBypass.process(out);
    out = tubeSaturate(out, v3b_gain);
    out = v3b_outputCoupling.process(out);
    out = v4a_cathodeBypass.process(out);
    out = tubeSaturate(out, v4a_gain);
    out = v4a_outputCoupling.process(out);
    return out;
}

// --- V2A ---
void V2AStage::prepare(double sampleRate) {
    make_hpf(cathodeBypassFilter, sampleRate, 10.3);
    make_hpf(outputCouplingFilter, sampleRate, 3.4);
    reset();
}

void V2AStage::reset() {
    cathodeBypassFilter.reset();
    outputCouplingFilter.reset();
}

void V2AStage::setMaster(double m) {
    masterVol = m * m;
}

double V2AStage::process(double in) {
    double out = cathodeBypassFilter.process(in);
    out = tubeSaturate(out, gain);
    out = outputCouplingFilter.process(out);
    return out * masterVol;
}

// --- V2B ---
void V2BStage::prepare(double sampleRate) {
    make_hpf(outputCouplingFilter, sampleRate, 144.0);
    reset();
}

void V2BStage::reset() {
    outputCouplingFilter.reset();
}

double V2BStage::process(double in) {
    double out = tubeSaturate(in, gain * 0.5);
    out = outputCouplingFilter.process(out);
    return out;
}
