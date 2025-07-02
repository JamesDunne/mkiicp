
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

#endif
