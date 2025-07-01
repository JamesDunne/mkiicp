#include "Preamp.h"
#include <cmath>
#include <numbers>

// --- Utility Functions ---
void make_lpf(IIRBiquad& filter, double sampleRate, double cutoff_freq) {
    if (cutoff_freq <= 0.0) return;
    double K = tan(std::numbers::pi * cutoff_freq / sampleRate);
    double norm = 1.0 / (1.0 + K);
    filter.setCoefficients(K * norm, K * norm, 0.0, 1.0, (K - 1.0) * norm, 0.0);
}

void make_hpf(IIRBiquad& filter, double sampleRate, double cutoff_freq) {
    if (cutoff_freq <= 0.0) return;
    double K = tan(std::numbers::pi * cutoff_freq / sampleRate);
    double norm = 1.0 / (K + 1.0);
    filter.setCoefficients(norm, -norm, 0.0, 1.0, (K - 1.0) * norm, 0.0);
}

// --- Preamp Implementation ---
Preamp::Preamp() {
    setParameters(0.8, 0.5, 0.25, 0.75, 0.75, 0.5);
}

void Preamp::setParameters(double treble, double mid, double bass, double vol1, double gain, double master) {
    toneStack.setParams(treble, mid, bass, vol1);
    lead_drive = gain * gain;
    master_vol = master * master;
}

void Preamp::prepare(double sampleRate) {
    // --- Prepare Tube Stages (with bypass cap info) ---
    v1a.prepare(sampleRate, 1.5e3, 150e3, 405.0, 22.47e-6);
    v1b.prepare(sampleRate, 1.5e3, 100e3, 405.0, 22.0e-6);
    v3b.prepare(sampleRate, 1.5e3, 82e3,  410.0, 2.2e-6);
    v4a.prepare(sampleRate, 3.3e3, 274e3, 410.0, 0.22e-6);
    v2a.prepare(sampleRate, 1.0e3, 120e3, 410.0, 15.47e-6);
    v2b.prepare(sampleRate, 1.5e3, 100e3, 410.0, 0.0);

    // --- Prepare Tone Stack (and prime its DC blocker) ---
    toneStack.prepare(sampleRate, v1a.getVpDC());

    // --- Prepare ALL Inter-stage Filters ---
    make_lpf(v1a_output_lpf, sampleRate, 12000.0);
    make_lpf(v1b_to_rhythm_lpf, sampleRate, 2900.0);
    make_hpf(ts_to_v3b_hpf, sampleRate, 5.3);
    make_hpf(v3b_to_v4a_hpf, sampleRate, 133.0);
    make_lpf(v3b_to_v4a_lpf, sampleRate, 2270.0);
    make_hpf(v4a_to_mixer_hpf, sampleRate, 22.0);
    make_lpf(mixer_output_lpf, sampleRate, 3400.0);
    make_hpf(v2b_output_hpf, sampleRate, 144.0);
}

double Preamp::processSample(double in) {
    double sample = in;

    // V1A -> Tone Stack
    sample = v1a.process(sample);
    mm_v1a.measureMinMax(sample);
#if 1
    return sample;
#else
    sample = v1a_output_lpf.process(sample);

    sample = toneStack.process(sample);
    mm_toneStack.measureMinMax(sample);
    sample = v1b.process(sample);
    mm_v1b.measureMinMax(sample);

    // Rhythm Path
    double rhythm_path = v1b_to_rhythm_lpf.process(sample);

    // Lead Path
    const double ATTENUATION_PRE_V3B = 0.2;
    double lead_path = sample * ATTENUATION_PRE_V3B;
    lead_path = ts_to_v3b_hpf.process(lead_path);
    lead_path *= lead_drive;
    lead_path = v3b.process(lead_path);
    lead_path = v3b_to_v4a_hpf.process(lead_path);
    lead_path *= 68e3 / (270e3 + 68e3);
    lead_path = v3b_to_v4a_lpf.process(lead_path);
    lead_path = v4a.process(lead_path);
    lead_path = v4a_to_mixer_hpf.process(lead_path);
    
    // Mixdown & Final Stages
    const double MIXER_ATTENUATION = 0.5;
    sample = (rhythm_path * 0.2 + lead_path * 0.8) * MIXER_ATTENUATION;
    sample = mixer_output_lpf.process(sample);
    sample = v2a.process(sample);
    sample *= master_vol;
    sample = v2b.process(sample);
    sample = v2b_output_hpf.process(sample);
    mm_output.measureMinMax(sample);
    
    return sample / 400.0;
#endif
}
