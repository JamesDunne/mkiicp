// ... includes ...

#include <numbers>

#include "Preamp.h"

constexpr double PI = std::numbers::pi;

void make_lpf(IIRBiquad& filter, double sampleRate, double cutoff_freq) {
    if (cutoff_freq >= sampleRate / 2.0) { // Avoid aliasing issues
        filter.setCoefficients(1.0, 0.0, 0.0, 1.0, 0.0, 0.0); // Pass-through
        return;
    }
    double wc = 2.0 * PI * cutoff_freq;
    double K = tan(wc / (2.0 * sampleRate));
    double norm = 1.0 / (K + 1.0);

    filter.setCoefficients(
        K * norm,     // b0
        K * norm,     // b1
        0.0,          // b2
        1.0,          // a0 (normalized)
        (K - 1.0) * norm, // a1
        0.0           // a2
    );
}

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

Preamp::Preamp() {
    setParameters(0.8, 0.5, 0.25, 0.75, 0.75, 0.5);
}

void Preamp::prepare(double sampleRate) {
    // --- 1. Prepare Tube Stages with Physical Components AND Bypass Caps ---
    // format: prepare(sampleRate, R_k,   R_L,   V_supply, C_k)
    v1a.prepare(sampleRate, 1.5e3, 150e3, 405.0, 22.47e-6); // XV1A
    v1b.prepare(sampleRate, 1.5e3, 100e3, 405.0, 22.0e-6);  // XV1B
    v3b.prepare(sampleRate, 1.5e3, 82e3,  410.0, 2.2e-6);   // XV3B
    v4a.prepare(sampleRate, 3.3e3, 274e3, 410.0, 0.22e-6);  // XV4A
    v2a.prepare(sampleRate, 1.0e3, 120e3, 410.0, 15.47e-6); // XV2A
    v2b.prepare(sampleRate, 1.5e3, 100e3, 410.0, 0.0);      // XV2B (Unbypassed)

    // --- 2. Prepare Passive Tone Stack ---
    toneStack.prepare(sampleRate);

    // --- 3. Prepare ALL Inter-stage Filters ---
    // These filters model the passive networks BETWEEN the tube stages.

    // V1A Output -> Tone Stack
    make_lpf(v1a_pre_ts_lpf, sampleRate, 12000.0); // Simulates Miller effect/stray capacitance

    // Rhythm Path Specific Filter (after V1B)
    make_lpf(v1b_to_rhythm_lpf, sampleRate, 2900.0); // R10/C10 network

    // Lead Path Coupling Filters
    make_hpf(ts_to_v3b_hpf, sampleRate, 5.3); // C21 into R21+pot network
    make_hpf(v3b_to_v4a_hpf, sampleRate, 133.0);
    make_lpf(v3b_to_v4a_lpf, sampleRate, 2270.0);
    make_hpf(v4a_to_mixer_hpf, sampleRate, 22.0);

    // Mixer and Final Stage Filters
    make_lpf(mixer_output_lpf, sampleRate, 3400.0); // R11/C11 network
    make_hpf(v2b_output_hpf, sampleRate, 144.0); // Final output coupling cap
}

void Preamp::setParameters(double treble, double mid, double bass, double vol1, double gain, double master) {
    toneStack.setParams(treble, mid, bass, vol1);
    lead_drive = gain * gain;
    master_vol = master * master;
}

double Preamp::processSample(double in) {
    // --- 1. Common Pre-Gain Section: Input -> V1A -> Tone Stack -> V1B ---
    double sample = in * 0.2; // Scale input to realistic guitar level

    sample = v1a.process(sample);
    sample = v1a_pre_ts_lpf.process(sample);
    mm_v1a.measureMinMax(sample);
    sample = toneStack.process(sample);
    mm_toneStack.measureMinMax(sample);

#if 1
    return sample / 80.0;
#else
    sample = v1b.process(sample);
    double v1b_out = sample; // Signal at the split point (N001)

    // --- 2. Rhythm & Lead Path Processing ---
    double rhythm_path = v1b_to_rhythm_lpf.process(v1b_out);

    const double ATTENUATION_PRE_V3B = 0.2;
    double lead_path = v1b_out * ATTENUATION_PRE_V3B;
    lead_path = ts_to_v3b_hpf.process(lead_path);
    lead_path *= lead_drive;

    lead_path = v3b.process(lead_path);

    lead_path = v3b_to_v4a_hpf.process(lead_path);
    const double R24 = 68e3, R25 = 270e3;
    lead_path *= R24 / (R25 + R24);
    lead_path = v3b_to_v4a_lpf.process(lead_path);

    lead_path = v4a.process(lead_path);
    lead_path = v4a_to_mixer_hpf.process(lead_path);

    // --- 3. Mixdown & Final Stages ---
    const double MIXER_ATTENUATION = 0.5;
    sample = (rhythm_path * (1.0 - 0.8) + lead_path * 0.8) * MIXER_ATTENUATION;
    sample = mixer_output_lpf.process(sample);

    sample = v2a.process(sample);
    sample *= master_vol;

    sample = v2b.process(sample); // V2B is unbypassed
    sample = v2b_output_hpf.process(sample);

    return sample * 0.00125;
#endif
}
