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
    // --- 1. Prepare Tube Stages (Simpler: just components) ---
    v1a.prepare(1.5e3, 150e3, 405.0);
    v1b.prepare(1.5e3, 100e3, 405.0);
    v3b.prepare(1.5e3, 82e3,  410.0);
    v4a.prepare(3.3e3, 274e3, 410.0);
    v2a.prepare(1.0e3, 120e3, 410.0);
    v2b.prepare(1.5e3, 100e3, 410.0);

    // --- 2. Prepare Passive Tone Stack ---
    toneStack.prepare(sampleRate);

    // --- 3. Prepare ALL Inter-stage Filters ---
    // Cathode bypass filters are modeled as HPF gain boosts *after* the stage.
    make_hpf(v1a_bypass_hpf, sampleRate, 4.7);
    make_hpf(v1b_bypass_hpf, sampleRate, 4.8);
    make_hpf(v3b_bypass_hpf, sampleRate, 48.2);
    make_hpf(v4a_bypass_hpf, sampleRate, 220.0);
    make_hpf(v2a_bypass_hpf, sampleRate, 10.3);

    // Coupling and shaping filters
    make_lpf(v1a_pre_ts_lpf, sampleRate, 12000.0);
    make_hpf(ts_to_v3b_hpf, sampleRate, 5.3); // C21 into R21+pot
    make_hpf(v3b_to_v4a_hpf, sampleRate, 133.0);
    make_lpf(v3b_to_v4a_lpf, sampleRate, 2270.0);
    make_hpf(v4a_to_mixer_hpf, sampleRate, 22.0);
    make_lpf(mixer_output_lpf, sampleRate, 3400.0);
    make_hpf(v2b_output_hpf, sampleRate, 144.0);

    // Rhythm path specific filter
    make_lpf(v1b_to_rhythm_lpf, sampleRate, 2900.0);
}

void Preamp::setParameters(double treble, double mid, double bass, double vol1, double gain, double master) {
    toneStack.setParams(treble, mid, bass, vol1);
    lead_drive = gain * gain;
    master_vol = master * master;
}

double Preamp::processSample(double in) {
    // --- 1. Define Load Resistances for Each Stage ---
    const double TONE_STACK_R_IN = 1e6;
    // V1B's load is the parallel combination of the rhythm path (R10+R11) and lead path (R9+...) impedances.
    // Let's approximate this complex load.
    const double V1B_COMBINED_LOAD_R_IN = 800e3;
    const double V3B_R_IN = 1.5e6;
    const double V4A_R_IN = 54.3e3;
    const double V2A_R_IN = 47e3;
    const double V2B_R_IN = 1e6;
    const double FINAL_R_OUT = 47e3;

    // --- 2. Common Pre-Gain Section: Input -> V1A -> Tone Stack -> V1B ---

    // Scale the raw audio signal to a realistic guitar level.
    double sample = in * 0.2;

    // Process through the first tube stage (V1A).
    sample = v1a.process(sample, TONE_STACK_R_IN);
    sample = v1a_bypass_hpf.process(sample);
    sample = v1a_pre_ts_lpf.process(sample);
    mm_v1a.measureMinMax(sample);

    // Process through the passive tone stack.
    sample = toneStack.process(sample);
    mm_toneStack.measureMinMax(sample);

    // Process through the second core tube stage (V1B).
    sample = v1b.process(sample, V1B_COMBINED_LOAD_R_IN);
    sample = v1b_bypass_hpf.process(sample);

    // `v1b_out` holds the signal at the crucial split point (node N001),
    // after V1B's amplification but before channel-specific filtering.
    double v1b_out = sample;
    mm_v1b.measureMinMax(sample);

    // --- 3. Rhythm Channel Path ---
    // The rhythm path takes the signal from V1B and applies its specific LPF.
    // This models the network of R10 and C10.
    double rhythm_path = v1b_to_rhythm_lpf.process(v1b_out);

    // --- 4. Lead Channel Path ---
    // The lead path also takes the signal from V1B and sends it to its own chain.

    // The coupling network before V3B (C21, R21) provides attenuation.
    const double ATTENUATION_PRE_V3B = 0.5; // Less attenuation than previously thought
    double lead_path = v1b_out * ATTENUATION_PRE_V3B;

    // Pass through the lead channel's input coupling capacitor filter.
    lead_path = ts_to_v3b_hpf.process(lead_path);

    // Apply the "Lead Drive" knob.
    lead_path *= lead_drive;

    // Process through V3B.
    lead_path = v3b.process(lead_path, V4A_R_IN);
    lead_path = v3b_bypass_hpf.process(lead_path);

    // Process through the V3B -> V4A coupling network.
    lead_path = v3b_to_v4a_hpf.process(lead_path);
    const double R24 = 68e3, R25 = 270e3;
    lead_path *= R24 / (R25 + R24); // Voltage divider attenuation.
    lead_path = v3b_to_v4a_lpf.process(lead_path); // "Fizz" filter.
    mm_v3b.measureMinMax(lead_path);

    // Process through V4A.
    lead_path = v4a.process(lead_path, V2A_R_IN);
    lead_path = v4a_bypass_hpf.process(lead_path);
    lead_path = v4a_to_mixer_hpf.process(lead_path);
    mm_v4a.measureMinMax(lead_path);

    // --- 5. Mixdown and Final Driver Stages ---

    // The resistive mixer network attenuates the signals before they are summed.
    const double MIXER_ATTENUATION = 0.5;

    // Blend the two paths. Here, we model being mostly in the lead channel.
    // The rhythm path signal is taken from its final point, the lead path from its.
    sample = (rhythm_path * (1.0 - 0.8) + lead_path * 0.8) * MIXER_ATTENUATION;

    // Apply the final post-mixer low-pass filter (R11/C11/etc.).
    sample = mixer_output_lpf.process(sample);

    // Process through V2A.
    sample = v2a.process(sample, V2B_R_IN);
    sample = v2a_bypass_hpf.process(sample);

    // Apply the Master Volume control.
    sample *= master_vol;
    mm_v2a.measureMinMax(sample);

    // Process through the final buffer/driver stage, V2B (unbypassed).
    sample = v2b.process(sample, FINAL_R_OUT);

    // Apply the final output coupling capacitor filter.
    sample = v2b_output_hpf.process(sample);
    mm_output.measureMinMax(sample);

    // --- 6. Final Output Attenuation ---
    return sample * 0.00125;
}
