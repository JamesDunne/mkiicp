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
    v1a.prepare(sampleRate, 1.5e3, 150e3, 405.0, 22.47e-6); // C1+C2
    v1b.prepare(sampleRate, 1.5e3, 100e3, 405.0, 22.0e-6);  // C13
    v3b.prepare(sampleRate, 1.5e3, 82e3,  410.0, 2.2e-6);   // C23
    v4a.prepare(sampleRate, 3.3e3, 274e3, 410.0, 0.22e-6);  // C29
    v2a.prepare(sampleRate, 1.0e3, 120e3, 410.0, 15.47e-6); // C15+C16
    v2b.prepare(sampleRate, 1.5e3, 100e3, 410.0, 0.0);      // No bypass cap

    // --- 2. Configure Inter-stage Coupling Filters ---
    // This logic remains the same, but now acts on the public filter members.
    make_hpf(v1a.outputFilter, sampleRate, 1.6);
    make_lpf(v1a.interStageLPF, sampleRate, 12000.0);

    make_hpf(v1b.inputFilter, sampleRate, 0.4);
    make_lpf(v1b.interStageLPF, sampleRate, 2900.0);

    make_hpf(v3b.inputFilter, sampleRate, 5.3);

    make_hpf(v4a.inputFilter, sampleRate, 133.0);
    make_lpf(v4a.interStageLPF, sampleRate, 3180.0);

    make_hpf(v2a.inputFilter, sampleRate, 72.0);

    make_hpf(v2b.inputFilter, sampleRate, 3.4);
    make_hpf(v2b.outputFilter, sampleRate, 144.0);

    // --- 3. Prepare Passive Tone Stack ---
    toneStack.prepare(sampleRate);
}

void Preamp::setParameters(double treble, double mid, double bass, double vol1, double gain, double master) {
    toneStack.setParams(treble, mid, bass, vol1);
    lead_drive = gain * gain;
    master_vol = master * master;
}

// In Preamp.cpp

double Preamp::processSample(double in) {
    // --- Define Load Resistances for Each Stage ---
    // These values are the input impedance of the *next* block in the chain.
    const double TONE_STACK_R_IN = 1e6; // Tone stack input impedance
    const double V1B_R_IN = 4e6;        // V1B input impedance (R10+R11)
    const double V3B_R_IN = 1e6 * (1.0 - lead_drive) + 680e3; // Approx. Z_in of drive pot + R21
    const double V4A_R_IN = 54.3e3;     // V4A input impedance (R25 || R24)
    const double V2A_R_IN = 47e3;       // V2A input impedance (R103)
    const double V2B_R_IN = 1e6;        // V2B input impedance (Master Pot)
    const double FINAL_R_OUT = 47e3;    // Final output load

    // --- 1. Input Stage -> Tone Stack -> V1B ---
    double sample = in * 0.2; // Scale input to realistic guitar level

    sample = v1a.process(sample, TONE_STACK_R_IN);
    sample = v1a.outputFilter.process(sample);
    sample = v1a.interStageLPF.process(sample);
    mm_v1a.measureMinMax(sample);

    sample = toneStack.process(sample);
    mm_toneStack.measureMinMax(sample);

    sample = v1b.inputFilter.process(sample);
    sample = v1b.process(sample, V3B_R_IN); // V1B drives the lead stage network
    sample = v1b.interStageLPF.process(sample);
    mm_v1b.measureMinMax(sample);

    // At this point, `sample` represents the signal at the output of V1B's
    // post-filtering network. This is the source for both rhythm and lead paths.
    double rhythm_path_out = sample;

    // --- 2. Lead Path with Correct Inter-stage Attenuation ---

    // --- V1B Output to V3B Input Network ---
    // The signal from V1B is attenuated by R9 (100k) and R10+R11 (4M) before C21.
    // And then AGAIN by R21 (680k) and the Lead Drive pot.
    // This is a significant signal loss, which we model here.
    const double ATTENUATION_PRE_DRIVE = 0.1; // Approx. 10x signal loss
    double lead_path = sample * ATTENUATION_PRE_DRIVE;

    // The Lead Drive knob itself is a variable attenuator.
    // A higher setting means LESS attenuation (more signal).
    lead_path *= lead_drive;

    // Process through V3B
    lead_path = v3b.inputFilter.process(lead_path);
    lead_path = v3b.process(lead_path, V4A_R_IN);
    mm_v3b.measureMinMax(lead_path);

    // --- V3B Output to V4A Input Network ---
    // This is a CRITICAL voltage divider formed by R25 and R24.
    // V(N030) = V(N018) * R24 / (R25 + R24)
    const double R24 = 68e3;
    const double R25 = 270e3;
    const double V4A_DIVIDER_ATTENUATION = R24 / (R25 + R24); // Approx 0.2

    lead_path *= V4A_DIVIDER_ATTENUATION;

    // Process through V4A
    lead_path = v4a.inputFilter.process(lead_path);
    lead_path = v4a.process(lead_path, V2A_R_IN);
    lead_path = v4a.interStageLPF.process(lead_path); // "Fizz" filter
    mm_v4a.measureMinMax(lead_path);

    double lead_path_out = lead_path;

    // --- 3. Mixdown and Final Stages ---

    // --- Mixer Simulation ---
    // The outputs of the rhythm path (at V1B) and lead path (at V4A) are mixed
    // by a resistive network (R10, R11, R31, R32). This also causes attenuation.
    const double MIXER_ATTENUATION = 0.5;

    // Blend the attenuated signals
    sample = (rhythm_path_out * (1.0 - 0.8) + lead_path_out * 0.8) * MIXER_ATTENUATION;

    // Process through V2A
    sample = v2a.inputFilter.process(sample);
    sample = v2a.process(sample, V2B_R_IN);
    sample *= master_vol;
    mm_v2a.measureMinMax(sample);

    // Process through V2B
    sample = v2b.inputFilter.process(sample);
    sample = v2b.process(sample, FINAL_R_OUT);
    sample = v2b.outputFilter.process(sample);
    mm_output.measureMinMax(sample);

    // --- 4. Final Output ---
    return sample * 0.00125; // Final trim may need adjustment
}
