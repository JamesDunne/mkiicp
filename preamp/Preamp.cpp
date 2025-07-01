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

double Preamp::processSample(double in) {
    double input_voltage = in * 0.2;

    // Input Stage -> Tone Stack
    double sample = v1a.process(input_voltage);
    sample = v1a.outputFilter.process(sample);
    sample = v1a.interStageLPF.process(sample);
    mm_v1a.measureMinMax(sample);
    sample = toneStack.process(sample);
    mm_toneStack.measureMinMax(sample);

    // Rhythm Path
    double rhythm_path = v1b.inputFilter.process(sample);
    rhythm_path = v1b.process(rhythm_path);
    rhythm_path = v1b.interStageLPF.process(rhythm_path);

    // Lead Path
    double lead_path = sample * lead_drive;
    lead_path = v3b.inputFilter.process(lead_path);
    lead_path = v3b.process(lead_path);

    lead_path = v4a.inputFilter.process(lead_path);
    lead_path = v4a.process(lead_path);
    lead_path = v4a.interStageLPF.process(lead_path);

    // Mixdown and Final Stages
    sample = lead_path * 0.8 + rhythm_path * 0.2;
    sample = v2a.inputFilter.process(sample);
    sample = v2a.process(sample) * master_vol;
    
    sample = v2b.inputFilter.process(sample);
    sample = v2b.process(sample);
    sample = v2b.outputFilter.process(sample);

    sample /= 600.0;
    mm_output.measureMinMax(sample);

    // Final attenuation to bring virtual voltage back to audio level
    return sample;
}
