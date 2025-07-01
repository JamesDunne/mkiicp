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
    // --- 1. Prepare Individual Tube Stages with Physical Components ---
    // Sets up the internal DC solver for each tube based on its unique environment.
    // format: prepare(sampleRate, R_cathode, R_plate_load, V_B+_supply)
    v1a.prepare(sampleRate, 1.5e3, 150e3, 405.0); // XV1A: R2=1.5k, R4=150k, VE=405V
    v1b.prepare(sampleRate, 1.5e3, 100e3, 405.0); // XV1B: R7=1.5k, R8=100k, VE=405V
    v3b.prepare(sampleRate, 1.5e3, 82e3,  410.0); // XV3B: R23=1.5k, R26=82k, VC=410V
    v4a.prepare(sampleRate, 3.3e3, 274e3, 410.0); // XV4A: R30=3.3k, R27=274k, VC=410V
    v2a.prepare(sampleRate, 1.0e3, 120e3, 410.0); // XV2A: R104=1k, R19=120k, VC2=410V
    v2b.prepare(sampleRate, 1.5e3, 100e3, 410.0); // XV2B: R16=1.5k, R13=100k, VC2=410V


    // --- 2. Configure Cathode Bypass Filters (HPF) for each stage ---
    // These create a low-frequency boost, increasing gain at non-bass frequencies.
    make_hpf(v1a.cathodeBypass, sampleRate, 4.7);   // C1+C2 = 22.47uF across R2=1.5k
    make_hpf(v1b.cathodeBypass, sampleRate, 4.8);   // C13 = 22uF across R7=1.5k
    make_hpf(v3b.cathodeBypass, sampleRate, 48.2);  // C23 = 2.2uF across R23=1.5k
    make_hpf(v4a.cathodeBypass, sampleRate, 220.0); // C29 = 0.22uF across R30=3.3k (Lead tone shaping)
    make_hpf(v2a.cathodeBypass, sampleRate, 10.3);  // C15+C16 = 15.47uF across R104=1k
    // V2B (XV2B) has an unbypassed cathode (R16), so it has no cathodeBypass filter.


    // --- 3. Configure Inter-stage Coupling and Filtering ---
    // These HPFs and LPFs shape the sound between tube stages.

    // V1A Output -> Tone Stack Input
    // Coupling cap C4=.1u into Tone Stack (~1Meg impedance)
    make_hpf(v1a.outputFilter, sampleRate, 1.6);
    // Gentle LPF for Miller effect simulation
    make_lpf(v1a.interStageLPF, sampleRate, 12000.0);

    // Tone Stack Output -> V1B (Rhythm Path)
    // C7=.1u into ~4Meg impedance (R10+R11)
    make_hpf(v1b.inputFilter, sampleRate, 0.4);
    // LPF on the grid of V2B (fed by V1B) from C11+C32 across R8
    make_lpf(v1b.interStageLPF, sampleRate, 2900.0);

    // Tone Stack Output -> V3B (Lead Path)
    // C21=.02u into >1.5Meg impedance
    make_hpf(v3b.inputFilter, sampleRate, 5.3);

    // V3B Output -> V4A Input
    // C25=.022u into ~54.3k impedance (R25 || R24) - Critical lead tone shaping
    make_hpf(v4a.inputFilter, sampleRate, 133.0);
    // "Fizz" LPF on V4A grid: C24=1000p across ~50k impedance
    make_lpf(v4a.interStageLPF, sampleRate, 3180.0);

    // V4A Output -> V2A Input
    // C30=.047u into R103=47k
    make_hpf(v2a.inputFilter, sampleRate, 72.0);

    // V2A Output -> V2B Input
    // This is a complex path. The main coupling is via C12=.047u into the Master pot.
    // For simplicity, we model the output coupling of V2A into the V2B stage.
    make_hpf(v2b.inputFilter, sampleRate, 3.4); // C12 into 1Meg Master pot

    // V2B Output -> Final Output Jack
    // C9=.047u into R105/R46 divider (~23.5k)
    make_hpf(v2b.outputFilter, sampleRate, 144.0);

    // --- 4. Prepare Passive Tone Stack ---
    toneStack.prepare(sampleRate);
}

void Preamp::setParameters(double treble, double mid, double bass, double vol1, double gain, double master) {
    toneStack.setParams(treble, mid, bass, vol1);
    lead_drive = gain * gain;
    master_vol = master * master;
}

double Preamp::processSample(double in) {
    // The signal flow is now a cascade of amplification and attenuation stages
    
    // Input Stage -> Tone Stack
    double sample = v1a.cathodeBypass.process(in);
    sample = v1a.process(sample);
    sample = v1a.outputFilter.process(sample);
    sample = v1a.interStageLPF.process(sample);
    mm_v1a.measureMinMax(sample);
    sample /= 1000.0;
#if 0

    sample = toneStack.process(sample);
    mm_toneStack.measureMinMax(sample);

    // Rhythm Path
    double rhythm_path = v1b.cathodeBypass.process(sample);
    rhythm_path = v1b.process(rhythm_path);
    rhythm_path = v1b.interStageLPF.process(rhythm_path);
    
    // Lead Path
    double lead_path = sample * lead_drive; // Drive knob attenuates signal before V3B
    lead_path = v3b.inputFilter.process(lead_path);
    lead_path = v3b.cathodeBypass.process(lead_path);
    lead_path = v3b.process(lead_path); // V3B adds gain/distortion
    
    lead_path = v4a.inputFilter.process(lead_path); // Models coupling network
    lead_path = v4a.cathodeBypass.process(lead_path);
    lead_path = v4a.process(lead_path); // V4A adds more gain/distortion

    // Mixdown and Final Stages
    sample = lead_path * 0.8 + rhythm_path * 0.2;
    sample = v2a.inputFilter.process(sample);
    sample = v2a.cathodeBypass.process(sample);
    sample = v2a.process(sample) * master_vol; // Master volume attenuates
    
    sample = v2b.inputFilter.process(sample);
    sample = v2b.process(sample);
    sample = v2b.outputFilter.process(sample);
#endif

    // Final attenuation to bring virtual voltage back to audio level
    return sample;
}
