
#include <numbers>

#include "Preamp.h"
#include "ToneStack.h"
#include "IIRFilter.h"

constexpr double PI = std::numbers::pi;

/**
 * @brief Final, correctly normalized 12AX7 waveshaper with variable plate voltage.
 *
 * @param x The input signal, [-1, 1].
 * @param drive Controls how much of the input signal is used.
 * @param bias The cathode bias point in Volts. Typical: -1.5V.
 * @param Vp The DC plate voltage of the specific tube stage.
 * @return The distorted, normalized output signal.
 */
inline double tubeSaturate12AX7(double x, double drive, double bias, double Vp) {
    // --- Parameters ---
    const double Mu = 96.20;
    const double Ex = 1.437;

    // 1. Calculate Grid Voltage (Vgk). 'drive' scales the input signal.
    double Vgk = x * drive + bias;

    // 2. Grid Conduction.
    if (Vgk > 0.0) {
        Vgk = 0.5 * std::tanh(Vgk * 2.0);
    }

    // 3. Plate Current Calculation using the new Vp parameter.
    double Vp_eff = (Vp / Mu) + Vgk;
    double plate_current = (Vp_eff > 0.0) ? std::pow(Vp_eff, Ex) : 0.0;

    // 4. Calculate Quiescent Current to make the output AC-coupled (DC-free)
    double Vp_eff_quiescent = (Vp / Mu) + bias;
    double quiescent_current = (Vp_eff_quiescent > 0.0) ? std::pow(Vp_eff_quiescent, Ex) : 0.0;

    // 5. Calculate output and normalize
    double output = -(plate_current - quiescent_current);

    return output * 0.4;
}

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

// --- V1A ---
void V1AStage::prepare(double sampleRate) {
    make_hpf(cathodeBypassFilter, sampleRate, 4.7);
    make_hpf(outputCouplingFilter, sampleRate, 1.6);
    make_lpf(interStageLPF, sampleRate, 12000.0); // Gentle roll-off
    gain = 50.0;
    reset();
}

void V1AStage::reset() {
    cathodeBypassFilter.reset();
    outputCouplingFilter.reset();
    interStageLPF.reset();
}

double V1AStage::process(double in) {
    double out = cathodeBypassFilter.process(in);
    out = tubeSaturate12AX7(out, 1.0, -1.8, 230.0);
    out *= gain;
    out = outputCouplingFilter.process(out);
    out = interStageLPF.process(out);
    return out;
}

// --- V1B ---
void V1BStage::prepare(double sampleRate) {
    make_hpf(cathodeBypassFilter, sampleRate, 4.8);
    make_hpf(outputCouplingFilter, sampleRate, 0.4);
    make_lpf(interStageLPF, sampleRate, 2900.0); // Initialize the LPF
    // This stage expects the post-tonestack signal, so its gain can be modest.
    gain = 25.0;
    reset();
}

void V1BStage::reset() {
    cathodeBypassFilter.reset();
    outputCouplingFilter.reset();
    interStageLPF.reset();
}

double V1BStage::process(double in) {
    double out = cathodeBypassFilter.process(in);
    out = tubeSaturate12AX7(out, 1.0, -1.8, 285.0);
    out = outputCouplingFilter.process(out);
    out *= gain;
    // Apply the LPF at the very end of the stage's output
    out = interStageLPF.process(out);
    return out;
}

// --- V3B and V4A (Lead Gain) ---
void V3BV4AStage::prepare(double sampleRate) {
    make_hpf(v3b_inputFilter, sampleRate, 5.3);
    make_hpf(v3b_cathodeBypass, sampleRate, 48.2);
    make_hpf(interStageHPF, sampleRate, 133.0);
    make_lpf(fizzFilter, sampleRate, 3180.0);
    make_hpf(v4a_cathodeBypass, sampleRate, 220.0);
    make_hpf(v4a_outputCoupling, sampleRate, 15.0);
    reset();
}

void V3BV4AStage::reset() {
    v3b_inputFilter.reset();
    v3b_cathodeBypass.reset();
    interStageHPF.reset();
    fizzFilter.reset();
    v4a_cathodeBypass.reset();
    v4a_outputCoupling.reset();
}

void V3BV4AStage::setGain(double g) {
    // Maps the [0,1] knob to a powerful gain factor [1x to 51x].
    // This is the core of the lead drive.
    drive = 1.0 + 150.0 * (g * g);
}

double V3BV4AStage::process(double in) {
    // V3B acts as a clean booster whose output level is controlled by the 'drive' knob.
    double out = v3b_inputFilter.process(in);
    out = v3b_cathodeBypass.process(out);
    out = interStageHPF.process(out);
    // v3b_plateVoltage = 300.0
    out *= drive; // Apply the main lead gain here.
    out = fizzFilter.process(out);

    // V4A gets slammed by the high-level signal from V3B.
    // This is where the heavy clipping happens.
    out = v4a_cathodeBypass.process(out);
    out = tubeSaturate12AX7(out, 1.0, -2.0, 218.0);
    out = v4a_outputCoupling.process(out);

    return out;
}

// --- V2A ---
void V2AStage::prepare(double sampleRate) {
    make_hpf(inputCouplingFilter, sampleRate, 72.0);
    make_hpf(cathodeBypassFilter, sampleRate, 10.3);
    make_hpf(outputCouplingFilter, sampleRate, 3.4);
    gain = 15.0;
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
    double out = inputCouplingFilter.process(in);

    out = cathodeBypassFilter.process(out);
    out = tubeSaturate12AX7(out, 1.0, -1.5, 230.0);
    out *= gain;
    out = outputCouplingFilter.process(out);
    return out * masterVol;
}

// --- V2B ---
void V2BStage::prepare(double sampleRate) {
    make_hpf(outputCouplingFilter, sampleRate, 144.0);
    gain = 1.0; // Unity gain buffer.
    reset();
}

void V2BStage::reset() {
    outputCouplingFilter.reset();
}

double V2BStage::process(double in) {
    double out = tubeSaturate12AX7(in, 1.0, -1.5, 290.0);
    out = outputCouplingFilter.process(out);
    return out * gain;
}

double Preamp::processSample(double in) {
    // Process through the chain
    double sample = v1a.process(in);

    sample = toneStack.process(sample);

#if 0
    double lead_path = v3b_v4a.process(sample);

    double rhythm_path = v1b.process(sample);

    // Mix the paths.
    sample = lead_path * 0.8 + rhythm_path * 0.2;

    sample = v2a.process(sample);
    sample = v2b.process(sample);

    sample *= 0.25;
#else
    sample /= 2.5;
#endif
    measureMinMax(sample);

    return sample;
}
