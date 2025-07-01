
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

/**
 * @brief An advanced 12AX7 waveshaper that outputs a virtual AC plate voltage.
 *
 * This function models a tube stage based on its physical components, calculating
 * the plate current and producing an output voltage swing based on the plate load resistor.
 * This correctly models the inherent gain of the stage.
 *
 * @param x The input AC "virtual voltage" at the grid.
 * @param bias The DC cathode bias point in Volts (e.g., -1.8V).
 * @param V_supply The B+ supply voltage for the stage.
 * @param R_L The plate load resistor value in Ohms.
 * @return The resulting AC "virtual voltage" at the plate.
 */
inline double tubeSaturateVirtualVoltage(double x, double bias, double V_supply, double R_L) {
    // --- Tube Parameters from SPICE Model ---
    const double Mu = 96.20;  // Amplification Factor
    const double Ex = 1.437;  // Power Law Exponent
    const double KG1 = 613.4; // Plate Current Scaling Factor

    // --- 1. Calculate Quiescent (DC) State ---
    // Effective grid voltage at idle (no signal)
    double Vgk_quiescent = bias;
    double Vp_eff_quiescent = (V_supply / Mu) + Vgk_quiescent;
    // Plate current at idle
    double Ip_quiescent = 0.0;
    if (Vp_eff_quiescent > 0.0) {
        Ip_quiescent = pow(Vp_eff_quiescent, Ex) / KG1;
    }

    // --- 2. Calculate Dynamic (AC Signal) State ---
    // Total grid-cathode voltage = AC input + DC bias
    double Vgk_total = x + bias;

    // Emulate grid conduction (clipping) when grid goes positive
    if (Vgk_total > 0.0) {
        // This models the grid acting like a diode and clamping the signal
        Vgk_total = 1.1 * std::tanh(Vgk_total / 1.1);
    }

    // Effective grid voltage with signal
    double Vp_eff_total = (V_supply / Mu) + Vgk_total;
    // Plate current with signal
    double Ip_total = 0.0;
    if (Vp_eff_total > 0.0) {
        Ip_total = pow(Vp_eff_total, Ex) / KG1;
    }

    // --- 3. Calculate Output Voltage ---
    // The AC output voltage is the change in plate current multiplied by the load resistor.
    // The negative sign represents the inverting nature of a common-cathode stage.
    double V_out_ac = -(Ip_total - Ip_quiescent) * R_L;

    return V_out_ac;
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
    bias = -1.8; V_supply = 405.0; R_L = 150e3;
    reset();
}

void V1AStage::reset() {
    cathodeBypassFilter.reset();
    outputCouplingFilter.reset();
    interStageLPF.reset();
}

double V1AStage::process(double in) {
    double out = cathodeBypassFilter.process(in);
    out = tubeSaturateVirtualVoltage(out, bias, V_supply, R_L);
    out = outputCouplingFilter.process(out);
    out = interStageLPF.process(out);
    return out;
}

// --- V1B ---
void V1BStage::prepare(double sampleRate) {
    make_hpf(cathodeBypassFilter, sampleRate, 4.8);
    make_hpf(outputCouplingFilter, sampleRate, 0.4);
    make_lpf(interStageLPF, sampleRate, 2900.0); // Initialize the LPF
    bias = -1.8; V_supply = 405.0; R_L = 100e3;
    reset();
}

void V1BStage::reset() {
    cathodeBypassFilter.reset();
    outputCouplingFilter.reset();
    interStageLPF.reset();
}

double V1BStage::process(double in) {
    double out = cathodeBypassFilter.process(in);
    out = tubeSaturateVirtualVoltage(out, bias, V_supply, R_L);
    out = outputCouplingFilter.process(out);
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
    // V3B Parameters
    v3b_bias = -2.0; v3b_V_supply = 410.0; v3b_R_L = 82e3;
    // V4A Parameters
    v4a_bias = -2.2; v4a_V_supply = 410.0; v4a_R_L = 274e3;
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
    drive = (g * g);
}

double V3BV4AStage::process(double in) {
    // Attenuate signal based on drive knob, then feed to V3B
    double out = v3b_inputFilter.process(in * drive);
    out = v3b_cathodeBypass.process(out);
    // V3B now adds gain and distortion
    out = tubeSaturateVirtualVoltage(out, v3b_bias, v3b_V_supply, v3b_R_L);

    // Signal passes through coupling network to V4A
    out = interStageHPF.process(out);
    out = fizzFilter.process(out);

    // V4A gets slammed by the high-level signal from V3B.
    out = v4a_cathodeBypass.process(out);
    out = tubeSaturateVirtualVoltage(out, v4a_bias, v4a_V_supply, v4a_R_L);
    out = v4a_outputCoupling.process(out);

    return out;
}

// --- V2A ---
void V2AStage::prepare(double sampleRate) {
    make_hpf(inputCouplingFilter, sampleRate, 72.0);
    make_hpf(cathodeBypassFilter, sampleRate, 10.3);
    make_hpf(outputCouplingFilter, sampleRate, 3.4);
    bias = -1.5; V_supply = 410.0; R_L = 120e3;
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
    out = tubeSaturateVirtualVoltage(out, bias, V_supply, R_L);
    out = outputCouplingFilter.process(out);
    // Master volume now acts as a simple attenuator on the virtual voltage
    return out * masterVol;
}

// --- V2B ---
void V2BStage::prepare(double sampleRate) {
    make_hpf(outputCouplingFilter, sampleRate, 144.0);
    bias = -1.5; V_supply = 410.0; R_L = 100e3;
    reset();
}

void V2BStage::reset() {
    outputCouplingFilter.reset();
}

double V2BStage::process(double in) {
    double out = tubeSaturateVirtualVoltage(in, bias, V_supply, R_L);
    out = outputCouplingFilter.process(out);
    return out;
}

double Preamp::processSample(double in) {
    // Process through the chain
    double sample = v1a.process(in);
    measureMinMax(sample);

    sample = toneStack.process(sample);

#if 1
    double lead_path = v3b_v4a.process(sample);

    double rhythm_path = v1b.process(sample);

    // Mix the paths.
    sample = lead_path * 0.8 + rhythm_path * 0.2;

    sample = v2a.process(sample);
    sample = v2b.process(sample);

    sample /= 5000.0;
#else
    sample /= 30.0;
#endif

    return sample;
}
