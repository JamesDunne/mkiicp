#include "Preamp.h"
#include <numbers>

constexpr double PI = std::numbers::pi;

/**
 * @brief Final, correctly normalized 12AX7 waveshaper.
 *
 * This function is designed to take a signal in a nominal [-1, 1] range
 * and return a distorted signal also in a nominal [-1, 1] range.
 * It is DC-free and all gain is handled externally.
 *
 * @param x The input signal, [-1, 1].
 * @param drive Controls how much of the input signal is used, pushing it into non-linearity.
 * @param bias The cathode bias point in Volts. Typical: -1.5V.
 * @return The distorted, normalized output signal.
 */
inline double tubeSaturate12AX7(double x, double drive, double bias = -1.5) {
    // --- Parameters ---
    const double Mu = 96.20;
    const double Ex = 1.437;
    const double Vp = 270.0;

    // 1. Calculate Grid Voltage (Vgk). 'drive' scales the input signal.
    double Vgk = x * drive + bias;

    // 2. Grid Conduction.
    if (Vgk > 0.0) {
        Vgk = 0.5 * std::tanh(Vgk * 2.0);
    }

    // 3. Plate Current Calculation (without artificial scaling)
    double Vp_eff = (Vp / Mu) + Vgk;
    double plate_current = (Vp_eff > 0.0) ? std::pow(Vp_eff, Ex) : 0.0;

    // 4. Calculate Quiescent Current to make the output AC-coupled (DC-free)
    double Vp_eff_quiescent = (Vp / Mu) + bias;
    double quiescent_current = (Vp_eff_quiescent > 0.0) ? std::pow(Vp_eff_quiescent, Ex) : 0.0;

    // 5. Calculate output and normalize
    // The difference is the AC part of the signal.
    double output = -(plate_current - quiescent_current);

    // 6. Normalize the output. This scaling factor is tuned to ensure that
    // for an input of x=1, drive=1, the output is close to 1.0. This prevents
    // the function itself from adding unwanted gain.
    return output * 0.4;
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
    make_hpf(outputCouplingFilter, sampleRate, 31.0);
    // Let's give V1A a more realistic "hot" gain before the tone stack.
    gain = 15.0;
    reset();
}

void V1AStage::reset() {
    cathodeBypassFilter.reset();
    outputCouplingFilter.reset();
}

double V1AStage::process(double in) {
    double out = cathodeBypassFilter.process(in);
    out = tubeSaturate12AX7(out, 1.0, -1.8);
    out = outputCouplingFilter.process(out);
    // Output a hot signal, as it would be at the tube's plate.
    return out * gain;
}

// --- Tone Stack ---
void ToneStack::prepare(double sampleRate) {
    this->sampleRate = sampleRate;
    calculateCoefficients();
    reset();
}

void ToneStack::reset() {
    m_s_c56 = m_s_c4 = m_s_c3 = m_s_c13b = 0.0;
}

void ToneStack::setParams(double treble, double mid, double bass, double volume) {
    // Clamp input parameters to the [0, 1] range
    treble_p = std::clamp(treble, 0.0, 1.0);
    bass_p = std::clamp(bass, 0.0, 1.0);
    mid_p = std::clamp(mid, 0.0, 1.0);
    volume_p = std::clamp(volume, 0.0, 1.0);
    // analog tapers:
    treble_p *= treble_p;
    bass_p *= bass_p;
    mid_p *= mid_p;
    volume_p *= volume_p;
    calculateCoefficients();
}

void ToneStack::calculateCoefficients() {
    // --- Stabilize Potentiometer Travel ---
    // To prevent the MNA matrix from becoming singular, we must prevent any
    // resistance from becoming a perfect zero. We map the user's [0,1] range
    // to a slightly smaller internal range, e.g., [0.001, 0.999].
    const double travel_min = 0.001;
    const double travel_max = 1.0 - travel_min;
    double treble_t = treble_p * travel_max + travel_min;
    double volume_t = volume_p * travel_max + travel_min;

    // --- Define Component Values and Epsilon for Stability ---
    const double epsilon = 1e-12;

    const double C56 = 1.0e-9;
    const double C4 = 1.0e-7;
    const double C3 = 4.7e-8;
    const double C13B = 1.8e-10;

    // --- Calculate Component Conductances (G = 1/R) ---
    m_g_r5 = 1.0 / 100000.0;
    const double g_r5a = 1.0 / 100000.0;

    double r_treble_a = 250000.0 * (1.0 - treble_t);
    double r_treble_c = 250000.0 * treble_t;
    double g_treble_a = 1.0 / (r_treble_a); // No epsilon needed due to travel clamp
    double g_treble_c = 1.0 / (r_treble_c);

    double r_bass_a = 250000.0 * bass_p;
    double g_bass_a = 1.0 / (r_bass_a + epsilon); // Epsilon needed as this can be 0

    double r_mid_a = 10000.0 * mid_p;
    double g_mid_a = 1.0 / (r_mid_a + epsilon); // Epsilon needed

    double r_vol_a = 1000000.0 * (1.0 - volume_t);
    double r_vol_c = 1000000.0 * volume_t;
    double g_vol_a = 1.0 / (r_vol_a);
    double g_vol_c = 1.0 / (r_vol_c);

    // --- Update Capacitor Conductances for Trapezoidal Integration ---
    m_g_c56 = 2.0 * C56 * sampleRate;
    m_g_c4 = 2.0 * C4 * sampleRate;
    m_g_c3 = 2.0 * C3 * sampleRate;
    m_g_c13b = 2.0 * C13B * sampleRate;

    // --- Build the 7x7 MNA System Matrix 'A' ---
    double A[7][7] = {0};

    // KCL at N005
    A[0][0] = g_treble_a + m_g_c56; A[0][1] = -g_treble_a;
    // KCL at N007
    A[1][0] = -g_treble_a; A[1][1] = g_treble_a + g_treble_c + g_r5a; A[1][3] = -g_treble_c; A[1][5] = -g_r5a;
    // KCL at N015
    A[2][2] = m_g_r5 + m_g_c4 + m_g_c3; A[2][3] = -m_g_c4; A[2][4] = -m_g_c3;
    // KCL at N016
    A[3][1] = -g_treble_c; A[3][2] = -m_g_c4; A[3][3] = g_treble_c + g_bass_a + m_g_c4; A[3][4] = -g_bass_a;
    // KCL at N024
    A[4][2] = -m_g_c3; A[4][3] = -g_bass_a; A[4][4] = g_bass_a + g_mid_a + m_g_c3;
    // KCL at N008
    A[5][1] = -g_r5a; A[5][5] = g_r5a + g_vol_a + m_g_c13b; A[5][6] = -g_vol_a - m_g_c13b;
    // KCL at N020
    A[6][5] = -g_vol_a - m_g_c13b; A[6][6] = g_vol_a + g_vol_c + m_g_c13b;

    // Pre-calculate the inverse of the system matrix
    invertMatrix(A, m_A_inv);
}

double ToneStack::process(double vin) {
    // Build the 'b' vector for Ax = b
    double b[7];
    b[0] = m_s_c56 + vin * m_g_c56;
    b[1] = 0.0;
    b[2] = vin * m_g_r5 + m_s_c4 + m_s_c3;
    b[3] = -m_s_c4;
    b[4] = -m_s_c3;
    b[5] = m_s_c13b;
    b[6] = -m_s_c13b;

    // Solve for the node voltages: v = A_inv * b
    double v[7];
    for (int i = 0; i < 7; ++i) {
        v[i] = 0.0;
        for (int j = 0; j < 7; ++j) {
            v[i] += m_A_inv[i][j] * b[j];
        }
    }

    // The output is the voltage at node N020 (index 6)
    double output = v[6];

    // Update capacitor state variables for the next time step
    m_s_c56  = 2.0 * m_g_c56  * (v[0] - vin) - m_s_c56;
    m_s_c4   = 2.0 * m_g_c4   * (v[2] - v[3]) - m_s_c4;
    m_s_c3   = 2.0 * m_g_c3   * (v[2] - v[4]) - m_s_c3;
    m_s_c13b = 2.0 * m_g_c13b * (v[5] - v[6]) - m_s_c13b;

    return output;
}

// --- V1B ---
void V1BStage::prepare(double sampleRate) {
    make_hpf(cathodeBypassFilter, sampleRate, 4.8);
    make_hpf(outputCouplingFilter, sampleRate, 0.4);
    // This stage expects the post-tonestack signal, so its gain can be modest.
    gain = 4.0;
    reset();
}

void V1BStage::reset() {
    cathodeBypassFilter.reset();
    outputCouplingFilter.reset();
}

double V1BStage::process(double in) {
    double out = cathodeBypassFilter.process(in);
    out = tubeSaturate12AX7(out, 1.0, -1.8);
    out = outputCouplingFilter.process(out);
    return out * gain;
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
    // Maps the [0,1] knob to a powerful gain factor [1x to 51x].
    // This is the core of the lead drive.
    drive = 1.0 + 50.0 * (g * g);
}

double V3BV4AStage::process(double in) {
    // V3B acts as a clean booster whose output level is controlled by the 'drive' knob.
    double out = v3b_inputFilter.process(in);
    out = v3b_cathodeBypass.process(out);
    out *= drive; // Apply the main lead gain here.

    // V4A gets slammed by the high-level signal from V3B.
    // This is where the heavy clipping happens.
    out = v4a_cathodeBypass.process(out);
    out = tubeSaturate12AX7(out, 1.0, -2.0); // Drive is 1.0 as the signal is already huge.
    out = v4a_outputCoupling.process(out);

    return out;
}

// --- V2A ---
void V2AStage::prepare(double sampleRate) {
    make_hpf(cathodeBypassFilter, sampleRate, 10.3);
    make_hpf(outputCouplingFilter, sampleRate, 3.4);
    gain = 1.5; // This stage is a recovery/driver stage.
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
    out = tubeSaturate12AX7(out, 1.0, -1.5);
    out *= gain;
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
    double out = tubeSaturate12AX7(in, 1.0, -1.5);
    out = outputCouplingFilter.process(out);
    return out * gain;
}

double Preamp::processSample(double in) {
    // Process through the chain
    double sample = v1a.process(in);

    sample = toneStack.process(sample);

    // 3. CRUCIAL: Apply Makeup Gain to compensate for tone stack insertion loss.
    // Your numbers showed a ~20x drop. We'll compensate for that here.
    // This brings the signal back to a healthy level for the next stages.
    sample *= 20.0;

    double lead_path = v3b_v4a.process(sample);

    double rhythm_path = v1b.process(sample);

    // Mix the paths.
    sample = lead_path * 0.8 + rhythm_path * 0.2;

    sample = v2a.process(sample);
    sample = v2b.process(sample);

    measureMinMax(sample);

    return sample;
}
