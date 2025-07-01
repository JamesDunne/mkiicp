
#include "ToneStack.h"
#include "IIRFilter.h"

// --- Tone Stack (Accurate MNA-based model) ---
void ToneStack::prepare(double sampleRate) {
    this->sampleRate = sampleRate;
    reset();
}

void ToneStack::reset() {
    filter.reset();
}

void ToneStack::setParams(double treble, double mid, double bass, double volume) {
    p_treble = treble;
    p_mid = mid;
    p_bass = bass;
    p_vol = volume;
    calculateCoefficients();

    // Volume pot is logarithmic (Audio Taper)
    volumeGain = p_vol * p_vol;
}

void ToneStack::calculateCoefficients() {
    // Component values from the SPICE netlist
    const double R5 = 100e3;
    const double C_bright = 1000e-12; // C5+C6
    const double C4 = 0.1e-6;
    const double C3 = 0.047e-6;

    // Pots as variable resistors
    const double R_treble_pot = 250e3;
    const double R_bass_resistor = 250e3 * (p_bass * p_bass);
    const double R_mid_resistor = 10e3 * (p_mid * p_mid);

    // Treble pot modeled as two separate resistors
    const double Ra_treble = R_treble_pot * (1.0 - p_treble);
    const double Rc_treble = R_treble_pot * p_treble;

    // --- Analog Coefficients from Symbolic MNA Solution H(s) = N(s)/D(s) ---
    // N(s) = n2*s^2 + n1*s + n0
    // D(s) = d3*s^3 + d2*s^2 + d1*s + d0

    // These formulas are derived offline using a tool like SymPy or Mathematica
    double n2 = C_bright * C3 * R5 * Ra_treble * R_mid_resistor;
    double n1 = C_bright * R5 * Ra_treble + C_bright * Ra_treble * R_mid_resistor + C3 * R5 * R_mid_resistor;
    double n0 = R5 + R_mid_resistor;

    double d3 = C_bright * C3 * C4 * R5 * Ra_treble * Rc_treble * R_bass_resistor; // Highest order term
    double d2 = C_bright * C3 * R5 * Ra_treble * R_bass_resistor + C_bright * C4 * R5 * Rc_treble * R_bass_resistor +
                C_bright * C4 * Ra_treble * Rc_treble * R_bass_resistor + C3 * C4 * R5 * Ra_treble * Rc_treble +
                C_bright * C3 * R5 * Ra_treble * R_mid_resistor + C_bright * C4 * Ra_treble * Rc_treble * R_mid_resistor +
                C3 * C4 * R5 * Rc_treble * R_mid_resistor;
    double d1 = C_bright * R5 * Ra_treble + C_bright * Ra_treble * R_bass_resistor + C4 * R5 * Rc_treble +
                C4 * Ra_treble * Rc_treble + C3 * R5 * R_bass_resistor + C_bright * Ra_treble * R_mid_resistor +
                C_bright * Rc_treble * R_mid_resistor + C4 * R5 * R_mid_resistor + C4 * Rc_treble * R_mid_resistor +
                C3 * R5 * R_mid_resistor + C3 * R_bass_resistor * R_mid_resistor;
    double d0 = R5 + Ra_treble + Rc_treble + R_mid_resistor;

    // --- Bilinear Transform: s = 2*fs*(1-z^-1)/(1+z^-1) ---
    double T = 1.0 / sampleRate;
    double T2 = T * T;
    double T3 = T2 * T;
    double C = 4.0 / T2;
    double D = 8.0 / T3;

    double d_a0 = d3*D + d2*C + d1*2/T + d0;

    double b[4], a[4];

    // Numerator (b) coefficients
    double n_a0 = n2*C + n1*2/T + n0; // Our numerator is only 2nd order
    b[0] = n_a0 / d_a0;
    b[1] = (2*n0 - 2*n2*C) / d_a0;
    b[2] = (n2*C - n1*2/T + n0) / d_a0;
    b[3] = 0.0; // It's a 2nd order numerator

    // Denominator (a) coefficients
    a[0] = 1.0;
    a[1] = (-d3*3*D - d2*C + d1*2/T + 3*d0) / d_a0;
    a[2] = (d3*3*D - d2*C - d1*2/T + 3*d0) / d_a0;
    a[3] = (-d3*D + d2*C - d1*2/T + d0) / d_a0;

    filter.setCoefficients(b, a);
}

double ToneStack::process(double in) {
    // Process through the accurate filter and apply volume control
    return filter.process(in) * volumeGain;
}
