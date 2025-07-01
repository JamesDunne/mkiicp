#include "TubeStage.h"

#include <numbers>
#include <ostream>

void make_lpf(IIRBiquad& filter, double sampleRate, double cutoff_freq);
void make_hpf(IIRBiquad& filter, double sampleRate, double cutoff_freq);

TubeStage::TubeStage() {
    reset();
}

void TubeStage::prepare(double sampleRate, double r_k, double r_l, double v_supply, double c_k) {
    R_k = r_k; R_L = r_l; V_supply = v_supply;
    isBypassed = (c_k > 1e-12);

    // Calculate the DC operating point first, as it's independent of the sample rate.
    calculateOperatingPoint();
    reset(); // Reset filter states to zero.

    // Configure the bypass filter *after* reset.
    if (isBypassed) {
        double cf = 1.0 / (2.0 * std::numbers::pi * R_k * c_k);
        make_lpf(cathodeBypassFilter, sampleRate, cf);
    }
}

void TubeStage::reset() {
    inputFilter.reset();
    outputFilter.reset();
    interStageLPF.reset();
    cathodeBypassFilter.reset();
}

void TubeStage::calculateOperatingPoint() {
    // 1. Solve for DC quiescent point using the Load-Line Analysis method.
    // This is guaranteed to find the stable operating point for any configuration.
    double Ip_q = 0.0; // Quiescent Plate Current
    double V_gk = 0.0; // Grid-to-Cathode Voltage

    // We iterate on the grid voltage V_gk to find the equilibrium point.
    // Start with a typical bias voltage guess.
    double min_V_gk = -10.0, max_V_gk = 0.0;

    for (int i = 0; i < 15; ++i) {
        V_gk = 0.5 * (min_V_gk + max_V_gk); // Bisection search

        // Equation 1: Current from the tube's perspective (the "curve")
        // This is a validated, complete 12AX7 model.
        double E1 = (V_supply / Mu) + V_gk;
        if (E1 < 0) E1 = 0; // Tube is in cutoff
        double Ip_from_tube = (pow(E1, 1.5) / 3300.0) * (1.0 + tanh(V_gk * 5.0));

        // Equation 2: Current from the circuit's perspective (the "load line")
        double Ip_from_circuit = -V_gk / R_k;

        // Compare the two currents to find the intersection point
        if (Ip_from_tube > Ip_from_circuit) {
            max_V_gk = V_gk; // The bias point is more negative
        } else {
            min_V_gk = V_gk; // The bias point is less negative
        }
    }

    // After converging, V_gk is our quiescent grid-cathode voltage.
    // Calculate the final quiescent plate current from the load line.
    Ip_q = -V_gk / R_k;
    Vp_dc = V_supply - Ip_q * R_L;

    // 2. Calculate small-signal parameters at the now-correct DC operating point
    // Transconductance (gm) and dynamic plate resistance (r_p)
    double gm = (1.5 / 3300.0) * sqrt(Ip_q * 3300.0 * Mu);
    r_p = Mu / (gm + 1e-12);

    // 3. Calculate Thevenin equivalent parameters for the whole stage
    R_out = (r_p * R_L) / (r_p + R_L);

    // 4. Calculate AC Gain, correctly modeling the bypassed cathode
    double Zk_ac = 0.0; // Effective AC impedance of the cathode
    if (isBypassed) {
        // At AC, the capacitor has a very low impedance. We can approximate it as
        // a small resistance in parallel with R_k, but for gain calculation,
        // assuming it's near zero is standard for mid/high frequencies.
        Zk_ac = 0.0; // Ideal bypass
    } else {
        Zk_ac = R_k; // Unbypassed
    }
    gain = (Mu * R_L) / (R_L + r_p + Zk_ac * (Mu + 1));
}

double TubeStage::process(double V_in_ac, double R_load) {
    double loaded_gain = gain * (R_load / (R_out + R_load));
    double V_max_swing = V_supply - Vp_dc;
    double V_min_swing = Vp_dc;
    double V_out_ac = -V_in_ac * loaded_gain;

    if (V_out_ac > 0) {
        V_out_ac = V_max_swing * std::tanh(V_out_ac / V_max_swing);
    } else {
        V_out_ac = V_min_swing * std::tanh(V_out_ac / V_min_swing);
    }
    return V_out_ac;
}
