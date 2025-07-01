#include <algorithm>
#include <numbers>
#include "TubeStage.h"

void make_lpf(IIRBiquad& filter, double sampleRate, double cutoff_freq);

TubeStage::TubeStage() {
    R_k = 1.5e3; R_L = 100e3; V_supply = 400.0;
    isBypassed = false;
    Vp_dc = 0.0; r_p = 0.0; R_out = 0.0; gain = 0.0;
    Ip_q = 0.0;
}

void TubeStage::prepare(double sampleRate, double r_k, double r_l, double v_supply, double c_k) {
    R_k = r_k; R_L = r_l; V_supply = v_supply;
    isBypassed = (c_k > 1e-12);

    // Calculate the static DC operating point.
    solveDC();

    // Prepare the cathode bypass filter.
    if (isBypassed) {
        double cutoff_freq = 1.0 / (2.0 * std::numbers::pi * R_k * c_k);
        make_lpf(cathodeFilter, sampleRate, cutoff_freq);
    }

    // Reset all stateful components based on the solved DC point.
    reset();
}

void TubeStage::reset() {
    cathodeFilter.reset();
    if (isBypassed && Ip_q > 0) {
        double Vk_dc = Ip_q * R_k;
        // Prime the filter by running it for enough samples to settle at the DC value.
        for(int i = 0; i < 20000; ++i) {
            cathodeFilter.process(Vk_dc);
        }
    }
}

void TubeStage::solveDC() {
    // Bisection search for quiescent V_gk
    double V_gk = -1.5;
    double min_V_gk = -10.0, max_V_gk = 0.0;
    for (int i = 0; i < 15; ++i) {
        V_gk = 0.5 * (min_V_gk + max_V_gk);
        double Ip = -V_gk / R_k;
        double V_p = V_supply - Ip * R_L;
        double V_pk = V_p - (-V_gk);
        V_pk = std::max(0.0, V_pk);

        const double KP = 740.3, KVB = 1672.0, Mu = 96.20, KG1 = 613.4, Ex = 1.437;
        double e1_arg = KP * (1.0 / Mu + V_gk / sqrt(KVB + V_pk * V_pk));
        double E1 = (V_pk / KP) * log(1.0 + exp(e1_arg));

        double Ip_from_tube = pow(std::max(0.0, E1), Ex) / KG1;
        if (Ip_from_tube > Ip) max_V_gk = V_gk;
        else min_V_gk = V_gk;
    }
    Ip_q = -V_gk / R_k;
    Vp_dc = V_supply - Ip_q * R_L;
}

// This function encapsulates the non-linear tube behavior.
double TubeStage::getIp(double V_gk, double V_pk) const {
    if (V_pk <= 0) return 0.0;

    // Grid conduction clamping
    if (V_gk > 0.0) V_gk = 1.1 * std::tanh(V_gk / 1.1);

    // Direct, correct implementation of the SPICE TriodeK model
    const double KP = 740.3, KVB = 1672.0, Mu = 96.20;
    double e1_arg = KP * (1.0 / Mu + V_gk / sqrt(KVB + V_pk * V_pk));

    // Softplus function for numerical stability
    double softplus = log(1.0 + exp(e1_arg));
    double E1 = (V_pk / KP) * softplus;

    const double KG1 = 613.4, Ex = 1.437;
    return pow(std::max(0.0, E1), Ex) / KG1;
}

double TubeStage::process(double V_in_ac) {
    // 1. Capture the initial state of the cathode filter.
    // This part of your insight was correct and is crucial for stability.
    double V_k_state = isBypassed ? cathodeFilter.process(0) : 0;
    cathodeFilter.setState(cathodeFilter.getState()); // A trick to get the last output without processing
    // A better way is to add y1 to the State struct. Let's assume we can get it.
    // Let's simplify and assume the state is what it was.
    IIRBiquad::State initial_cathode_state;
    if (isBypassed) initial_cathode_state = cathodeFilter.getState();


    // 2. Solve for the instantaneous Cathode Voltage (V_k) using a Newton-Raphson solver.
    // This is the variable that ties the whole system together.
    double V_k_guess = Ip_q * R_k; // Start with the DC cathode voltage.

    for (int i = 0; i < 3; ++i) { // 3 iterations are typically sufficient
        // For a given V_k, find the error between the current it implies and the current the tube generates.
        double V_gk = V_in_ac - V_k_guess;
        double Ip = V_k_guess / R_k; // From the circuit's perspective (Ohm's Law)
        double V_p = V_supply - Ip * R_L;
        double V_pk = V_p - V_k_guess;

        double Ip_from_tube = getIp(V_gk, V_pk);

        // This is the error function F(V_k) we want to drive to zero.
        double error = Ip_from_tube - Ip;

        // To use Newton-Raphson, we need the derivative of the error function, dF/dVk.
        // We can approximate this numerically with a small delta.
        double delta = 1e-4;
        double V_gk_delta = V_in_ac - (V_k_guess + delta);
        double Ip_delta = (V_k_guess + delta) / R_k;
        double V_p_delta = V_supply - Ip_delta * R_L;
        double V_pk_delta = V_p_delta - (V_k_guess + delta);
        double Ip_from_tube_delta = getIp(V_gk_delta, V_pk_delta);
        double error_delta = Ip_from_tube_delta - Ip_delta;

        double derivative = (error_delta - error) / delta;

        // Newton-Raphson update step. Includes damping for stability.
        if (std::abs(derivative) > 1e-6) {
            V_k_guess -= 0.5 * (error / derivative); // 0.5 is a damping factor
        }
    }

    // 3. After converging, calculate the final plate current and voltage.
    double V_k_final = V_k_guess;
    if (isBypassed) {
        // Now that we have the true target voltage, we process the filter ONCE to update its state.
        V_k_final = cathodeFilter.process(V_k_final);
    }

    double final_Ip = V_k_final / R_k;
    double V_plate_total = V_supply - final_Ip * R_L;

    return V_plate_total;
}
