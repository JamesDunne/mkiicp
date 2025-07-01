#include "TubeStage.h"

#include <numbers>
#include <ostream>

void make_lpf(IIRBiquad& filter, double sampleRate, double cutoff_freq);
void make_hpf(IIRBiquad& filter, double sampleRate, double cutoff_freq);

TubeStage::TubeStage() {
    reset();
}

void TubeStage::prepare(double sampleRate, double r_k, double r_l, double v_supply, double c_k) {
    R_k = r_k;
    R_L = r_l;
    V_supply = v_supply;
    C_k = c_k;
    isBypassed = (C_k > 1e-12); // Is there a bypass cap?

    if (isBypassed) {
        // Calculate cutoff freq of the bypass RC network
        double cutoff_freq = 1.0 / (2.0 * std::numbers::pi * R_k * C_k);
        make_lpf(cathodeBypassFilter, sampleRate, cutoff_freq);
    }

    solveDC();
    reset();
}

void TubeStage::reset() {
    inputFilter.reset();
    outputFilter.reset();
    interStageLPF.reset();
    cathodeBypassFilter.reset();
    if (Ip_q > 0) {
        // Prime the LPF state with the DC cathode voltage to prevent startup thumps
        cathodeBypassFilter.process(Ip_q * R_k);
        cathodeBypassFilter.process(Ip_q * R_k);
    }
}

// Calculates the quiescent (DC) plate current using an iterative solver.
void TubeStage::solveDC() {
    double current_Ip = 0.001; // Start with a reasonable guess (1mA)

    for (int i = 0; i < 10; ++i) {
        // Direct, stable formulation for inter-electrode voltages
        double V_pk = V_supply - current_Ip * (R_L + R_k);
        double V_gk = -current_Ip * R_k;

        V_pk = std::max(0.0, V_pk); // Ensure physicality

        // Calculate target current using the SPICE model
        double e1_arg = KP * (1.0 / Mu + V_gk / sqrt(KVB + V_pk * V_pk));
        double E1 = (e1_arg > 30.0) ? (V_pk / KP * e1_arg) : (V_pk / KP * log(1.0 + exp(e1_arg)));
        double target_Ip = pow(std::max(0.0, E1), Ex) / KG1;

        // Damped update for smooth convergence
        current_Ip += (target_Ip - current_Ip) * 0.5;
    }

    Ip_q = current_Ip;
    Vp_dc = V_supply - Ip_q * R_L; // Store the DC plate voltage for reference
}

// Processes one sample using the stateful SPICE model
double TubeStage::process(double V_in_ac) {
    double current_Ip = Ip_q;

    for (int i = 0; i < 5; ++i) {
        // 1. Calculate the raw, unfiltered voltage that would appear at the cathode
        double Vk_unfiltered = current_Ip * R_k;

        // 2. Determine the actual cathode voltage by applying the bypass filter
        double V_k;
        if (isBypassed) {
            // The LPF smooths the cathode voltage, shorting AC components to ground
            V_k = cathodeBypassFilter.process(Vk_unfiltered);
        } else {
            // If unbypassed, the full voltage is present
            V_k = Vk_unfiltered;
        }

        // 3. Calculate inter-electrode voltages based on this correct V_k
        double V_pk = V_supply - current_Ip * R_L - V_k;
        double V_gk = V_in_ac - V_k;

        V_pk = std::max(0.0, V_pk);
        if (V_gk > 0.0) { V_gk = 1.1 * std::tanh(V_gk / 1.1); }

        // 4. Calculate target plate current
        double e1_arg = KP * (1.0 / Mu + V_gk / sqrt(KVB + V_pk * V_pk));
        double E1 = (e1_arg > 30.0) ? (V_pk / KP * e1_arg) : (V_pk / KP * log(1.0 + exp(e1_arg)));
        double target_Ip = pow(std::max(0.0, E1), Ex) / KG1;

        // 5. Damped Update
        current_Ip += (target_Ip - current_Ip) * 0.5;
    }

    // 6. Numerically Stable Output Calculation
    double I_ac = current_Ip - Ip_q;
    double V_out_ac = -I_ac * R_L;

    return V_out_ac;
}
