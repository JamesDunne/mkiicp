#include "TubeStage.h"
#include <ostream>

TubeStage::TubeStage() {
    reset();
}

void TubeStage::prepare(double sampleRate, double r_k, double r_l, double v_supply) {
    R_k = r_k;
    R_L = r_l;
    V_supply = v_supply;
    solveDC(); // Calculate the operating point once
    reset();
}

void TubeStage::reset() {
    // Reset filters
    inputFilter.reset();
    outputFilter.reset();
    cathodeBypass.reset();
    interStageLPF.reset();
    // Reset state to the solved DC point
    last_Ip = Ip_q;
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
    // Start our guess for the current sample's plate current with the quiescent one.
    double current_Ip = Ip_q;

    // Iterative solver to find the stable plate current for this sample
    for (int i = 0; i < 5; ++i) {
        // 1. Direct and stable formulation for inter-electrode voltages
        double V_pk = V_supply - current_Ip * (R_L + R_k);
        double V_gk = V_in_ac - current_Ip * R_k;

        // 2. Physical Limits and Grid Conduction
        V_pk = std::max(0.0, V_pk);

        if (V_gk > 0.0) {
            V_gk = 1.1 * std::tanh(V_gk / 1.1);
        }

        // 3. Calculate the target plate current
        double e1_arg = KP * (1.0 / Mu + V_gk / sqrt(KVB + V_pk * V_pk));
        double E1 = (e1_arg > 30.0) ? (V_pk / KP * e1_arg) : (V_pk / KP * log(1.0 + exp(e1_arg)));
        double target_Ip = pow(std::max(0.0, E1), Ex) / KG1;

        // 4. Damped Update
        current_Ip += (target_Ip - current_Ip) * 0.5;
    }

    // 5. Update the state for the next sample.
    // NOTE: This was a bug in the previous version. We must use the *converged* current.
    // The `last_Ip` state is not needed inside the loop.
    // last_Ip = current_Ip; // This was wrong logic. The converged current is what matters.

    // --- 6. The Numerically Stable Output Calculation ---
    // Calculate the AC component of the current (Total - DC)
    double I_ac = current_Ip - Ip_q;

    // The AC output voltage is the result of the AC current flowing through the plate resistor.
    // The negative sign represents the inverting nature of the common-cathode stage.
    double V_out_ac = -I_ac * R_L;

    return V_out_ac;
}
