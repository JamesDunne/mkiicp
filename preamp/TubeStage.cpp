#include "TubeStage.h"
#include <cmath>
#include <numbers>
#include <algorithm>

// Assumes make_lpf utility function is available
void make_lpf(IIRBiquad& filter, double sampleRate, double cutoff_freq);

TubeStage::TubeStage() { reset(); }

void TubeStage::prepare(double sampleRate, double r_k, double r_l, double v_supply, double c_k) {
    R_k = r_k;
    R_L = r_l;
    V_supply = v_supply;
    isBypassed = (c_k > 1e-12);

    if (isBypassed) {
        double cutoff_freq = 1.0 / (2.0 * std::numbers::pi * R_k * c_k);
        make_lpf(cathodeFilter, sampleRate, cutoff_freq);
    }

    solveDC(); // Calculate the static DC operating point
    reset();   // Reset all stateful components
}

void TubeStage::reset() {
    // Prime the LPF state with the DC cathode voltage to prevent startup thumps.
    // This must be done after solveDC() has calculated Ip_q.
    if (isBypassed && Ip_q > 0) {
        double Vk_dc = Ip_q * R_k;
        cathodeFilter.reset();

        // To prime a Direct Form 1 IIR filter, we set its state variables (x1, y1)
        // such that for a constant input of Vk_dc, the output is also Vk_dc.
        // This is achieved by setting the past states equal to the DC value.
        // We'll add a simple `prime` method to IIRBiquad for this.
        // For now, let's process it for a long time to settle it.
        for(int i=0; i<10000; ++i) {
            cathodeFilter.process(Vk_dc);
        }
    }
}

void TubeStage::solveDC() {
    double Ip_q = 0.0008; // Start with a reasonable guess.

    // This bisection search is guaranteed to converge on the correct Vgk.
    double min_V_gk = -10.0, max_V_gk = 0.0;
    double V_gk = -1.5; // A typical starting point for Vgk

    for (int i = 0; i < 15; ++i) {
        // From V_gk, find the current the tube "wants" to draw.
        // This must use the full plate voltage equation.
        double V_p = V_supply - (-V_gk / R_k) * R_L;
        double V_pk = V_p + V_gk; // Vp - Vk = Vp - (-Vgk)
        V_pk = std::max(0.0, V_pk);

        const double KP = 740.3, KVB = 1672.0;
        double e1_arg = KP * (1.0 / Mu + V_gk / sqrt(KVB + V_pk * V_pk));
        double E1 = (e1_arg > 30.0) ? (V_pk / KP * e1_arg) : (V_pk / KP * log(1.0 + exp(e1_arg)));

        const double KG1 = 613.4, Ex = 1.437;
        double Ip_from_tube = pow(std::max(0.0, E1), Ex) / KG1;

        // From V_gk, find the current the circuit "allows" via the load line.
        double Ip_from_circuit = -V_gk / R_k;

        // Adjust the search range.
        if (Ip_from_tube > Ip_from_circuit) {
            max_V_gk = V_gk;
        } else {
            min_V_gk = V_gk;
        }
        V_gk = 0.5 * (min_V_gk + max_V_gk);
    }

    // After converging, V_gk is our quiescent grid-cathode voltage.
    this->Ip_q = -V_gk / R_k;
    this->Vp_dc = V_supply - this->Ip_q * R_L;
}

double TubeStage::process(double V_in_ac) {
    // For a zero input signal, the tube should be at its DC resting point.
    // The iterative solver must converge to Ip_q.
    if (std::abs(V_in_ac) < 1e-9) {
        // If the input is silent, we don't need to solve.
        // We also must ensure the cathode filter state remains at DC.
        if (isBypassed) {
            cathodeFilter.process(Ip_q * R_k);
        }
        return Vp_dc;
    }

    double current_Ip = Ip_q; // Start solver at the known DC quiescent current.

    for (int i = 0; i < 5; ++i) {
        // Calculate the raw, unfiltered voltage that *would* appear at the cathode.
        double Vk_unfiltered = current_Ip * R_k;

        // Determine the actual cathode voltage by applying the bypass LPF.
        double V_k = isBypassed ? cathodeFilter.process(Vk_unfiltered) : Vk_unfiltered;

        // Calculate inter-electrode voltages using the correct, filtered V_k.
        double V_pk = V_supply - (current_Ip * R_L) - V_k;
        double V_gk = V_in_ac - V_k;

        V_pk = std::max(0.0, V_pk);
        if (V_gk > 0.0) { V_gk = 1.1 * std::tanh(V_gk / 1.1); }

        // Calculate target plate current using the full, correct SPICE formula.
        const double KP = 740.3, KVB = 1672.0;
        double e1_arg = KP * (1.0 / Mu + V_gk / sqrt(KVB + V_pk * V_pk));
        double E1 = (e1_arg > 30.0) ? (V_pk / KP * e1_arg) : (V_pk / KP * log(1.0 + exp(e1_arg)));

        const double KG1 = 613.4, Ex = 1.437;
        double target_Ip = pow(std::max(0.0, E1), Ex) / KG1;

        // Damped Update
        current_Ip += (target_Ip - current_Ip) * 0.5;
    }

    // Return the total, instantaneous plate voltage.
    return V_supply - current_Ip * R_L;
}
