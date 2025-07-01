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

    solveDC(); // Calculate the static DC operating point
    reset();   // Reset all stateful components

    if (isBypassed) {
        double cutoff_freq = 1.0 / (2.0 * std::numbers::pi * R_k * c_k);
        make_lpf(cathodeFilter, sampleRate, cutoff_freq);
    }
}

void TubeStage::reset() {
    cathodeFilter.reset();
    // Prime the LPF state with the DC cathode voltage to prevent startup thumps.
    // This must be done after solveDC() has calculated Ip_q.
    if (isBypassed && Ip_q > 0) {
        double Vk_dc = Ip_q * R_k;
        // Process a few times to settle the filter's state
        for(int i=0; i<100; ++i) {
            cathodeFilter.process(Vk_dc);
        }
    }
}

void TubeStage::solveDC() {
    // ... same robust load-line solver from the previous version ...
    double current_Ip = 0.0008;
    for (int i = 0; i < 15; ++i) {
        double V_gk = -current_Ip * R_k;
        double E1 = (V_supply / Mu) + V_gk;
        if (E1 <= 0) {
            current_Ip = 0;
            break;
        }
        double target_Ip = pow(E1, 1.5) / 3300.0;
        current_Ip += (target_Ip - current_Ip) * 0.5;
    }
    Ip_q = current_Ip; // Store the result in the member variable
    Vp_dc = V_supply - Ip_q * R_L;
}

double TubeStage::process(double V_in_ac) {
    // ... same definitive process method from the previous version ...
    double current_Ip = Ip_q;
    for (int i = 0; i < 5; ++i) {
        double Vk_unfiltered = current_Ip * R_k;
        double V_k = isBypassed ? cathodeFilter.process(Vk_unfiltered) : Vk_unfiltered;
        double V_pk = V_supply - (current_Ip * R_L) - V_k;
        double V_gk = V_in_ac - V_k;
        V_pk = std::max(0.0, V_pk);
        if (V_gk > 0.0) { V_gk = 1.1 * std::tanh(V_gk / 1.1); }
        double E1 = (V_pk / Mu) + V_gk;
        if (E1 < 0) E1 = 0;
        double target_Ip = pow(E1, 1.5) * (1.0 / 3300.0);
        current_Ip += (target_Ip - current_Ip) * 0.5;
    }
    double I_ac = current_Ip - Ip_q;
    double V_out_ac = -I_ac * R_L;
    return V_out_ac;
}
