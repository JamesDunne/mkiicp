#include "TubeStage.h"
#include <algorithm>

TubeStage::TubeStage() {
    R_k = 1.5e3; R_L = 100e3; V_supply = 400.0;
    isBypassed = false;
    Vp_dc = 0.0; r_p = 0.0; R_out = 0.0; gain = 0.0;
}

void TubeStage::prepare(double sampleRate, double r_k, double r_l, double v_supply, double c_k) {
    R_k = r_k; R_L = r_l; V_supply = v_supply;
    isBypassed = (c_k > 1e-12);

    calculateOperatingPoint();
}

void TubeStage::reset() {
    // This model is stateless, so reset does nothing.
}

// Calculates DC operating point using robust load-line analysis.
void TubeStage::calculateOperatingPoint() {
    constexpr double KP = 740.3, KVB = 1672.0, KG1 = 613.4, Ex = 1.437;
    double Ip_q = 0.0008;
    double V_gk = -1.5;
    double min_V_gk = -10.0, max_V_gk = 0.0;

    for (int i = 0; i < 15; ++i) {
        V_gk = 0.5 * (min_V_gk + max_V_gk);
        double V_p = V_supply - (-V_gk / R_k) * R_L;
        double V_pk = V_p + V_gk;
        V_pk = std::max(0.0, V_pk);

        double e1_arg = KP * (1.0 / Mu + V_gk / sqrt(KVB + V_pk * V_pk));
        double E1 = (e1_arg > 30.0) ? (V_pk/KP*e1_arg) : (V_pk/KP*log(1.0+exp(e1_arg)));
        double Ip_from_tube = pow(std::max(0.0, E1), Ex) / KG1;
        double Ip_from_circuit = -V_gk / R_k;

        if (Ip_from_tube > Ip_from_circuit) max_V_gk = V_gk;
        else min_V_gk = V_gk;
    }

    Ip_q = -V_gk / R_k;
    Vp_dc = V_supply - Ip_q * R_L;

    double gm = (Ip_q > 0) ? (Ex * Ip_q / std::max(0.001, (V_supply / Mu) + V_gk)) : 0.0;
    r_p = Mu / (gm + 1e-12);
    R_out = (r_p * R_L) / (r_p + R_L);

    double Zk_ac = isBypassed ? 0.0 : R_k;
    gain = (Mu * R_L) / (R_L + r_p + Zk_ac * (Mu + 1.0));
}

// Processes a sample using the calculated Thevenin model and soft saturation.
double TubeStage::process(double V_in_ac, double R_load) {
    double loaded_gain = gain * (R_load / (R_out + R_load));
    double V_out_linear = -V_in_ac * loaded_gain;

    double V_headroom = V_supply - Vp_dc;
    double V_footroom = Vp_dc;
    double k = 20.0; // Softness factor: higher is softer

    double V_out_ac = 0.0;
    if (V_out_linear > 0) {
        V_out_ac = V_headroom * (V_out_linear / sqrt(V_out_linear * V_out_linear + V_headroom * V_headroom / (k * k)));
    } else {
        V_out_ac = V_footroom * (V_out_linear / sqrt(V_out_linear * V_out_linear + V_footroom * V_footroom / (k * k)));
    }

    return Vp_dc + V_out_ac;
}
