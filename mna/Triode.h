#pragma once
#include <cmath>
#include <algorithm> // For std::max

// Model parameters for a 12AX7 triode based on the SPICE model
// .SUBCKT 12AX7 1 2 3; A G C;
// X1 1 2 3 TriodeK MU=96.20 EX=1.437 KG1=613.4 KP=740.3 KVB=1672. RGI=2000
class Triode {
private:
    // Koren Tube Model Parameters
    static constexpr double MU = 96.20;  // Amplification Factor
    static constexpr double EX = 1.437;  // Exponent
    static constexpr double KG1 = 613.4; // Grid 1 Factor
    static constexpr double KP = 740.3;  // Plate Factor
    static constexpr double KVB = 1672.0;// Cathode-Plate Factor
    static constexpr double RGI = 2000.0;// Grid Input Resistance

    // Grid Diode Parameters (simplified)
    static constexpr double IS = 1e-9;   // Saturation Current
    static constexpr double VT = 0.026;  // Thermal Voltage

public:
    struct State {
        double ip; // Plate Current
        double ig; // Grid Current
        double g_p; // Plate Conductance (dIp/dVp)
        double g_g; // Transconductance (dIp/dVg)
        double g_ig; // Grid Conductance (dIg/dVg)
    };

    // Calculates the state of the triode given plate-cathode and grid-cathode voltages
    static State calculate(double v_p, double v_g) {
        State s = {0.0, 0.0, 1e-12, 1e-12, 1e-12}; // Initialize with tiny conductances for stability

        // --- Plate Current Calculation (Koren Model) ---
        // To avoid floating point issues, ensure v_p is slightly positive.
        // The sqrt(KVB + v_p*v_p) term handles negative v_p gracefully, but v_p=0 can be problematic.
        double v_p_safe = std::max(v_p, 1e-6);

        double e1_den_sqrt = sqrt(KVB + v_p_safe * v_p_safe);
        double e1_den = v_g / e1_den_sqrt;

        // The 'soft-rectifier' log(1 + exp(x)) function
        double e1_exp_arg = KP * (1.0 / MU + e1_den);
        double e1 = 0;
        if (e1_exp_arg < 100) { // Prevent overflow in exp()
            e1 = (v_p_safe / KP) * log(1.0 + exp(e1_exp_arg));
        } else {
            e1 = (v_p_safe / KP) * e1_exp_arg; // Linear approximation for large x
        }

        // The PWRS function from SPICE is 0 for negative base.
        if (e1 < 0) e1 = 0;
        s.ip = pow(e1, EX) / KG1;

        // --- Derivatives for Newton-Raphson (Corrected) ---
        if (s.ip > 1e-12) {
            double exp_val = exp(e1_exp_arg);
            double log_val = log(1.0 + exp_val);
            double common_term1 = exp_val / (1.0 + exp_val);

            // Derivative of e1 with respect to v_p
            double de1_dvp = (1.0 / KP) * log_val - (v_p_safe * v_g * common_term1) / pow(e1_den_sqrt, 3);

            // Derivative of e1 with respect to v_g
            double de1_dvg = (v_p_safe / e1_den_sqrt) * common_term1;

            // Final derivatives of plate current Ip
            double common_term2 = EX * pow(e1, EX - 1.0) / KG1;
            s.g_p = common_term2 * de1_dvp;
            s.g_g = common_term2 * de1_dvg;
        }

        // --- Grid Current Calculation ---
        if (v_g > 0) {
            s.ig = IS * (exp(v_g / VT) - 1.0);
            s.g_ig = (IS / VT) * exp(v_g / VT);
        }
        s.ig += v_g / RGI;
        s.g_ig += 1.0 / RGI;

        // Clamp conductances to prevent instability
        s.g_p = std::max(1e-12, s.g_p);
        s.g_g = std::max(1e-12, s.g_g);
        s.g_ig = std::max(1e-12, s.g_ig);

        return s;
    }
};
