#pragma once
#include <algorithm>
#include <cmath>

// Model parameters for a 12AX7 triode based on the SPICE model
// .SUBCKT 12AX7 1 2 3; A G C;
// X1 1 2 3 TriodeK MU=96.20 EX=1.437 KG1=613.4 KP=740.3 KVB=1672. RGI=2000
// Note: This implementation simplifies grid current for stability, using a simple diode model.

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
        State s = {0.0, 0.0, 1e-9, 1e-9, 1e-9}; // Initialize with tiny conductances for stability

        // Plate Current Calculation (Koren Model)
        double e1_num = v_p / KP;
        double e1_den = 1.0 / MU + v_g / sqrt(KVB + v_p * v_p);

        double e1 = 0.0;
        if (e1_den > -0.999) { // Avoid log of negative/zero
             e1 = e1_num * log(1.0 + exp(KP * e1_den));
        } else {
             e1 = e1_num * KP * e1_den; // Linear approximation for large negative args
        }

        if (e1 < 0) e1 = 0; // Current doesn't flow backward

        s.ip = (2.0 * pow(e1, EX)) / (2.0 * KG1); // Simplified from 0.5*(pwr+pwrs) since e1>=0

        // Derivatives for Newton-Raphson
        if (s.ip > 1e-12) { // Only calculate if there is current
            double exp_val = exp(KP * e1_den);
            double common_term = EX * pow(e1, EX - 1.0) / KG1;
            double d_e1_dvp = (e1 / v_p) - (v_p * v_g * KP * exp_val) / ( (1.0+exp_val) * sqrt(pow(KVB + v_p*v_p, 3)) );
            double d_e1_dvg = (v_p * KP * exp_val) / (KP * (1.0 + exp_val) * sqrt(KVB + v_p*v_p));

            s.g_p = common_term * d_e1_dvp;
            s.g_g = common_term * d_e1_dvg;
        }

        // Grid Current Calculation (simplified diode model + RGI)
        if (v_g > 0) {
            double v_diode = v_g; // Simplified: assume voltage is directly on diode
            s.ig = IS * (exp(v_diode / VT) - 1.0);
            s.g_ig = IS / VT * exp(v_diode / VT);
        }

        // Add RGI contribution
        s.ig += v_g / RGI;
        s.g_ig += 1.0 / RGI;

        // Clamp conductances to prevent instability
        s.g_p = std::max(1e-9, s.g_p);
        s.g_g = std::max(1e-9, s.g_g);
        s.g_ig = std::max(1e-9, s.g_ig);

        return s;
    }
};
