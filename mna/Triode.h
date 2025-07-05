#pragma once
#include <cmath>
#include <algorithm> // For std::max
#include <unordered_set>

// Model parameters for a 12AX7 triode based on the SPICE model
// .SUBCKT 12AX7 1 2 3; A G C;
// X1 1 2 3 TriodeK MU=96.20 EX=1.437 KG1=613.4 KP=740.3 KVB=1672. RGI=2000
class Triode {
    // Koren Tube Model Parameters
    static constexpr double MU = 96.20;  // Amplification Factor
    static constexpr double EX = 1.437;  // Exponent
    static constexpr double KG1 = 613.4; // Grid 1 Factor
    static constexpr double KP = 740.3;  // Plate Factor
    static constexpr double KVB = 1672.0;// Cathode-Plate Factor
    static constexpr double RGI = 2000.0;// Grid Input Resistance

    // Grid Diode Parameters (simplified)
    static constexpr double IS = 1e-9;   // Saturation Current
    static constexpr double VT = 25.852e-3; // 0.026;  // Thermal Voltage

public:
    struct State {
        double ip; // Plate Current
        double ig; // Grid Current
        double g_p; // Plate Conductance (dIp/dVp)
        double g_g; // Transconductance (dIp/dVg)
        double g_ig; // Grid Conductance (dIg/dVg)
    };

private:
    // --- Linear Interpolation Helpers ---
    static inline double lerp(double a, double b, double t) {
        return a + t * (b - a);
    }

    static State interpolateState(const State& s1, const State& s2, double t) {
        State result;
        result.ip = lerp(s1.ip, s2.ip, t);
        result.ig = lerp(s1.ig, s2.ig, t);
        result.g_p = lerp(s1.g_p, s2.g_p, t);
        result.g_g = lerp(s1.g_g, s2.g_g, t);
        result.g_ig = lerp(s1.g_ig, s2.g_ig, t);
        return result;
    }

    // --- LUT Configuration and State ---
    static constexpr int VG_STEPS = 1000;
    static constexpr int VP_STEPS = 1000;
    static constexpr double MIN_VG = -50.0, MAX_VG = 2.0;
    static constexpr double MIN_VP = 0.0, MAX_VP = 1000.0;

    // Using 'inline' for C++17+ allows definition in the header
    inline static std::vector<State> lut;
    inline static double vg_step, vp_step, inv_vg_step, inv_vp_step;
    inline static bool lut_initialized = false;

    // Calculates the state of the triode given plate-cathode and grid-cathode voltages
    static State calculateDirectly(double v_p, double v_g) {
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

        // --- DEFINITIVE, STABLE GRID CURRENT MODEL ---
        // This version prevents numerical overflow by linearizing the diode at high voltages.
        if (v_g > 0) {
            // Smooth transition point
            const double V_crit = 0.5; // A good point for the transition

            // Calculate the linear part's "y-intercept" based on the exponential curve at V_crit
            // This ensures a smooth connection.
            const double i_crit = IS * (exp(V_crit / VT) - 1.0);
            const double g_crit = (IS / VT) * exp(V_crit / VT);
            const double y_intercept = i_crit - g_crit * V_crit;

            // This is a common "softplus"-like log-sum-exp trick for smooth blending.
            // It smoothly blends between exp(v_g/VT) and g_crit.
            double soft_exp_arg = (v_g - V_crit) / (2.0 * VT);
            if (soft_exp_arg > 100) soft_exp_arg = 100; // Prevent overflow

            double blend_factor = 1.0 + exp(soft_exp_arg);
            double log_blend = VT * log(blend_factor);

            s.ig = y_intercept + g_crit * V_crit + g_crit * 2.0 * log_blend;
            s.g_ig = g_crit * (exp(soft_exp_arg) / blend_factor);

            // RGI is always in parallel with this diode behavior
            s.ig += v_g / RGI;
            s.g_ig += 1.0 / RGI;
        }

        // Clamp conductances to prevent instability
        s.g_p = std::max(1e-12, s.g_p);
        s.g_g = std::max(1e-12, s.g_g);
        s.g_ig = std::max(1e-12, s.g_ig);

        return s;
    }

public:
#ifdef USE_TRIODE_CACHE
    /**
     * @brief Initializes the lookup table. Must be called once before any audio processing.
     */
    static void initializeLUT() {
        if (lut_initialized) return;

        std::cout << "Initializing 12AX7 lookup table..." << std::endl;

        lut.resize(VG_STEPS * VP_STEPS);
        vg_step = (MAX_VG - MIN_VG) / (VG_STEPS - 1);
        vp_step = (MAX_VP - MIN_VP) / (VP_STEPS - 1);
        inv_vg_step = 1.0 / vg_step;
        inv_vp_step = 1.0 / vp_step;

        for (int i = 0; i < VG_STEPS; ++i) {
            double v_g = MIN_VG + i * vg_step;
            for (int j = 0; j < VP_STEPS; ++j) {
                double v_p = MIN_VP + j * vp_step;
                lut[i * VP_STEPS + j] = calculateDirectly(v_p, v_g);
            }
        }
        lut_initialized = true;
        std::cout << "12AX7 lookup table initialized." << std::endl;
    }

    /**
     * @brief Calculates the triode state using fast bilinear interpolation from the LUT.
     */
    static State calculate(double v_p, double v_g) {
        if (v_p < min_v_p) min_v_p = v_p;
        if (v_p > max_v_p) max_v_p = v_p;
        if (v_g < min_v_g) min_v_g = v_g;
        if (v_g > max_v_g) max_v_g = v_g;

        // Clamp inputs to the defined range of the LUT
        v_g = std::clamp(v_g, MIN_VG, MAX_VG);
        v_p = std::clamp(v_p, MIN_VP, MAX_VP);

        // Find grid position and interpolation weight
        double vg_pos = (v_g - MIN_VG) * inv_vg_step;
        int vg_idx = std::clamp(static_cast<int>(vg_pos), 0, VG_STEPS - 2);
        double vg_t = vg_pos - vg_idx;

        // Find plate position and interpolation weight
        double vp_pos = (v_p - MIN_VP) * inv_vp_step;
        int vp_idx = std::clamp(static_cast<int>(vp_pos), 0, VP_STEPS - 2);
        double vp_t = vp_pos - vp_idx;

        // Get the four corner points from the LUT for interpolation
        const State& s00 = lut[vg_idx * VP_STEPS + vp_idx];
        const State& s10 = lut[(vg_idx + 1) * VP_STEPS + vp_idx];
        const State& s01 = lut[vg_idx * VP_STEPS + (vp_idx + 1)];
        const State& s11 = lut[(vg_idx + 1) * VP_STEPS + (vp_idx + 1)];

        // Bilinear interpolation
        State s0 = interpolateState(s00, s01, vp_t); // Lerp along plate voltage axis for lower grid voltage
        State s1 = interpolateState(s10, s11, vp_t); // Lerp along plate voltage axis for upper grid voltage
        return interpolateState(s0, s1, vg_t);       // Lerp along grid voltage axis
    }
#else
    static State calculate(double v_p, double v_g) {
        if (v_p < min_v_p) min_v_p = v_p;
        if (v_p > max_v_p) max_v_p = v_p;
        if (v_g < min_v_g) min_v_g = v_g;
        if (v_g > max_v_g) max_v_g = v_g;
        return calculateDirectly(v_p, v_g);
    }

    static void initializeLUT() {}
#endif
    inline static double min_v_p = 1e6, max_v_p = -1e6;
    inline static double min_v_g = 1e6, max_v_g = -1e6;

    static void printStats() {
        std::cout << "v_p min " << min_v_p << " max " << max_v_p << std::endl;
        std::cout << "v_g min " << min_v_g << " max " << max_v_g << std::endl;
    }
};
