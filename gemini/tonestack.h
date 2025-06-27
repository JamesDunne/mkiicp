#ifndef TONESTACK_FILTER_NODAL_H
#define TONESTACK_FILTER_NODAL_H

#include <cmath>

/**
 * @class PassiveToneStack
 * @brief A realtime digital filter simulating a passive tone stack via Nodal Analysis.
 *
 * This version uses the Modified Nodal Analysis (MNA) approach, solving the circuit
 * equations at runtime for maximum accuracy and robustness. The circuit is solved
 * as a system of 7 linear equations representing the 7 key unknown node voltages.
 *
 * The input samples passed to processSample() are expected to be in the
 * standard floating-point audio range of [-1.0, 1.0]. The output will be in the
 * same range and is guaranteed to be passive (no amplification).
 */
class PassiveToneStack {
public:
    PassiveToneStack() {
        reset();
    }

    /**
     * @brief Resets the filter's internal state variables to zero.
     */
    void reset() {
        m_s_c56 = m_s_c4 = m_s_c3 = m_s_c13b = 0.0;
    }

    /**
     * @brief Sets all filter parameters and recalculates the system matrix.
     *
     * This method is called at setup and whenever a control parameter changes.
     * It prepares the solver for the new component values.
     *
     * @param sampleRate The audio processing sample rate in Hz.
     * @param treble Control for the treble potentiometer (0.0 to 1.0).
     * @param bass Control for the bass variable resistor (0.0 to 1.0).
     * @param mid Control for the mid variable resistor (0.0 to 1.0).
     * @param volume Control for the volume potentiometer (0.0 to 1.0).
     */
    void setParameters(double sampleRate, double treble, double bass, double mid, double volume) {
        // Clamp input parameters to the [0, 1] range
        treble = std::clamp(treble, 0.0, 1.0);
        bass = std::clamp(bass, 0.0, 1.0);
        mid = std::clamp(mid, 0.0, 1.0);
        volume = std::clamp(volume, 0.0, 1.0);

        // --- Stabilize Potentiometer Travel ---
        // To prevent the MNA matrix from becoming singular, we must prevent any
        // resistance from becoming a perfect zero. We map the user's [0,1] range
        // to a slightly smaller internal range, e.g., [0.001, 0.999].
        const double travel_min = 0.001;
        const double travel_max = 1.0 - travel_min;
        double treble_t = treble * travel_max + travel_min;
        double volume_t = volume * travel_max + travel_min;

        // --- Define Component Values and Epsilon for Stability ---
        const double epsilon = 1e-12;

        const double C56 = 1.0e-9;
        const double C4 = 1.0e-7;
        const double C3 = 4.7e-8;
        const double C13B = 1.8e-10;

        // --- Calculate Component Conductances (G = 1/R) ---
        m_g_r5 = 1.0 / 100000.0;
        const double g_r5a = 1.0 / 100000.0;

        double r_treble_a = 250000.0 * (1.0 - treble_t);
        double r_treble_c = 250000.0 * treble_t;
        double g_treble_a = 1.0 / (r_treble_a); // No epsilon needed due to travel clamp
        double g_treble_c = 1.0 / (r_treble_c);

        double r_bass_a = 250000.0 * bass;
        double g_bass_a = 1.0 / (r_bass_a + epsilon); // Epsilon needed as this can be 0

        double r_mid_a = 10000.0 * mid;
        double g_mid_a = 1.0 / (r_mid_a + epsilon); // Epsilon needed

        double r_vol_a = 1000000.0 * (1.0 - volume_t);
        double r_vol_c = 1000000.0 * volume_t;
        double g_vol_a = 1.0 / (r_vol_a);
        double g_vol_c = 1.0 / (r_vol_c);

        // --- Update Capacitor Conductances for Trapezoidal Integration ---
        m_g_c56 = 2.0 * C56 * sampleRate;
        m_g_c4 = 2.0 * C4 * sampleRate;
        m_g_c3 = 2.0 * C3 * sampleRate;
        m_g_c13b = 2.0 * C13B * sampleRate;

        // --- Build the 7x7 MNA System Matrix 'A' ---
        double A[7][7] = {0};

        // KCL at N005
        A[0][0] = g_treble_a + m_g_c56; A[0][1] = -g_treble_a;
        // KCL at N007
        A[1][0] = -g_treble_a; A[1][1] = g_treble_a + g_treble_c + g_r5a; A[1][3] = -g_treble_c; A[1][5] = -g_r5a;
        // KCL at N015
        A[2][2] = m_g_r5 + m_g_c4 + m_g_c3; A[2][3] = -m_g_c4; A[2][4] = -m_g_c3;
        // KCL at N016
        A[3][1] = -g_treble_c; A[3][2] = -m_g_c4; A[3][3] = g_treble_c + g_bass_a + m_g_c4; A[3][4] = -g_bass_a;
        // KCL at N024
        A[4][2] = -m_g_c3; A[4][3] = -g_bass_a; A[4][4] = g_bass_a + g_mid_a + m_g_c3;
        // KCL at N008
        A[5][1] = -g_r5a; A[5][5] = g_r5a + g_vol_a + m_g_c13b; A[5][6] = -g_vol_a - m_g_c13b;
        // KCL at N020
        A[6][5] = -g_vol_a - m_g_c13b; A[6][6] = g_vol_a + g_vol_c + m_g_c13b;

        // Pre-calculate the inverse of the system matrix
        invertMatrix(A, m_A_inv);
    }

    /**
     * @brief Processes one audio sample. The input is assumed to be Vin.
     * @param vin The input sample, typically in range [-1.0, 1.0].
     * @return The filtered output sample from node N020.
     */
    double processSample(double vin) {
        // Build the 'b' vector for Ax = b
        double b[7];
        b[0] = m_s_c56 + vin * m_g_c56;
        b[1] = 0.0;
        b[2] = vin * m_g_r5 + m_s_c4 + m_s_c3;
        b[3] = -m_s_c4;
        b[4] = -m_s_c3;
        b[5] = m_s_c13b;
        b[6] = -m_s_c13b;

        // Solve for the node voltages: v = A_inv * b
        double v[7];
        for (int i = 0; i < 7; ++i) {
            v[i] = 0.0;
            for (int j = 0; j < 7; ++j) {
                v[i] += m_A_inv[i][j] * b[j];
            }
        }

        // The output is the voltage at node N020 (index 6)
        double output = v[6];

        // Update capacitor state variables for the next time step
        m_s_c56  = 2.0 * m_g_c56  * (v[0] - vin) - m_s_c56;
        m_s_c4   = 2.0 * m_g_c4   * (v[2] - v[3]) - m_s_c4;
        m_s_c3   = 2.0 * m_g_c3   * (v[2] - v[4]) - m_s_c3;
        m_s_c13b = 2.0 * m_g_c13b * (v[5] - v[6]) - m_s_c13b;

        return static_cast<float>(output);
    }

private:
    // System matrix inverse
    double m_A_inv[7][7];

    // Component conductances needed in the process loop
    double m_g_r5;
    double m_g_c56, m_g_c4, m_g_c3, m_g_c13b;

    // State variables for capacitors
    double m_s_c56, m_s_c4, m_s_c3, m_s_c13b;

    // A simple 7x7 matrix inverter using Gaussian elimination.
    // For production code, a more robust numerical library (e.g., Eigen) is recommended.
    static void invertMatrix(double a[7][7], double inv[7][7]) {
        double temp[7][14];
        for (int i = 0; i < 7; i++) {
            for (int j = 0; j < 7; j++) {
                temp[i][j] = a[i][j];
                temp[i][j + 7] = (i == j) ? 1.0 : 0.0;
            }
        }

        for (int i = 0; i < 7; i++) {
            double pivot = temp[i][i];
            for (int j = i; j < 14; j++) {
                temp[i][j] /= pivot;
            }
            for (int k = 0; k < 7; k++) {
                if (i != k) {
                    double factor = temp[k][i];
                    for (int j = i; j < 14; j++) {
                        temp[k][j] -= factor * temp[i][j];
                    }
                }
            }
        }

        for (int i = 0; i < 7; i++) {
            for (int j = 0; j < 7; j++) {
                inv[i][j] = temp[i][j + 7];
            }
        }
    }
};

#endif // TONESTACK_FILTER_NODAL_H
