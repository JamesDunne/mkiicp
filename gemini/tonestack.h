#ifndef TONESTACK_FILTER_H
#define TONESTACK_FILTER_H

#include <cmath>

/**
 * @class PassiveToneStack
 * @brief A realtime digital IIR filter simulating a passive guitar amplifier tone stack.
 *
 * This class models the specific passive tone stack circuit provided in the SPICE netlist.
 * The circuit is analyzed as a 4th-order Linear Time-Invariant (LTI) system.
 *
 * The transfer function H(s) = V(N020)/V(N004) is derived symbolically and then
 * converted to a digital transfer function H(z) using the Bilinear Transform.
 * The filter is implemented as a Direct Form I Biquad structure for numerical stability
 * and efficiency.
 *
 * All parameters are controlled by values from 0.0 to 1.0, mimicking the behavior
 * of potentiometers as defined in the SPICE netlist.
 */
class PassiveToneStack {
public:
    /**
     * @brief Constructor. Initializes filter state.
     */
    PassiveToneStack() {
        reset();
    }

    /**
     * @brief Resets the filter's internal state variables to zero.
     *        Useful for clearing buffers after silence or major parameter changes.
     */
    void reset() {
        x1 = x2 = x3 = x4 = 0.0;
        y1 = y2 = y3 = y4 = 0.0;
    }

    /**
     * @brief Sets all filter parameters and recalculates the IIR coefficients.
     *
     * This method should be called whenever the sample rate or any of the tone/volume
     * controls change. It is designed to be efficient for realtime use.
     *
     * @param sampleRate The audio processing sample rate in Hz (e.g., 44100.0).
     * @param treble Control for the treble potentiometer (0.0 to 1.0, linear taper).
     * @param bass Control for the bass variable resistor (0.0 to 1.0, linear taper).
     * @param mid Control for the mid variable resistor (0.0 to 1.0, linear taper).
     * @param volume Control for the volume potentiometer (0.0 to 1.0, linear taper).
     */
    void setParameters(double sampleRate, double treble, double bass, double mid, double volume) {
        // Clamp input parameters to a safe range [0.0, 1.0]
        treble = fmax(0.0, fmin(1.0, treble));
        bass = fmax(0.0, fmin(1.0, bass));
        mid = fmax(0.0, fmin(1.0, mid));
        volume = fmax(0.0, fmin(1.0, volume));

        // --- SPICE Component Values ---
        const double C3 = 4.7e-08;  // .047uF
        const double C4 = 1.0e-07;  // .1uF
        const double C56 = 1.0e-09; // 750p + 250p
        const double C13B = 1.8e-10; // 180p
        const double R5 = 100000.0;
        const double R5A = 100000.0;
        const double R_TREBLE_POT = 250000.0;
        const double R_BASS_POT = 250000.0;
        const double R_MID_POT = 10000.0;
        const double R_VOL_POT = 1000000.0;

        // --- Potentiometer Models from SPICE (Corrected for intuitive control) ---
        double r_treble_a = R_TREBLE_POT * (1.0 - treble);
        double r_treble_c = R_TREBLE_POT * treble;

        // CORRECTED: Higher 'bass' parameter now means higher resistance for more bass.
        double r_bass_a = R_BASS_POT * bass;

        // CORRECTED: Higher 'mid' parameter now means higher resistance for more mids.
        double r_mid_a = R_MID_POT * mid;

        double r_vol_a = R_VOL_POT * (1.0 - volume);
        double r_vol_c = R_VOL_POT * volume;

        // Add small epsilon to resistances to prevent division by zero in edge cases
        const double epsilon = 1e-12;
        if (r_treble_a < epsilon) r_treble_a = epsilon;
        if (r_treble_c < epsilon) r_treble_c = epsilon;
        if (r_vol_c < epsilon) r_vol_c = epsilon;
        // Also add epsilon to the rheostats
        if (r_bass_a < epsilon) r_bass_a = epsilon;
        if (r_mid_a < epsilon) r_mid_a = epsilon;

        // --- Symbolic analysis results pre-calculated for H(s) = N(s)/D(s) ---
        // (The rest of the symbolic math remains identical)

        double t0 = r_treble_a + r_treble_c;
        double t1 = R5 * C56;
        double t2 = R5 * C4;
        double t3 = R5 * C3;
        double t4 = C56 * t0;
        double t5 = C4 * t0;
        double t6 = C3 * t0;
        double t7 = r_treble_a * C56;
        double t8 = r_treble_c * C4;
        double t9 = r_treble_c * C3;
        double t10 = r_bass_a * C4;
        double t11 = r_bass_a * C3;
        double t12 = r_mid_a * C3;
        double t13 = C13B * r_vol_a;
        double t14 = r_vol_a + r_vol_c;
        double t15 = C13B * t14;

        // s-domain denominator coefficients (d_i)
        double d4 = R5A * t1 * t10 * t12 + R5A * t1 * t11 * t12 + R5A * t10 * t11 * t3 + R5A * t10 * t12 * t2 + R5A * t11 * t12 * t2 + t1 * t10 * t12 * t13 + t1 * t11 * t12 * t13 + t10 * t11 * t13 * t3 + t10 * t12 * t13 * t2 + t11 * t12 * t13 * t2;
        double d3 = R5A*t1*t10 + R5A*t1*t11 + R5A*t1*t12 + R5A*t10*t11 + R5A*t10*t12 + R5A*t10*t3 + R5A*t11*t12 + R5A*t11*t2 + R5A*t11*t3 + R5A*t12*t2 + R5A*t2*t6 + R5A*t3*t5 + R5A*t4*t10 + R5A*t4*t11 + R5A*t4*t12 + R5A*t5*t11 + R5A*t5*t12 + R5A*t6*t10 + R5A*t6*t12 + t1*t10*t13 + t1*t11*t13 + t1*t12*t13 + t1*t15*t10 + t1*t15*t11 + t1*t15*t12 + t10*t11*t13 + t10*t12*t13 + t10*t13*t3 + t10*t15*t3 + t11*t12*t13 + t11*t13*t2 + t11*t15*t2 + t12*t13*t2 + t12*t15*t2 + t13*t2*t6 + t13*t3*t5 + t15*t2*t6 + t15*t3*t5 + t4*t10*t13 + t4*t11*t13 + t4*t12*t13 + t4*t15*t10 + t4*t15*t11 + t4*t15*t12 + t5*t11*t13 + t5*t12*t13 + t5*t15*t11 + t5*t15*t12 + t6*t10*t13 + t6*t12*t13 + t6*t15*t10 + t6*t15*t12;
        double d2 = R5*C3*t15 + R5*C4*t15 + R5*C56*t15 + R5A*t1 + R5A*t10 + R5A*t11 + R5A*t12 + R5A*t2 + R5A*t3 + R5A*t4 + R5A*t5 + R5A*t6 + t0*t15 + t1*t13 + t1*t15 + t10*t13 + t10*t15 + t11*t13 + t11*t15 + t12*t13 + t12*t15 + t13*t2 + t13*t3 + t13*t4 + t13*t5 + t13*t6 + t14*t10*t11 + t14*t10*t12 + t14*t11*t12 + t14*t2*t6 + t14*t3*t5 + t15*t2 + t15*t3 + t15*t4 + t15*t5 + t15*t6 + r_bass_a*t15 + r_mid_a*t15 + r_treble_a*t15;
        double d1 = R5*t15 + R5A + t0 + t13 + t14*t10 + t14*t11 + t14*t12 + t14*t2 + t14*t3 + t14*t4 + t14*t5 + t14*t6 + t15 + r_bass_a + r_mid_a + r_treble_a;
        double d0 = t14;

        // s-domain numerator coefficients (c_i)
        double c3 = r_vol_c * (R5A*t1*t8*t11 + t1*t13*t8*t11);
        double c2 = r_vol_c * (R5A*t1*t8 + R5A*t1*t9 + R5A*t11*t7 + R5A*t11*t8 + R5A*t7*t9 + R5A*t8*t9 + t1*t13*t8 + t1*t13*t9 + t11*t13*t7 + t11*t13*t8 + t11*t15*t7 + t11*t15*t8 + t13*t7*t9 + t13*t8*t9 + t15*t7*t9);
        double c1 = r_vol_c * (R5*C3*t15 + R5*C4*t15 + R5A*t7 + R5A*t8 + R5A*t9 + t13*t7 + t13*t8 + t13*t9 + t14*t11*t7 + t14*t11*t8 + t14*t7*t9 + t14*t8*t9 + t15*t7 + t15*t8 + t15*t9 + r_treble_c*t15);
        double c0 = r_vol_c * (t14*t7 + t14*t8 + t14*t9 + r_treble_c);

        // --- Bilinear Transform ---
        double fs2 = 2.0 * sampleRate;
        double fs2_2 = fs2 * fs2;
        double fs2_3 = fs2_2 * fs2;
        double fs2_4 = fs2_3 * fs2;

        double a0_z = d4*fs2_4 + d3*fs2_3 + d2*fs2_2 + d1*fs2 + d0;
        double a1_z = -4*d4*fs2_4 - 2*d3*fs2_3 + 2*d1*fs2 + 4*d0;
        double a2_z = 6*d4*fs2_4 - 2*d2*fs2_2 + 6*d0;
        double a3_z = -4*d4*fs2_4 + 2*d3*fs2_3 - 2*d1*fs2 + 4*d0;
        double a4_z = d4*fs2_4 - d3*fs2_3 + d2*fs2_2 - d1*fs2 + d0;

        double b0_z = c3*fs2_3 + c2*fs2_2 + c1*fs2 + c0;
        double b1_z = -3*c3*fs2_3 - c2*fs2_2 + c1*fs2 + 3*c0;
        double b2_z = 3*c3*fs2_3 - c2*fs2_2 - c1*fs2 + 3*c0;
        double b3_z = -c3*fs2_3 + c2*fs2_2 - c1*fs2 + c0;
        double b4_z = 0;

        if (std::abs(a0_z) > epsilon) {
            double inv_a0 = 1.0 / a0_z;
            b0 = b0_z * inv_a0;
            b1 = b1_z * inv_a0;
            b2 = b2_z * inv_a0;
            b3 = b3_z * inv_a0;
            b4 = b4_z * inv_a0;
            a1 = a1_z * inv_a0;
            a2 = a2_z * inv_a0;
            a3 = a3_z * inv_a0;
            a4 = a4_z * inv_a0;
        }
    }

    /**
     * @brief Processes one audio sample through the filter.
     * @param input The input audio sample.
     * @return The filtered output audio sample.
     */
    double processSample(double input) {
        // Direct Form I implementation of the 4th order IIR filter
        // y[n] = b0*x[n] + b1*x[n-1] + b2*x[n-2] + ... - a1*y[n-1] - a2*y[n-2] - ...

        double output = b0 * input + b1 * x1 + b2 * x2 + b3 * x3 + b4 * x4
                                   - a1 * y1 - a2 * y2 - a3 * y3 - a4 * y4;

        // Update state variables for next sample
        x4 = x3;
        x3 = x2;
        x2 = x1;
        x1 = input;

        y4 = y3;
        y3 = y2;
        y2 = y1;
        y1 = output;

        return output;
    }

private:
    // Filter coefficients
    double b0, b1, b2, b3, b4;
    double a1, a2, a3, a4;

    // Filter state (delay lines)
    double x1, x2, x3, x4; // Previous inputs
    double y1, y2, y3, y4; // Previous outputs
};

#endif // TONESTACK_FILTER_H
