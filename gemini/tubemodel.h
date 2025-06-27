#pragma once
#include <cmath> // For std::fma if you want to use it, though not required here.

/**
 * @class TubeModel
 * @brief Implements a 5-segment piecewise polynomial model of a 12AX7 triode tube.
 *
 * This class is designed for high-performance, real-time audio processing.
 * All model parameters are compile-time constants, and the processSample
 * method is allocation-free and suitable for use in an audio thread.
 *
 * The model maps an input grid voltage (Vg) to an output plate voltage (Vp).
 */
class TubeModel
{
public:
    /**
     * @brief Processes a single input sample (Vg) to produce an output sample (Vp).
     * @param vg The input grid voltage.
     * @return The corresponding plate voltage based on the model.
     */
    static double processSample(double vg)
    {
        // 1. Handle the clamping regions outside the main model range
        if (vg < BOUNDARY_0) {
            return CLAMP_HIGH_V;
        }
        else if (vg >= BOUNDARY_5) {
            return CLAMP_LOW_V;
        }

        // 2. Select the correct polynomial based on the segment
        else if (vg < BOUNDARY_1) { // Segment 1
            return evalQuadratic(vg, SEG1_COEFFS);
        }
        else if (vg < BOUNDARY_2) { // Segment 2
            return evalCubic(vg, SEG2_COEFFS);
        }
        else if (vg < BOUNDARY_3) { // Segment 3
            return evalQuadratic(vg, SEG3_COEFFS);
        }
        else if (vg < BOUNDARY_4) { // Segment 4
            return evalQuadratic(vg, SEG4_COEFFS);
        }
        else { // Segment 5 is the final remaining case
            return evalQuadratic(vg, SEG5_COEFFS);
        }
    }

private:
    // --- Model Parameters (Compile-time constants) ---

    // Clamping Voltages
    static constexpr double CLAMP_HIGH_V = 405.0;
    static constexpr double CLAMP_LOW_V  = 15.682660;

    // Optimized Segment Boundaries
    static constexpr double BOUNDARY_0 = -8.0;       // Start of Segment 1
    static constexpr double BOUNDARY_1 = -5.677862;  // Start of Segment 2
    static constexpr double BOUNDARY_2 = -2.481141;  // Start of Segment 3
    static constexpr double BOUNDARY_3 = 1.722849;   // Start of Segment 4
    static constexpr double BOUNDARY_4 = 6.116238;   // Start of Segment 5
    static constexpr double BOUNDARY_5 = 17.628460;  // End of Segment 5

    // Polynomial Coefficients (in descending power order: a, b, c, d)

    // Segment 1: Vg in [-8.000000, -5.677862) -> Quadratic
    static constexpr double SEG1_COEFFS[] = { -4.95201870e-01, -7.56116037e+00, 3.76137455e+02 };

    // Segment 2: Vg in [-5.677862, -2.481141) -> Cubic
    static constexpr double SEG2_COEFFS[] = { 1.23478003e+00, 8.90103893e+00, -2.02810020e+01, 2.27017396e+02 };

    // Segment 3: Vg in [-2.481141, 1.722849) -> Quadratic
    static constexpr double SEG3_COEFFS[] = { -2.89953394e-01, -4.30851555e+01, 2.08157285e+02 };

    // Segment 4: Vg in [1.722849, 6.116238) -> Quadratic
    static constexpr double SEG4_COEFFS[] = { 4.66301432e+00, -6.01515846e+01, 2.22858724e+02 };

    // Segment 5: Vg in [6.116238, 17.628460] -> Quadratic
    static constexpr double SEG5_COEFFS[] = { 1.66816510e-01, -5.15195138e+00, 5.46632961e+01 };

    // --- Efficient Polynomial Evaluation using Horner's Method ---

    /**
     * @brief Evaluates ax^2 + bx + c
     */
    static inline double evalQuadratic(double x, const double* c) {
        // (a*x + b)*x + c
        return (c[0] * x + c[1]) * x + c[2];
    }

    /**
     * @brief Evaluates ax^3 + bx^2 + cx + d
     */
    static inline double evalCubic(double x, const double* c) {
        // ((a*x + b)*x + c)*x + d
        return ((c[0] * x + c[1]) * x + c[2]) * x + c[3];
    }
};
