#pragma once

#include <array>    // For std::array

/**
 * @struct TubeStageV1A
 * @brief A data-only struct containing all the compile-time parameters for our 12AX7 model.
 *
 * This struct will be passed as a template argument to the TubeModel class.
 * To create a new tube model, simply copy this struct, rename it, and replace the values.
 */
struct TubeStageV1A
{
    // --- Clamping Voltages ---
    static constexpr double CLAMP_HIGH_V = 405.0;
    static constexpr double CLAMP_LOW_V  = 15.682660;

    // --- Optimized Segment Boundaries ---
    static constexpr double BOUNDARY_0 = -8.0;
    static constexpr double BOUNDARY_1 = -5.677862;
    static constexpr double BOUNDARY_2 = -2.481141;
    static constexpr double BOUNDARY_3 = 1.722849;
    static constexpr double BOUNDARY_4 = 6.116238;
    static constexpr double BOUNDARY_5 = 17.628460;

    // --- Polynomial Coefficients ---
    // Using std::array for type safety and to bundle size information.

    // Segment 1: Vg in [-8.0, -5.677862) -> Quadratic
    static constexpr std::array<double, 3> SEG1_COEFFS = { -4.95201870e-01, -7.56116037e+00, 3.76137455e+02 };

    // Segment 2: Vg in [-5.677862, -2.481141) -> Cubic
    static constexpr std::array<double, 4> SEG2_COEFFS = { 1.23478003e+00, 8.90103893e+00, -2.02810020e+01, 2.27017396e+02 };

    // Segment 3: Vg in [-2.481141, 1.722849) -> Quadratic
    static constexpr std::array<double, 3> SEG3_COEFFS = { -2.89953394e-01, -4.30851555e+01, 2.08157285e+02 };

    // Segment 4: Vg in [1.722849, 6.116238) -> Quadratic
    static constexpr std::array<double, 3> SEG4_COEFFS = { 4.66301432e+00, -6.01515846e+01, 2.22858724e+02 };

    // Segment 5: Vg in [6.116238, 17.628460] -> Quadratic
    static constexpr std::array<double, 3> SEG5_COEFFS = { 1.66816510e-01, -5.15195138e+00, 5.46632961e+01 };
};

/**
 * @struct TubeStageV1B
 * @brief A data-only struct for the second 12AX7 triode model (V1B).
 *
 * This struct contains a different set of parameters but is fully compatible
 * with the generic TubeModel<T> processor class.
 */
struct TubeStageV1B
{
    // --- Clamping Voltages ---
    static constexpr double CLAMP_HIGH_V = 405.0;
    static constexpr double CLAMP_LOW_V  = 19.543140;

    // --- Optimized Segment Boundaries ---
    static constexpr double BOUNDARY_0 = -9.699668;
    static constexpr double BOUNDARY_1 = -5.672815;
    static constexpr double BOUNDARY_2 = -2.072053;
    static constexpr double BOUNDARY_3 = 2.602954;
    static constexpr double BOUNDARY_4 = 9.289066;
    static constexpr double BOUNDARY_5 = 21.894080;

    // --- Polynomial Coefficients ---

    // Segment 1: Vg in [-9.699668, -5.672815) -> Quadratic
    static constexpr std::array<double, 3> SEG1_COEFFS = { -9.30546152e-02, -1.68233448e+00, 3.97396323e+02 };

    // Segment 2: Vg in [-5.672815, -2.072053) -> Cubic
    static constexpr std::array<double, 4> SEG2_COEFFS = { 1.10799712e+00, 8.33020329e+00, -1.30839726e+01, 2.63920916e+02 };

    // Segment 3: Vg in [-2.072053, 2.602954) -> Quadratic
    static constexpr std::array<double, 3> SEG3_COEFFS = { -1.52942219e-01, -3.39677900e+01, 2.47213187e+02 };

    // Segment 4: Vg in [2.602954, 9.289066) -> Quadratic
    static constexpr std::array<double, 3> SEG4_COEFFS = { 2.39975190e+00, -4.72568819e+01, 2.64508636e+02 };

    // Segment 5: Vg in [9.289066, 21.894080] -> Quadratic
    static constexpr std::array<double, 3> SEG5_COEFFS = { 1.29938811e-01, -5.08799549e+00, 6.86538555e+01 };
};

/**
 * @struct TubeStageV3B
 * @brief A data-only struct for a third tube model (V3B).
 *
 * This struct contains a different set of parameters but is fully compatible
 * with the generic TubeModel<T> processor class.
 */
struct TubeStageV3B
{
    // --- Clamping Voltages ---
    static constexpr double CLAMP_HIGH_V = 410.0;
    static constexpr double CLAMP_LOW_V  = 22.570510;

    // --- Optimized Segment Boundaries ---
    static constexpr double BOUNDARY_0 = -9.887168;
    static constexpr double BOUNDARY_1 = -5.939045;
    static constexpr double BOUNDARY_2 = -1.566184;
    static constexpr double BOUNDARY_3 = 3.612823;
    static constexpr double BOUNDARY_4 = 12.800175;
    static constexpr double BOUNDARY_5 = 25.222780;

    // --- Polynomial Coefficients ---

    // Segment 1: Vg in [-9.887168, -5.939045) -> Quadratic
    static constexpr std::array<double, 3> SEG1_COEFFS = { -5.31602679e-02, -9.66479730e-01, 4.05607294e+02 };

    // Segment 2: Vg in [-5.939045, -1.566184) -> Cubic
    static constexpr std::array<double, 4> SEG2_COEFFS = { 4.11523852e-01, 1.57922933e+00, -2.51229307e+01, 2.90770443e+02 };

    // Segment 3: Vg in [-1.566184, 3.612823) -> Quadratic
    static constexpr std::array<double, 3> SEG3_COEFFS = { -3.54365215e-01, -2.81513399e+01, 2.89189405e+02 };

    // Segment 4: Vg in [3.612823, 12.800175) -> Quadratic
    static constexpr std::array<double, 3> SEG4_COEFFS = { 1.56190133e+00, -4.19976045e+01, 3.14201458e+02 };

    // Segment 5: Vg in [12.800175, 25.222780] -> Quadratic
    static constexpr std::array<double, 3> SEG5_COEFFS = { 9.74324848e-02, -4.50669058e+00, 7.42563358e+01 };
};

/**
 * @class TubeModel
 * @brief A generic, templated processor for a piecewise polynomial tube model.
 *
 * This class takes a struct (like TubeStageV1A) as a template parameter
 * and uses its static members to perform the processing. This separates the
 * algorithm from the data, allowing for flexible and efficient use with
 * different tube models.
 */
template <typename TubeParams>
class TubeModel
{
public:
    using Params = TubeParams;

    /**
     * @brief Processes a single input sample (Vg) to produce an output sample (Vp).
     * @param vg The input grid voltage.
     * @return The corresponding plate voltage based on the provided TubeParams.
     */
    static double processSample(double vg)
    {
        // 1. Handle the clamping regions
        if (vg < TubeParams::BOUNDARY_0) {
            return TubeParams::CLAMP_HIGH_V;
        }
        else if (vg >= TubeParams::BOUNDARY_5) {
            return TubeParams::CLAMP_LOW_V;
        }

        // 2. Select the correct polynomial based on the segment
        else if (vg < TubeParams::BOUNDARY_1) {
            return evalQuadratic(vg, TubeParams::SEG1_COEFFS);
        }
        else if (vg < TubeParams::BOUNDARY_2) {
            return evalCubic(vg, TubeParams::SEG2_COEFFS);
        }
        else if (vg < TubeParams::BOUNDARY_3) {
            return evalQuadratic(vg, TubeParams::SEG3_COEFFS);
        }
        else if (vg < TubeParams::BOUNDARY_4) {
            return evalQuadratic(vg, TubeParams::SEG4_COEFFS);
        }
        else {
            return evalQuadratic(vg, TubeParams::SEG5_COEFFS);
        }
    }

private:
    // --- Efficient Polynomial Evaluation using Horner's Method ---

    /**
     * @brief Evaluates ax^2 + bx + c
     */
    static inline double evalQuadratic(double x, const std::array<double, 3>& c) {
        return (c[0] * x + c[1]) * x + c[2];
    }

    /**
     * @brief Evaluates ax^3 + bx^2 + cx + d
     */
    static inline double evalCubic(double x, const std::array<double, 4>& c) {
        return ((c[0] * x + c[1]) * x + c[2]) * x + c[3];
    }
};
