#include <iostream>
#include <vector>
#include <cmath> // For std::pow in the test, not in the optimized model

class TubeAmpModel {
public:
    // Constructor
    TubeAmpModel() = default;

    // Process a single sample of the grid voltage to get the plate voltage.
    // This is designed for real-time audio processing.
    double process(double gridVoltage) const {
        // Region 1: Cutoff (tube is off)
        if (gridVoltage < B1_CUTOFF) {
            return 405.0;
        }
        // Region 2: First polynomial segment (approaching linear region)
        else if (gridVoltage < B2_INFLECTION1) {
            return evalPoly(gridVoltage, P1_A, P1_B, P1_C, P1_D, P1_E, P1_F);
        }
        // Region 3: Second polynomial segment (most linear region)
        else if (gridVoltage < B3_INFLECTION2) {
            return evalPoly(gridVoltage, P2_A, P2_B, P2_C, P2_D, P2_E, P2_F);
        }
        // Region 4: Third polynomial segment (approaching saturation)
        else if (gridVoltage < B4_SATURATION) {
            return evalPoly(gridVoltage, P3_A, P3_B, P3_C, P3_D, P3_E, P3_F);
        }
        // Region 5: Saturation (grid conducting heavily)
        else {
            return 16.4;
        }
    }

private:
    // --- Boundary Points ---
    // These constants define the grid voltage ranges for each polynomial.
    static constexpr double B1_CUTOFF      = -8.0;
    static constexpr double B2_INFLECTION1 = -4.38;
    static constexpr double B3_INFLECTION2 = -0.73;
    static constexpr double B4_SATURATION  = 0.96;

    // --- Polynomial 1 Coefficients: Vg in [-8.0, -4.38] ---
    static constexpr double P1_A = -0.063718;
    static constexpr double P1_B = -1.542918;
    static constexpr double P1_C = -13.849313;
    static constexpr double P1_D = -55.882995;
    static constexpr double P1_E = -96.393165;
    static constexpr double P1_F = -1.916161;

    // --- Polynomial 2 Coefficients: Vg in [-4.38, -0.73] ---
    static constexpr double P2_A = 0.536341;
    static constexpr double P2_B = 3.651717;
    static constexpr double P2_C = 3.591434;
    static constexpr double P2_D = -19.462319;
    static constexpr double P2_E = -82.611183;
    static constexpr double P2_F = 175.717978;

    // --- Polynomial 3 Coefficients: Vg in [-0.73, 0.96] ---
    static constexpr double P3_A = 8.139139;
    static constexpr double P3_B = -4.321855;
    static constexpr double P3_C = -21.053424;
    static constexpr double P3_D = -1.419223;
    static constexpr double P3_E = -83.693156;
    static constexpr double P3_F = 171.328325;

    // Efficiently evaluate a 5th-degree polynomial using Horner's method.
    // P(x) = f + x*(e + x*(d + x*(c + x*(b + x*a))))
    inline double evalPoly(double x, double a, double b, double c, double d, double e, double f) const {
        return f + x * (e + x * (d + x * (c + x * (b + x * a))));
    }
};
