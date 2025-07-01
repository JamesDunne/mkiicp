#pragma once

#include <cmath>

// A standard Direct Form 1 IIR Biquad filter.
class IIRBiquad {
public:
    IIRBiquad() {
        reset();
    }

    void reset() {
        a0 = 1.0; a1 = 0.0; a2 = 0.0;
        b0 = 1.0; b1 = 0.0; b2 = 0.0;
        x1 = x2 = y1 = y2 = 0.0;
    }

    // Set coefficients for the transfer function:
    // H(z) = (b0 + b1*z^-1 + b2*z^-2) / (a0 + a1*z^-1 + a2*z^-2)
    void setCoefficients(double c_b0, double c_b1, double c_b2, double c_a0, double c_a1, double c_a2) {
        b0 = c_b0 / c_a0;
        b1 = c_b1 / c_a0;
        b2 = c_b2 / c_a0;
        a1 = c_a1 / c_a0;
        a2 = c_a2 / c_a0;
    }

    // Process a single sample
    inline double process(double in) {
        double out = b0 * in + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
        x2 = x1;
        x1 = in;
        y2 = y1;
        y1 = out;
        return out;
    }

private:
    double a0, a1, a2, b0, b1, b2;
    double x1, x2, y1, y2; // state
};
