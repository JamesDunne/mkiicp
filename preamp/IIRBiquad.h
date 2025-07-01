#pragma once

// A standard Direct Form 1 IIR Biquad filter.
class IIRBiquad {
public:
    IIRBiquad() {
        reset();
    }

    void reset() {
        a1 = a2 = 0.0;
        b0 = 1.0; b1 = b2 = 0.0;
        x1 = x2 = y1 = y2 = 0.0;
    }

    // Set coefficients for the transfer function:
    // H(z) = (b0 + b1*z^-1 + b2*z^-2) / (1 + a1*z^-1 + a2*z^-2)
    void setCoefficients(double c_b0, double c_b1, double c_b2, double c_a0, double c_a1, double c_a2) {
        double a0_recip = 1.0 / c_a0;
        b0 = c_b0 * a0_recip;
        b1 = c_b1 * a0_recip;
        b2 = c_b2 * a0_recip;
        a1 = c_a1 * a0_recip;
        a2 = c_a2 * a0_recip;
    }

    // Prime the filter's state for a given DC value to prevent startup pops.
    void prime(double value) {
        // For a constant input, the output should be the same constant after settling.
        // We can solve for the steady-state output y = x * (b_sum / a_sum).
        double dc_gain = (b0 + b1 + b2) / (1.0 + a1 + a2);
        x1 = x2 = value;
        y1 = y2 = value * dc_gain;
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
    double a1, a2, b0, b1, b2;
    double x1, x2, y1, y2; // state
};
