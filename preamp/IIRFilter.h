#pragma once

#include <cstring> // for memset

template <int ORDER>
class DirectFormIIR {
public:
    DirectFormIIR() {
        reset();
    }

    void reset() {
        memset(x, 0, sizeof(x));
        memset(y, 0, sizeof(y));
    }

    // H(z) = (b0 + b1*z^-1 + ... + bN*z^-N) / (a0 + a1*z^-1 + ... + aN*z^-N)
    void setCoefficients(const double b_coeffs[ORDER + 1], const double a_coeffs[ORDER + 1]) {
        double a0_recip = 1.0 / a_coeffs[0];
        for (int i = 0; i <= ORDER; ++i) {
            b[i] = b_coeffs[i] * a0_recip;
            a[i] = a_coeffs[i] * a0_recip;
        }
    }

    inline double process(double in) {
        // Shift delay lines
        for (int i = ORDER; i > 0; --i) {
            x[i] = x[i - 1];
            y[i] = y[i - 1];
        }
        x[0] = in;

        // Calculate output
        double out = b[0] * x[0];
        for (int i = 1; i <= ORDER; ++i) {
            out += b[i] * x[i] - a[i] * y[i];
        }

        y[0] = out;
        return out;
    }

private:
    double b[ORDER + 1]; // Numerator coefficients
    double a[ORDER + 1]; // Denominator coefficients
    double x[ORDER + 1]; // Input delay line
    double y[ORDER + 1]; // Output delay line
};
