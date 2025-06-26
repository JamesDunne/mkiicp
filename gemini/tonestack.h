#pragma once

#include <array>
#include <vector>

// Define fixed-size matrix and vector types for our 7-node system.
// Using std::array provides compile-time bounds checking and avoids dynamic allocation.
constexpr int MATRIX_SIZE = 7;
using Vector7 = std::array<double, MATRIX_SIZE>;
using Matrix7x7 = std::array<std::array<double, MATRIX_SIZE>, MATRIX_SIZE>;

class ToneStack {
public:
    /**
     * @brief Constructs the tone stack filter.
     * @param sampleRate The audio sample rate in Hz.
     */
    ToneStack(double sampleRate);

    /**
     * @brief Sets the audio sample rate.
     * @param sr The new sample rate in Hz.
     */
    void setSampleRate(double sr);

    /**
     * @brief Resets the internal state of the filter (clears memory).
     */
    void reset();

    /**
     * @brief Sets the control parameters for the tone stack.
     * @param treble Control (0.0 to 1.0).
     * @param bass Control (0.0 to 1.0).
     * @param mid Control (0.0 to 1.0).
     * @param volume Control (0.0 to 1.0).
     */
    void setParams(double treble, double bass, double mid, double volume);

    /**
     * @brief Processes a single audio sample.
     * @param input The input sample. The model expects the AC component of the signal.
     * @return The processed output sample, scaled to a nominal -1 to +1 range.
     */
    double process(double input);

private:
    // Parameters
    double sampleRate;
    double g; // 2.0 * sampleRate, for Bilinear Transform

    // Control parameters (0.0 to 1.0)
    double p_treble, p_bass, p_mid, p_volume;

    // State variables
    Vector7 V_nm1; // Previous node voltages [v5, v7, v8, v15, v16, v20, v24]
    double u_nm1;  // Previous input sample

    // DC Blocker state for input conditioning
    double dc_x1 = 0.0;
    double dc_y1 = 0.0;

    // Matrices for the state-space equation
    Matrix7x7 A_inv;
    Matrix7x7 B;
    Vector7 C_u;
    Vector7 D_u;

    bool params_changed;

    // Recalculates the matrices when parameters change.
    void updateCoefficients();

    /**
     * @brief Applies a standard audio taper (logarithmic feel) to a linear control.
     * @param x Linear value from 0.0 to 1.0.
     * @return Tapered value.
     */
    inline double audioTaper(double x) {
        // An x^2 curve is a common and effective approximation for an audio taper.
        return x * x;
    }
};

#include <cmath>
#include <algorithm> // For std::clamp
#include <stdexcept> // For std::runtime_error

// A small number to prevent division by zero in conductance calculations
constexpr double epsilon = 1e-12;

// Anonymous namespace for internal linear algebra helper functions
namespace {

// Inverts a 7x7 matrix using Gauss-Jordan elimination with partial pivoting.
// Returns false if the matrix is singular (non-invertible).
bool invertMatrix(const Matrix7x7& in, Matrix7x7& out) {
    Matrix7x7 a = in;

    // Initialize 'out' as the identity matrix
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            out[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int i = 0; i < MATRIX_SIZE; ++i) {
        // Partial Pivoting: Find the row with the largest pivot
        int pivot_row = i;
        for (int j = i + 1; j < MATRIX_SIZE; ++j) {
            if (std::abs(a[j][i]) > std::abs(a[pivot_row][i])) {
                pivot_row = j;
            }
        }

        // Swap rows if necessary
        if (pivot_row != i) {
            std::swap(a[i], a[pivot_row]);
            std::swap(out[i], out[pivot_row]);
        }

        double pivot = a[i][i];
        if (std::abs(pivot) < epsilon) {
            return false; // Matrix is singular or nearly singular
        }

        // Normalize the pivot row
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            a[i][j] /= pivot;
            out[i][j] /= pivot;
        }

        // Eliminate other entries in the current column
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            if (i != j) {
                double factor = a[j][i];
                for (int k = 0; k < MATRIX_SIZE; ++k) {
                    a[j][k] -= factor * a[i][k];
                    out[j][k] -= factor * out[i][k];
                }
            }
        }
    }
    return true;
}

// Matrix-Vector multiplication
Vector7 multiplyMatrixVector(const Matrix7x7& m, const Vector7& v) {
    Vector7 result{};
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            result[i] += m[i][j] * v[j];
        }
    }
    return result;
}

// Vector-scalar multiplication
Vector7 multiplyVectorScalar(const Vector7& v, double s) {
    Vector7 result{};
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        result[i] = v[i] * s;
    }
    return result;
}

// Vector-vector addition
Vector7 addVectors(const Vector7& v1, const Vector7& v2) {
    Vector7 result{};
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        result[i] = v1[i] + v2[i];
    }
    return result;
}

// Vector-vector subtraction
Vector7 subtractVectors(const Vector7& v1, const Vector7& v2) {
    Vector7 result{};
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        result[i] = v1[i] - v2[i];
    }
    return result;
}

// Matrix-matrix addition
Matrix7x7 addMatrices(const Matrix7x7& m1, const Matrix7x7& m2) {
    Matrix7x7 result{};
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            result[i][j] = m1[i][j] + m2[i][j];
        }
    }
    return result;
}

// Matrix-matrix subtraction
Matrix7x7 subtractMatrices(const Matrix7x7& m1, const Matrix7x7& m2) {
    Matrix7x7 result{};
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            result[i][j] = m1[i][j] - m2[i][j];
        }
    }
    return result;
}

// Matrix-scalar multiplication
Matrix7x7 multiplyMatrixScalar(const Matrix7x7& m, double s) {
    Matrix7x7 result{};
    for (int i = 0; i < MATRIX_SIZE; ++i) {
        for (int j = 0; j < MATRIX_SIZE; ++j) {
            result[i][j] = m[i][j] * s;
        }
    }
    return result;
}

} // end anonymous namespace

inline ToneStack::ToneStack(double sr) {
    reset();
    params_changed = true; // Force initial calculation
    setSampleRate(sr);
    setParams(0.5, 0.5, 0.5, 0.5); // Default to mid settings
}

inline void ToneStack::setSampleRate(double sr) {
    sampleRate = sr;
    g = 2.0 * sampleRate;
    params_changed = true;
}

inline void ToneStack::reset() {
    V_nm1.fill(0.0);
    u_nm1 = 0.0;
    dc_x1 = 0.0;
    dc_y1 = 0.0;
}

inline void ToneStack::setParams(double treble, double bass, double mid, double volume) {
    p_treble = std::clamp(treble, 0.0, 1.0);
    p_bass = std::clamp(bass, 0.0, 1.0);
    p_mid = std::clamp(mid, 0.0, 1.0);
    p_volume = std::clamp(volume, 0.0, 1.0);
    params_changed = true;
}

inline void ToneStack::updateCoefficients() {
    if (!params_changed) {
        return;
    }

    // Apply Audio Taper to controls
    double treble_tapered = audioTaper(p_treble);
    double bass_tapered = audioTaper(p_bass);
    double mid_tapered = audioTaper(p_mid);
    double volume_tapered = audioTaper(p_volume);

    // Define Component Values from Netlist
    double R5 = 100e3;
    double R5A = 100e3;

    double R_treble_C = 250e3 * treble_tapered + epsilon;
    double R_treble_A = 250e3 * (1.0 - treble_tapered) + epsilon;
    double R_bass = 250e3 * (1.0 - bass_tapered) + epsilon;
    double R_mid = 10e3 * (1.0 - mid_tapered) + epsilon;
    double R_vol_C = 1e6 * volume_tapered + epsilon;
    double R_vol_A = 1e6 * (1.0 - volume_tapered) + epsilon;

    double C_comb = 1e-9;
    double C4 = 0.1e-6;
    double C3 = 0.047e-6;
    double C13B = 180e-12;

    // Calculate Conductances
    double g5 = 1.0 / R5;
    double g5a = 1.0 / R5A;
    double g_tA = 1.0 / R_treble_A;
    double g_tC = 1.0 / R_treble_C;
    double g_b = 1.0 / R_bass;
    double g_m = 1.0 / R_mid;
    double g_vA = 1.0 / R_vol_A;
    double g_vC = 1.0 / R_vol_C;

    // Build G (Conductance) and C (Capacitance) Matrices
    // Based on KCL equations for each of the 7 unknown nodes:
    // v = [v5, v7, v8, v15, v16, v20, v24]^T
    Matrix7x7 G{}, C{}; // Zero-initialized

    // Node 5 (idx 0)
    G[0][0] = g_tA; G[0][1] = -g_tA;
    C[0][0] = C_comb;
    // Node 7 (idx 1)
    G[1][0] = -g_tA; G[1][1] = g_tA + g_tC + g5a; G[1][2] = -g5a; G[1][4] = -g_tC;
    // Node 8 (idx 2)
    G[2][1] = -g5a; G[2][2] = g5a + g_vA; G[2][5] = -g_vA;
    C[2][2] = C13B; C[2][5] = -C13B;
    // Node 15 (idx 3)
    G[3][3] = g5;
    C[3][3] = C4 + C3; C[3][4] = -C4; C[3][6] = -C3;
    // Node 16 (idx 4)
    G[4][1] = -g_tC; G[4][4] = g_tC + g_b; G[4][6] = -g_b;
    C[4][3] = -C4; C[4][4] = C4;
    // Node 20 (Output, idx 5)
    G[5][2] = -g_vA; G[5][5] = g_vA + g_vC;
    C[5][2] = -C13B; C[5][5] = C13B;
    // Node 24 (idx 6)
    G[6][4] = -g_b; G[6][6] = g_b + g_m;
    C[6][3] = -C3; C[6][6] = C3;

    // Build Source-related Vectors
    Vector7 D_G{}, D_C{}; // Zero-initialized
    D_G[3] = g5;
    D_C[0] = C_comb;

    // Form Discrete-Time Matrices using Bilinear Transform
    Matrix7x7 A_system = addMatrices(G, multiplyMatrixScalar(C, g));
    B = subtractMatrices(G, multiplyMatrixScalar(C, g));
    C_u = addVectors(D_G, multiplyVectorScalar(D_C, g));
    D_u = subtractVectors(D_G, multiplyVectorScalar(D_C, -g)); // D_u = D_G - g*D_C

    if (!invertMatrix(A_system, A_inv)) {
        // This should not happen with this circuit, but as a safeguard:
        // throw std::runtime_error("ToneStack matrix is singular. Cannot compute filter coefficients.");
        // In a real-time context, you might instead hold the previous coefficients.
        // For now, we'll just skip the update.
        return;
    }

    params_changed = false;
}

inline double ToneStack::process(double input) {
    if (params_changed) {
        updateCoefficients();
    }

    // Input Conditioning: Remove DC offset with a 1-pole HPF (~10 Hz)
    constexpr double R = 0.9997;
    double ac_input = input - dc_x1 + R * dc_y1;
    dc_x1 = input;
    dc_y1 = ac_input;

    // State-Variable Update Equation: V[n] = A⁻¹ * (-B * V[n-1] + Cᵤ * u[n] + Dᵤ * u[n-1])
    Vector7 term1 = multiplyMatrixVector(B, V_nm1);
    Vector7 term2 = multiplyVectorScalar(C_u, ac_input);
    Vector7 term3 = multiplyVectorScalar(D_u, u_nm1);

    Vector7 V_n = subtractVectors(addVectors(term2, term3), term1); // (term2 + term3) - term1
    V_n = multiplyMatrixVector(A_inv, V_n);

    // Update State for Next Sample
    V_nm1 = V_n;
    u_nm1 = ac_input;

    // Output is voltage at node N020, the 6th element (index 5) of our state vector
    double output = V_n[5];

    if (!std::isfinite(output)) {
        reset();
        return 0.0;
    }

    return output;
}
