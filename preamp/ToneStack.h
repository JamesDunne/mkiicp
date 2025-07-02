#pragma once

#ifndef USE_TONESTACK_MNA

#include "IIRBiquad.h"

class ToneStack {
public:
    ToneStack();
    void prepare(double sampleRate, double v1a_Vp_dc);
    void reset();
    void setParams(double treble, double mid, double bass, double volume);
    double process(double in);

private:
    void calculateCoefficients();

    double sampleRate;
    double p_treble, p_mid, p_bass, p_vol;

    IIRBiquad filter;
    IIRBiquad dcBlocker;
    double vol_gain;
};

#else

#include <vector>
#include <cmath>
#include <stdexcept>
#include <numeric>

// A simple fixed-size vector class
template<int N>
class Vector {
public:
    double data[N] = {0.0};

    double& operator[](int i) { return data[i]; }
    const double& operator[](int i) const { return data[i]; }
};

// A simple fixed-size square matrix class
template<int N>
class Matrix {
public:
    double data[N][N] = {{0.0}};

    double* operator[](int i) { return data[i]; }
    const double* operator[](int i) const { return data[i]; }

    // LU Decomposition based Inversion
    static Matrix<N> invert(Matrix<N> A) {
        Matrix<N> I; // Identity matrix
        for (int i = 0; i < N; ++i) I[i][i] = 1.0;

        for (int i = 0; i < N; ++i) {
            // Pivoting
            int pivot = i;
            for (int j = i + 1; j < N; ++j) {
                if (std::abs(A[j][i]) > std::abs(A[pivot][i])) {
                    pivot = j;
                }
            }
            if (pivot != i) {
                for(int k = 0; k < N; ++k) {
                    std::swap(A[i][k], A[pivot][k]);
                    std::swap(I[i][k], I[pivot][k]);
                }
            }

            if (std::abs(A[i][i]) < 1e-12) {
                 // This should not happen in a well-defined passive circuit
                throw std::runtime_error("Matrix is singular and cannot be inverted.");
            }

            // Normalize row
            double div = A[i][i];
            for (int j = 0; j < N; ++j) {
                A[i][j] /= div;
                I[i][j] /= div;
            }

            // Eliminate column
            for (int j = 0; j < N; ++j) {
                if (i != j) {
                    double mult = A[j][i];
                    for (int k = 0; k < N; ++k) {
                        A[j][k] -= mult * A[i][k];
                        I[j][k] -= mult * I[i][k];
                    }
                }
            }
        }
        return I;
    }
};

// Operator overloads for matrix/vector math
template<int N>
Vector<N> operator+(const Vector<N>& a, const Vector<N>& b) {
    Vector<N> result;
    for (int i = 0; i < N; ++i) result[i] = a[i] + b[i];
    return result;
}

template<int N>
Vector<N> operator-(const Vector<N>& a, const Vector<N>& b) {
    Vector<N> result;
    for (int i = 0; i < N; ++i) result[i] = a[i] - b[i];
    return result;
}

template<int N>
Vector<N> operator*(double s, const Vector<N>& v) {
    Vector<N> result;
    for (int i = 0; i < N; ++i) result[i] = s * v[i];
    return result;
}

template<int N>
Matrix<N> operator+(const Matrix<N>& A, const Matrix<N>& B) {
    Matrix<N> result;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            result[i][j] = A[i][j] + B[i][j];
    return result;
}

template<int N>
Matrix<N> operator-(const Matrix<N>& A, const Matrix<N>& B) {
    Matrix<N> result;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            result[i][j] = A[i][j] - B[i][j];
    return result;
}

template<int N>
Matrix<N> operator*(double s, const Matrix<N>& A) {
    Matrix<N> result;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            result[i][j] = s * A[i][j];
    return result;
}

template<int N>
Vector<N> operator*(const Matrix<N>& A, const Vector<N>& v) {
    Vector<N> result;
    for (int i = 0; i < N; ++i) {
        result[i] = 0;
        for (int j = 0; j < N; ++j) {
            result[i] += A[i][j] * v[j];
        }
    }
    return result;
}

// The number of unknown nodes in our circuit analysis
constexpr int NUM_NODES = 7;

class ToneStackFilter {
public:
    // Constructor: initializes the filter with a default sample rate
    ToneStackFilter();

    // Sets the sample rate. Call this if the audio device changes rate.
    void setSampleRate(double newSampleRate);

    // Sets the tone and volume controls (values from 0.0 to 1.0)
    void setParams(double treble, double mid, double bass, double vol1);

    // Processes a single audio sample and returns the filtered output
    double process(double inputSample);

private:
    // Re-calculates internal matrices when parameters or sample rate change
    void updateCoefficients();

    double sampleRate;
    double T; // Sample period (1 / sampleRate)

    // User-controlled parameters
    double p_treble, p_mid, p_bass, p_vol1;

    // State-space matrices and vectors for the discrete-time simulation
    Matrix<NUM_NODES> A_inv;
    Matrix<NUM_NODES> B_mat;
    Vector<NUM_NODES> C_vec;
    Vector<NUM_NODES> D_vec;

    // State variables (previous values)
    Vector<NUM_NODES> v_prev; // Previous node voltages v[n-1]
    double u_prev;            // Previous input sample u[n-1]
};

using ToneStack = ToneStackFilter;

#endif
