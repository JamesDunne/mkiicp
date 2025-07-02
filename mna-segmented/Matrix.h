#pragma once
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
