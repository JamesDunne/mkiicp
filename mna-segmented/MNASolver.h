#pragma once

#include "Matrix.h" // Assumes our previous Matrix.h utility
#include <iostream>
#include <vector>

// N = Number of unknown nodes in the MNA system
template <int N>
class MNASolver {
public:
    // Constructor: specifies sample rate and whether the circuit is non-linear
    MNASolver(double sampleRate, bool isNonLinear)
        : isNonLinear(isNonLinear), needsRecalculation(true) {
        setSampleRate(sampleRate);
        for(int i = 0; i < N; ++i) x_prev[i] = 0.0;
    }

    virtual ~MNASolver() = default;

    void setSampleRate(double newSampleRate) {
        if (newSampleRate > 0 && sampleRate != newSampleRate) {
            sampleRate = newSampleRate;
            g_s = 2.0 * sampleRate;
            needsRecalculation = true;
        }
    }

    // Call this from derived class when a parameter (pot, etc.) changes
    void flagForRecalculation() {
        needsRecalculation = true;
    }

    // The main processing function. Solves the system for one sample.
    // The input vector can be used for multi-input systems.
    double process(const std::vector<double>& inputs) {
        if (needsRecalculation) {
            recalculateMatrices();
            needsRecalculation = false;
        }

        if (isNonLinear) {
            // Solve using one iteration of Newton-Raphson
            // Linearized system is J * x_n = b
            
            // Start with the linear part of the system
            Matrix<N> J = G + g_s * C;
            Vector<N> b = (g_s * C) * x_prev;
            
            // Let the derived class add non-linear contributions to J and b
            stampNonlinear(J, b, x_prev, inputs);

            // Solve for the current state x_n
            try {
                Matrix<N> J_inv = Matrix<N>::invert(J);
                x_prev = J_inv * b;
            } catch(const std::runtime_error& e) {
                 // Solver failed, hold previous state to prevent artifacts
                 // In a real product, you might log this error.
                std::cerr << e.what() << std::endl;
            }

        } else {
            // Solve using pre-calculated matrices for the linear case
            // x_n = A_inv * (C_vec*u_n + D_vec*u_n-1 - B_mat*x_n-1)
            Vector<N> input_contrib = {0.0};
            for(size_t i = 0; i < inputs.size() && i < input_vectors.size(); ++i) {
                input_contrib = input_contrib + (inputs[i] * input_vectors[i].first) + (prev_inputs[i] * input_vectors[i].second);
            }

            Vector<N> x_n = A_inv * (input_contrib - B_mat * x_prev);
            x_prev = x_n;
        }
        
        // Store current inputs for next sample (z^-1 term)
        prev_inputs = inputs;
        
        // Return voltage of the designated output node
        return x_prev[getOutputNodeIndex()];
    }

    // Convenience overload for single-input systems
    double process(double input) {
        return process(std::vector<double>{input});
    }


protected:
    // --- Pure virtual functions to be implemented by the derived circuit class ---

    // Stamps the linear components (R, C) onto the G and C matrices.
    // Called only when recalculation is needed.
    virtual void stampLinear(Matrix<N>& G, Matrix<N>& C) = 0;

    // For linear circuits: Stamps the input source connections.
    // Each pair is (Fg + g_s*Fc) and (Fg - g_s*Fc) for an input.
    virtual void stampInputs(std::vector<std::pair<Vector<N>, Vector<N>>>& inputs) = 0;

    // For non-linear circuits: Stamps non-linear elements and their derivatives
    // onto the Jacobian J and RHS vector b. Called every sample.
    virtual void stampNonlinear(Matrix<N>& J, Vector<N>& b, const Vector<N>& x_k, const std::vector<double>& u_n) = 0;

    // Tells the solver which node index is the output.
    virtual int getOutputNodeIndex() const = 0;


private:
    void recalculateMatrices() {
        G = Matrix<N>();
        C = Matrix<N>();
        stampLinear(G, C);

        if (!isNonLinear) {
            Matrix<N> A_mat = G + g_s * C;
            B_mat = G - g_s * C;
            
            try {
                A_inv = Matrix<N>::invert(A_mat);
            } catch(const std::runtime_error& e) {
                std::cerr << "MNA Warning: Matrix is singular. A_inv set to identity." << std::endl;
                A_inv = Matrix<N>(); // Reset to identity or safe matrix
                for(int i=0; i<N; ++i) A_inv[i][i] = 1.0;
            }
            
            input_vectors.clear();
            stampInputs(input_vectors);
            prev_inputs.assign(input_vectors.size(), 0.0);
        }
    }

    double sampleRate;
    double g_s; // 2.0 * sampleRate
    bool isNonLinear;
    bool needsRecalculation;

    // State
    Vector<N> x_prev; // Previous node voltages x[n-1]
    std::vector<double> prev_inputs;

    // System Matrices
    Matrix<N> G, C; // Linear part

    // Pre-calculated matrices for linear systems
    Matrix<N> A_inv;
    Matrix<N> B_mat;
    std::vector<std::pair<Vector<N>, Vector<N>>> input_vectors;
};
