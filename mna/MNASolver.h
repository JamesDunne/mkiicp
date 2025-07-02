#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <cstddef>
#include <numeric>

template <size_t NumNodes, size_t NumVsrc>
class MNASolver {
public:
    static constexpr size_t NumUnknowns = NumNodes + NumVsrc;

    MNASolver() : sampleRate(48000.0), T(1.0 / 48000.0) {
        x.fill(0.0);
        A.fill({});
        b.fill(0.0);
        pivot.fill(0);
    }

    virtual ~MNASolver() = default;

    void setup(double sr) {
        sampleRate = sr;
        T = 1.0 / sampleRate;
        invT_2 = 2.0 / T;
        x.fill(0.0);
    }

protected:
    // MNA state
    std::array<std::array<double, NumUnknowns>, NumUnknowns> A;
    std::array<double, NumUnknowns> x; // Solution vector [node voltages..., vsrc currents...]
    std::array<double, NumUnknowns> b; // Sources vector

    // LU factorization state
    std::array<std::array<double, NumUnknowns>, NumUnknowns> A_lu;
    std::array<int, NumUnknowns> pivot;

    // Simulation parameters
    double sampleRate;
    double T; // Sample period
    double invT_2; // 2.0 / T

    // --- Component Stamping ---
    // Nodes n1 and n2 are 0-based indices. Use -1 for ground.

    void stampResistor(int n1, int n2, double R) {
        double G = 1.0 / R;
        if (n1 != -1) A[n1][n1] += G;
        if (n2 != -1) A[n2][n2] += G;
        if (n1 != -1 && n2 != -1) {
            A[n1][n2] -= G;
            A[n2][n1] -= G;
        }
    }

    void stampCapacitor(int n1, int n2, double C, double& v_hist) {
        double Gc = invT_2 * C;
        double Ic = Gc * v_hist;

        if (n1 != -1) {
            A[n1][n1] += Gc;
            b[n1] += Ic;
        }
        if (n2 != -1) {
            A[n2][n2] += Gc;
            b[n2] -= Ic;
        }
        if (n1 != -1 && n2 != -1) {
            A[n1][n2] -= Gc;
            A[n2][n1] -= Gc;
        }
    }

    // Note: This function operates on voltages, so ground is simply 0.0. No change needed.
    void updateCapacitorHistory(double v_n1, double v_n2, double& v_hist) {
        double vc = v_n1 - v_n2;
        v_hist = invT_2 * vc - v_hist;
    }

    void stampVoltageSource(int n_p, int n_n, int v_idx, double voltage) {
        int idx = NumNodes + v_idx;
        if (n_p != -1) {
            A[n_p][idx] += 1.0;
            A[idx][n_p] += 1.0;
        }
        if (n_n != -1) {
            A[n_n][idx] -= 1.0;
            A[idx][n_n] -= 1.0;
        }
        b[idx] += voltage;
    }

    void stampConductance(int n1, int n2, double g) {
        if (n1 != -1) A[n1][n1] += g;
        if (n2 != -1) A[n2][n2] += g;
        if (n1 != -1 && n2 != -1) {
            A[n1][n2] -= g;
            A[n2][n1] -= g;
        }
    }

    void stampCurrentSource(int n1, int n2, double i) {
        if (n1 != -1) b[n1] -= i;
        if (n2 != -1) b[n2] += i;
    }

    // --- Solver (No changes needed below this line) ---
    bool lu_decompose() {
        A_lu = A;
        for (size_t i = 0; i < NumUnknowns; ++i) pivot[i] = i;

        for (size_t i = 0; i < NumUnknowns; ++i) {
            double max_val = 0.0;
            size_t max_row = i;
            for (size_t k = i; k < NumUnknowns; ++k) {
                if (std::abs(A_lu[k][i]) > max_val) {
                    max_val = std::abs(A_lu[k][i]);
                    max_row = k;
                }
            }

            if (max_val < 1e-12) return false;

            std::swap(A_lu[i], A_lu[max_row]);
            std::swap(pivot[i], pivot[max_row]);

            for (size_t j = i + 1; j < NumUnknowns; ++j) {
                A_lu[j][i] /= A_lu[i][i];
                for (size_t k = i + 1; k < NumUnknowns; ++k) {
                    A_lu[j][k] -= A_lu[j][i] * A_lu[i][k];
                }
            }
        }
        return true;
    }

    void lu_solve(std::array<double, NumUnknowns>& result) {
        std::array<double, NumUnknowns> temp_b;
        for (size_t i = 0; i < NumUnknowns; ++i) temp_b[i] = b[pivot[i]];

        for (size_t i = 0; i < NumUnknowns; ++i) {
            result[i] = temp_b[i];
            for (size_t j = 0; j < i; ++j) {
                result[i] -= A_lu[i][j] * result[j];
            }
        }

        for (int i = NumUnknowns - 1; i >= 0; --i) {
            for (size_t j = i + 1; j < NumUnknowns; ++j) {
                result[i] -= A_lu[i][j] * result[j];
            }
            result[i] /= A_lu[i][i];
        }
    }

    void resetMatrices() {
        for(auto& row : A) row.fill(0.0);
        b.fill(0.0);
    }
};
