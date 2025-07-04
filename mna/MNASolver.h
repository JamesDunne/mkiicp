#pragma once

#include <vector>
#include <array>
#include <cmath>
#include <cstddef>
#include <numeric>

#if defined(__GNUC__) || defined(__clang__)
#define RESTRICT __restrict__
#elif defined(_MSC_VER)
#define RESTRICT __declspec(restrict)
#else
#define RESTRICT
#endif

inline static long long converged = 0, diverged = 0, failed = 0;

template <size_t Ns, size_t Vs>
class MNASolver {
public:
    static constexpr size_t NumUnknowns = Ns + Vs;
    static constexpr size_t NumNodes = Ns;

    MNASolver() : sampleRate(48000.0), T(1.0 / 48000.0), isDirty(true) {
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
        isDirty = true;
    }

protected:
    // MNA state
    std::array<std::array<double, NumUnknowns>, NumUnknowns> A;
    std::array<double, NumUnknowns> x; // Solution vector [node voltages..., vsrc currents...]
    std::array<double, NumUnknowns> b; // Sources vector

    // A matrix containing only the contributions from linear components (R, C conductance)
    std::array<std::array<double, NumUnknowns>, NumUnknowns> A_linear;
    std::array<double, NumUnknowns> b_linear;
    bool isDirty;

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

    // Helper to stamp into the linear matrix
    void stampResistorLinear(int n1, int n2, double R) {
        double G = 1.0 / R;
        if (n1 != -1) A_linear[n1][n1] += G;
        if (n2 != -1) A_linear[n2][n2] += G;
        if (n1 != -1 && n2 != -1) {
            A_linear[n1][n2] -= G;
            A_linear[n2][n1] -= G;
        }
    }

    // Stamps the capacitor's conductance and its history-dependent current source.
    void stampCapacitor(int n1, int n2, double C, const double& z_state) {
        double Gc = invT_2 * C;

        // Stamp the conductance Gc = 2C/T
        if (n1 != -1) A[n1][n1] += Gc;
        if (n2 != -1) A[n2][n2] += Gc;
        if (n1 != -1 && n2 != -1) {
            A[n1][n2] -= Gc;
            A[n2][n1] -= Gc;
        }

        // Stamp the history current source I_eq = z_state
        // Note the direction: current flows from n2 to n1
        if (n1 != -1) b[n1] += z_state;
        if (n2 != -1) b[n2] -= z_state;
    }

    /**
     * @brief Stamps the LINEAR part (conductance) of a capacitor into A_linear.
     * This should be called from within a derived class's stampLinear() override.
     */
    void stampCapacitor_A(int n1, int n2, double C) {
        double Gc = invT_2 * C;
        if (n1 != -1) A_linear[n1][n1] += Gc;
        if (n2 != -1) A_linear[n2][n2] += Gc;
        if (n1 != -1 && n2 != -1) {
            A_linear[n1][n2] -= Gc;
            A_linear[n2][n1] -= Gc;
        }
    }

    /**
     * @brief Stamps the DYNAMIC part (history current source) of a capacitor into b.
     * This should be called from within a derived class's stampDynamic() override.
     */
    void stampCapacitor_b(int n1, int n2, const double& z_state) {
        // The z_state *is* the history current source I_eq.
        // Current flows from n2 to n1.
        if (n1 != -1) b[n1] += z_state;
        if (n2 != -1) b[n2] -= z_state;
    }

    // Updates the capacitor's history state 'z' for the next time step.
    void updateCapacitorState(double v_n1, double v_n2, double C, double& z_state) {
        double Gc = invT_2 * C;
        double vc = v_n1 - v_n2; // Current voltage across the capacitor
        z_state = 2.0 * Gc * vc - z_state;
    }

    void stampVoltageSource(int n_p, int n_n, int v_idx, double voltage) {
        int idx = Ns + v_idx;
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

    /**
     * @brief Stamps the LINEAR part (topology) of a voltage source into A_linear.
     * This should be called from within a derived class's stampLinear() override.
     */
    void stampVoltageSource_A(int n_p, int n_n, int v_idx) {
        int idx = NumNodes + v_idx;
        if (n_p != -1) {
            A_linear[n_p][idx] += 1.0;
            A_linear[idx][n_p] += 1.0;
        }
        if (n_n != -1) {
            A_linear[n_n][idx] -= 1.0;
            A_linear[idx][n_n] -= 1.0;
        }
    }

    /**
     * @brief Stamps the DYNAMIC part (voltage value) of a voltage source into the main b vector.
     * This should be called from within a derived class's stampDynamic() override.
     */
    void stampVoltageSource_b(int v_idx, double voltage) {
        b[NumNodes + v_idx] = voltage; // Use '=' instead of '+=' as it's the only source here
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

    // --- NEW: Optimized LU Decomposition for Vectorization ---
    bool lu_decompose() {
        A_lu = A;
        for (size_t i = 0; i < NumUnknowns; ++i) pivot[i] = i;

        for (size_t i = 0; i < NumUnknowns; ++i) {
            // Find pivot (unchanged)
            double max_val = 0.0;
            size_t max_row = i;
            for (size_t k = i; k < NumUnknowns; ++k) {
                if (std::abs(A_lu[k][i]) > max_val) {
                    max_val = std::abs(A_lu[k][i]);
                    max_row = k;
                }
            }
            if (max_val < 1e-12) return false;
            if (i != max_row) {
                std::swap(A_lu[i], A_lu[max_row]);
                std::swap(pivot[i], pivot[max_row]);
            }

            // --- Vectorized numerical loops ---
            // Invert the diagonal element once
            const double inv_diag = 1.0 / A_lu[i][i];

            // This outer loop is still serial
            for (size_t j = i + 1; j < NumUnknowns; ++j) {
                // Isolate the multiplier to help the compiler
                const double multiplier = A_lu[j][i] * inv_diag;
                A_lu[j][i] = multiplier;

                // Get raw pointers and apply the RESTRICT hint
                double* RESTRICT row_j = &A_lu[j][0];
                const double* RESTRICT row_i = &A_lu[i][0];

                // This inner loop is the primary target for SIMD vectorization.
                // It performs a vector operation: row_j -= multiplier * row_i
                for (size_t k = i + 1; k < NumUnknowns; ++k) {
                    row_j[k] -= multiplier * row_i[k];
                }
            }
        }
        return true;
    }

    void lu_solve(std::array<double, NumUnknowns>& result) {
        //lu_solve_residual(b, result);
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

    // This helper now takes the residual 'f(x)' and calculates 'dx = J_inv * f(x)'
    void lu_solve_residual(const std::array<double, NumUnknowns>& rhs, std::array<double, NumUnknowns>& result) {
        std::array<double, NumUnknowns> temp_b;
        for (size_t i = 0; i < NumUnknowns; ++i) { temp_b[i] = rhs[pivot[i]]; }
        for (size_t i = 0; i < NumUnknowns; ++i) {
            double sum = temp_b[i];
            for (size_t j = 0; j < i; ++j) { sum -= A_lu[i][j] * result[j]; }
            result[i] = sum;
        }
        for (int i = NumUnknowns - 1; i >= 0; --i) {
            double sum = result[i];
            for (size_t j = i + 1; j < NumUnknowns; ++j) { sum -= A_lu[i][j] * result[j]; }
            result[i] = sum / A_lu[i][i];
        }
    }

    void resetMatrices() {
        for(auto& row : A) row.fill(0.0);
        b.fill(0.0);
    }

    /** @brief Stamps constant components (resistors, fixed sources) into A_linear and b_linear. */
    virtual void stampLinear() {}

    /** @brief Stamps components that change each sample (inputs, capacitors) into A and b. */
    virtual void stampDynamic(double in) {}

    /** @brief Stamps non-linear components (tubes, diodes) into A and b. */
    virtual void stampNonLinear(const std::array<double, NumUnknowns>& current_x) {}

    void setDirty() { isDirty = true; }

    // --- NEW: Generic Non-Linear Solver ---
    void solveNonlinear(double in) {
        if (isDirty) {
            A_linear.fill({});
            b_linear.fill({});
            stampLinear();
            isDirty = false;
        }

        const int MAX_ITER = 25;
        const double CONVERGENCE_THRESH = 1e-6;
        const double DAMPING_LIMIT = 1.0;

        std::array<double, NumUnknowns> current_x = x;
        for (int i = 0; i < MAX_ITER; ++i) {
            A = A_linear;
            b = b_linear;

            stampDynamic(in);
            stampNonLinear(current_x);

            if (!lu_decompose()) {
                x.fill(0.0); // Recover from singular matrix
                return;
            }

            std::array<double, NumUnknowns> next_x;
            lu_solve(next_x);

            double max_delta = 0.0;
            for (size_t j = 0; j < NumNodes; ++j) {
                max_delta = std::max(max_delta, std::abs(next_x[j] - current_x[j]));
            }

            if (max_delta < CONVERGENCE_THRESH) {
                x = next_x;
                return; // Converged
            }

            if (max_delta > DAMPING_LIMIT) {
                double scale = DAMPING_LIMIT / max_delta;
                for (size_t j = 0; j < NumUnknowns; ++j) {
                    current_x[j] += scale * (next_x[j] - current_x[j]);
                }
            } else {
                current_x = next_x;
            }
        }
        x = current_x; // Store result even if not fully converged
    }

    // --- NEW: Simplified Newton Method for High Performance ---
    void solveNonlinear_Simplified(double in) {
        if (isDirty) {
            A_linear.fill({});
            b_linear.fill({});
            stampLinear();
            isDirty = false;
        }

        // --- Iteration Loop ---
        const int MAX_ITER = 40; // Allow more iterations since convergence is slower
        const double CONVERGENCE_THRESH = 1e-6;
        const double DAMPING_LIMIT = 1.0;

        int i;
        double max_delta = 0.0;
        std::array<double, NumUnknowns> current_x = x;
        for (i = 0; i < MAX_ITER; ++i) {
            // Re-build the right-hand side 'b' and the full matrix 'A' at each step.
            // This is still much faster than re-decomposing.
            A = A_linear;
            b = b_linear;
            stampDynamic(in);
            stampNonLinear(current_x);

            // Decompose the matrix. If this fails, we can't continue.
            if (((i & 7) == 0) && !lu_decompose()) {
                x.fill(0.0);
                failed++;
                return;
            }

            // Solve using the PRE-CALCULATED LU factorization.
            std::array<double, NumUnknowns> next_x;
            lu_solve(next_x);

            // --- Dampening and Convergence Check (unchanged) ---
            max_delta = 0.0;
            for (size_t j = 0; j < NumNodes; ++j) {
                max_delta = std::max(max_delta, std::abs(next_x[j] - current_x[j]));
            }

            if (max_delta < CONVERGENCE_THRESH) {
                x = next_x;
                converged++;
                return; // Converged
            }

            if (max_delta > DAMPING_LIMIT) {
                double scale = DAMPING_LIMIT / max_delta;
                for (size_t j = 0; j < NumUnknowns; ++j) {
                    current_x[j] += scale * (next_x[j] - current_x[j]);
                }
            } else {
                current_x = next_x;
            }
        }

        diverged++;
        std::cerr << "diverged max_delta = " << max_delta << std::endl;

        x = current_x;
    }

    void solveNonlinear_Adaptive(double in) {
        if (isDirty) {
            A_linear.fill({}); b_linear.fill({});
            stampLinear();
            isDirty = false;
        }

        const int MAX_ITER = 60;
        const double REL_TOL = 1e-6;
        const double ABS_TOL = 1e-9;
        const double DAMPING_LIMIT = 1.0;
        const double STALL_THRESHOLD = 0.9;
        // NEW: Proactively refresh the Jacobian every N iterations
        const int JACOBIAN_REFRESH_INTERVAL = 5;

        std::array<double, NumUnknowns> current_x = x;
        double last_norm = 1e9;
        bool recompute_jacobian = true;
        bool just_dampened = false;

        for (int i = 0; i < MAX_ITER; ++i) {
            // --- Step 1: Proactive & Reactive Jacobian Update ---
            // Recompute if: first run, stalled, just dampened, OR at the periodic interval.
            if (recompute_jacobian || just_dampened || (i > 0 && ((i % JACOBIAN_REFRESH_INTERVAL) == 0))) {
                A = A_linear;
                b = b_linear;
                stampDynamic(in);
                stampNonLinear(current_x);
                if (!lu_decompose()) {
                    failed++;
                    x.fill(0.0);
                    return;
                }
                recompute_jacobian = false;
                just_dampened = false;
            }

            // --- Step 2: Solve A*x_next = b ---
            // Always solve for the absolute next state using the current LU factorization.
            std::array<double, NumUnknowns> next_x;
            lu_solve(next_x);

            // --- Step 3: Check for Convergence ---
            double delta_norm = 0.0, x_norm = 0.0;
            for (size_t j = 0; j < NumNodes; ++j) {
                double delta = next_x[j] - current_x[j];
                delta_norm += delta * delta;
                x_norm += next_x[j] * next_x[j];
            }
            delta_norm = std::sqrt(delta_norm);
            x_norm = std::sqrt(x_norm);

            if (delta_norm < (x_norm * REL_TOL + ABS_TOL)) {
                x = current_x;
                converged++;
                return; // Converged!
            }

            if (delta_norm > last_norm * STALL_THRESHOLD && i > 1) {
                recompute_jacobian = true;
            }
            last_norm = delta_norm;

            // --- Step 5: Apply Dampened Update x_new = x - scale*dx ---
            double scale = 1.0;
            if (delta_norm > DAMPING_LIMIT) {
                scale = DAMPING_LIMIT / delta_norm;
                just_dampened = true;
            }
            for (size_t j = 0; j < NumUnknowns; ++j) {
                current_x[j] += scale * (next_x[j] - current_x[j]);
            }
        }
        diverged++;
        x = current_x;
    }

    virtual void stampNonLinear_b_only(const std::array<double, NumUnknowns>& current_x) {
        // Default implementation just calls the full stampNonLinear.
        // Derived classes can provide a more optimized version.
        stampNonLinear(current_x);
    }
};
