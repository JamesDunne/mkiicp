#pragma once

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <algorithm> // For std::max/min

/**
 * @class RealtimeTubeSim
 * @brief A single-header C++ class for realtime audio simulation of a tube amplifier circuit.
 *
 * This class simulates a single-stage 12AX7 triode preamp with a tone stack and master volume.
 * It is based on the Modified Nodal Analysis (MNA) method for circuit simulation.
 *
 * --- Core Simulation Method ---
 * 1.  **SPICE Netlist to MNA**: The circuit topology is hard-coded into the `build_system()` method,
 *     which populates the MNA matrices (G, C, b) based on the components.
 * 2.  **Time-Domain Simulation**: The simulation uses the Backward Euler method for discretizing
 *     time-dependent components (capacitors), which is a stable implicit method suitable for
 *     stiff systems like audio circuits. The resulting system equation is:
 *     (G + C/dt) * x_n = b_n + (C/dt) * x_{n-1}
 * 3.  **Non-Linear Solver**: The triode is a non-linear component. A Newton-Raphson (NR)
 *     iterative solver is used at each time step to find the node voltages. The NR method
 *     solves J * dx = -f(x), where J is the Jacobian matrix and f(x) is the residual.
 * 4.  **Matrix Solver**: At each NR iteration, a linear system of equations must be solved.
 *     This class implements a dense LU decomposition solver with pivoting for this purpose.
 *
 * --- Simulated Circuit (Conceptual SPICE Netlist) ---
 *
 * // Nodes: 0 (GND), in, g, p, k, V_P, ts1, ts2, ts3, out, out_final
 *
 * // --- Input Stage & Volume ---
 * V_in   in 0 DC 0.0      ; Input signal, represented by b vector update. 'Volume' is a digital gain on this.
 * R_in   in 0 1M          ; Input impedance (for realism, though V_in source dominates)
 * C_in   in g 22n         ; Input coupling capacitor
 * R_g    g  0 1M          ; Grid leak resistor
 *
 * // --- Triode Stage (12AX7/ECC83) ---
 * X_T1   p g k Triode12AX7 ; Plate, Grid, Cathode
 * R_p    V_P p 100k        ; Plate load resistor
 * R_k    k 0 1.5k          ; Cathode bias resistor
 * C_k    k 0 22u          ; Cathode bypass capacitor
 * V_P    V_P 0 DC 300.0    ; Plate supply voltage
 *
 * // --- Tone Stack (Fender-style) ---
 * C_ts_in p ts1 22n       ; Coupling cap to tone stack
 * R_t    ts1 ts2 250k     ; Treble Potentiometer (Pot1)
 * R_b    ts2 0   1M       ; Bass Potentiometer (Pot2)
 * C_t    ts1 ts3 250p     ; Treble capacitor
 * R_m    ts3 0   25k      ; Mid Potentiometer (Pot3)
 * C_m    ts2 ts3 22n      ; Mid capacitor
 *
 * // --- Master Volume & Output ---
 * R_master ts3 out 1M     ; Master Volume Potentiometer (Pot4)
 * C_out  out out_final 100n   ; Output coupling capacitor
 * R_load out_final 0 1M       ; Final load resistor
 *
 * // Fixed DC voltages (VE, VC, VC2 are conceptual placeholders, V_P is used)
 */
class RealtimeTubeSim {
public:
    // For readability
    using Matrix = std::vector<std::vector<double>>;
    using Vector = std::vector<double>;

    /**
     * @brief Constructor.
     * @param sample_rate The audio sample rate in Hz.
     */
    explicit RealtimeTubeSim(double sample_rate) {
        if (sample_rate <= 0) {
            throw std::invalid_argument("Sample rate must be positive.");
        }
        m_sample_rate = sample_rate;
        m_dt = 1.0 / sample_rate;

        build_system();

        // Initialize state vectors
        m_x.assign(m_num_nodes, 0.0);
        m_x_prev.assign(m_num_nodes, 0.0);

        // Set initial potentiometer values
        set_volume(0.5);
        set_bass(0.5);
        set_mid(0.5);
        set_treble(0.5);
        set_master(1.0);

        // Initial DC analysis to find the quiescent operating point
        solve_dc();
    }

    /**
     * @brief Sets the input 'volume' control (gain before the tube).
     * @param value A value from 0.0 to 1.0. This is modeled as a simple gain
     *              multiplier on the input signal for simplicity.
     */
    void set_volume(double value) { m_params.volume = std::max(0.0, std::min(1.0, value)); }

    void set_bass(double value) { m_params.bass = std::max(0.0, std::min(1.0, value)); update_pots(); }
    void set_mid(double value) { m_params.mid = std::max(0.0, std::min(1.0, value)); update_pots(); }
    void set_treble(double value) { m_params.treble = std::max(0.0, std::min(1.0, value)); update_pots(); }
    void set_master(double value) { m_params.master = std::max(0.0, std::min(1.0, value)); update_pots(); }

    /**
     * @brief Processes a single audio sample.
     * @param input_sample The input audio sample.
     * @return The processed output audio sample.
     */
    double process_sample(double input_sample) {
        // --- 1. Update time-varying sources ---
        m_b[m_node_map["in"]] = input_sample * m_params.volume;

        // --- 2. Newton-Raphson Iteration for this time step ---
        // Initial guess for this timestep's solution is the previous one.
        m_x = m_x_prev;

        // Pre-calculate the linear parts of the system for this timestep.
        // The system is A*x = z, where A = G + C/dt and z = b + (C/dt)*x_prev
        Matrix A_linear = m_G_dynamic;
        Vector z = m_b;
        for (int r = 0; r < m_num_nodes; ++r) {
            for (int c = 0; c < m_num_nodes; ++c) {
                if (m_C(r, c) != 0.0) {
                    double c_val_dt = m_C(r, c) / m_dt;
                    A_linear[r][c] += c_val_dt;
                    z[r] += c_val_dt * m_x_prev[c];
                }
            }
        }

        for (int i = 0; i < m_nr_max_iter; ++i) {
            // At each NR step, we solve J*dx = -f(x)
            // where J is the Jacobian and f is the residual.

            // Start with the linear system parts.
            Matrix J = A_linear;
            Vector f = matrix_vector_mult(A_linear, m_x);
            for(size_t j = 0; j < f.size(); ++j) f[j] -= z[j];

            // Add non-linear contributions to J and f.
            add_triode_nr_stamps(J, f, m_x);

            // Now, f is the residual f(x_k). We solve J*dx = -f
            for(size_t j = 0; j < f.size(); ++j) {
                f[j] = -f[j];
            }

            // Solve the linear system J * dx = -f for dx
            Vector dx = lu_solve(J, lu_decompose(J), f);

            // Update the solution: x_{k+1} = x_k + dx
            double norm = 0.0;
            for (int j = 0; j < m_num_nodes; ++j) {
                m_x[j] += dx[j];
                norm += dx[j] * dx[j];
            }

            // Check for convergence
            if (sqrt(norm) < m_nr_tolerance) {
                break;
            }
        }

        // --- 3. Update state for next sample ---
        m_x_prev = m_x;

        // --- 4. Return output voltage ---
        return m_x[m_node_map["out_final"]];
    }


private:
    // --- Typedefs for matrix access convenience ---
    struct MatrixView {
        Matrix& mat;
        double& operator()(size_t r, size_t c) { return mat.at(r).at(c); }
        const double& operator()(size_t r, size_t c) const { return mat.at(r).at(c); }
    };

    // --- MNA System Variables ---
    Matrix m_G_static;
    Matrix m_G_dynamic;
    MatrixView m_C { m_C_mat };
    Matrix m_C_mat;
    Vector m_b;
    Vector m_x;
    Vector m_x_prev;
    int m_num_nodes = 0;
    std::map<std::string, int> m_node_map;

    double m_sample_rate;
    double m_dt;

    const int m_nr_max_iter = 15;
    const double m_nr_tolerance = 1e-6;

    struct Potentiometer {
        int n1, n2, n_wiper;
        double total_resistance;
    };
    Potentiometer m_pot_treble, m_pot_bass, m_pot_mid, m_pot_master;

    struct ControlParams {
        double volume, bass, mid, treble, master;
    } m_params;

    // --- Triode Model Parameters (for a 12AX7 / ECC83) ---
    const double TUBE_MU = 100.0;      // Amplification factor
    const double TUBE_EX = 1.35;       // Exponent
    const double TUBE_KG = 1060.0;     // Koren's gain factor
    const double TUBE_KP = 600.0;      // Plate resistance factor
    const double TUBE_K_GRID = 2000.0; // Grid current factor (increased for more effect)
    const double TUBE_V_GRID_THRESH = 0.5; // Grid current threshold voltage

    int m_p_node, m_g_node, m_k_node;

private:
    // --- Core Setup Methods ---

    void build_system() {
        int node_idx = 0;
        auto add_node = [&](const std::string& name) {
            m_node_map[name] = node_idx++;
        };
        add_node("in"); add_node("g"); add_node("p"); add_node("k");
        add_node("V_P"); add_node("ts1"); add_node("ts2"); add_node("ts3");
        add_node("out"); add_node("out_final");

        m_num_nodes = node_idx;

        m_G_static.assign(m_num_nodes, Vector(m_num_nodes, 0.0));
        m_C_mat.assign(m_num_nodes, Vector(m_num_nodes, 0.0));
        m_b.assign(m_num_nodes, 0.0);

        stamp_resistor(m_G_static, "in", "gnd", 1e6);
        stamp_resistor(m_G_static, "g", "gnd", 1e6);
        stamp_resistor(m_G_static, "V_P", "p", 100e3);
        stamp_resistor(m_G_static, "k", "gnd", 1.5e3);
        stamp_resistor(m_G_static, "out_final", "gnd", 1e6);

        stamp_capacitor("in", "g", 22e-9);
        stamp_capacitor("k", "gnd", 22e-6);
        stamp_capacitor("p", "ts1", 22e-9);
        stamp_capacitor("ts1", "ts3", 250e-12);
        stamp_capacitor("ts2", "ts3", 22e-9);
        stamp_capacitor("out", "out_final", 100e-9);

        stamp_voltage_source("V_P", "gnd", 300.0);

        m_pot_treble = {m_node("ts1"), m_node("ts2"), -1, 250e3};
        m_pot_bass = {m_node("ts2"), m_node("gnd"), -1, 1e6};
        m_pot_mid = {m_node("ts3"), m_node("gnd"), -1, 25e3};
        m_pot_master = {m_node("ts3"), m_node("gnd"), m_node("out"), 1e6};

        m_p_node = m_node("p");
        m_g_node = m_node("g");
        m_k_node = m_node("k");
    }

    void update_pots() {
        m_G_dynamic = m_G_static;
        // Audio taper approximation
        auto taper = [](double val) { return val * val; };

        double treble_res = m_pot_treble.total_resistance * taper(m_params.treble);
        stamp_resistor(m_G_dynamic, m_pot_treble.n1, m_pot_treble.n2, std::max(1.0, treble_res));

        // Reverse audio taper
        double bass_res = m_pot_bass.total_resistance * taper(1.0 - m_params.bass);
        stamp_resistor(m_G_dynamic, m_pot_bass.n1, m_pot_bass.n2, std::max(1.0, bass_res));

        // Linear taper
        double mid_res = m_pot_mid.total_resistance * m_params.mid;
        stamp_resistor(m_G_dynamic, m_pot_mid.n1, m_pot_mid.n2, std::max(1.0, mid_res));

        // Master Volume (voltage divider)
        double r1_val = m_pot_master.total_resistance * taper(m_params.master);
        double r2_val = m_pot_master.total_resistance - r1_val;
        stamp_resistor(m_G_dynamic, m_pot_master.n1, m_pot_master.n_wiper, std::max(1.0, r1_val));
        stamp_resistor(m_G_dynamic, m_pot_master.n_wiper, m_pot_master.n2, std::max(1.0, r2_val));
    }

    void solve_dc() {
        m_x.assign(m_num_nodes, 0.0);

        for (int i = 0; i < m_nr_max_iter * 2; ++i) {
            Matrix J = m_G_dynamic;
            Vector f = matrix_vector_mult(m_G_dynamic, m_x);
            for(size_t j = 0; j < f.size(); ++j) f[j] -= m_b[j];

            add_triode_nr_stamps(J, f, m_x);
            for(size_t j = 0; j < f.size(); ++j) f[j] = -f[j];

            Vector dx = lu_solve(J, lu_decompose(J), f);
            double norm = 0.0;
            for (int j = 0; j < m_num_nodes; ++j) {
                m_x[j] += dx[j];
                norm += dx[j] * dx[j];
            }
            if (sqrt(norm) < m_nr_tolerance) break;
        }
        m_x_prev = m_x;
    }

    // --- MNA Stamping Helpers ---
    int m_node(const std::string& name) {
        if (name == "gnd" || name == "0") return -1;
        return m_node_map.at(name);
    }

    void stamp_resistor(Matrix& G, int n1, int n2, double R) {
        if (R <= 1e-9) R = 1e-9;
        double g = 1.0 / R;
        if (n1 != -1) G[n1][n1] += g;
        if (n2 != -1) G[n2][n2] += g;
        if (n1 != -1 && n2 != -1) { G[n1][n2] -= g; G[n2][n1] -= g; }
    }
    void stamp_resistor(Matrix& G, const std::string& n1, const std::string& n2, double R) {
        stamp_resistor(G, m_node(n1), m_node(n2), R);
    }

    void stamp_capacitor(int n1, int n2, double C_val) {
        if (n1 != -1) m_C(n1, n1) += C_val;
        if (n2 != -1) m_C(n2, n2) += C_val;
        if (n1 != -1 && n2 != -1) { m_C(n1, n2) -= C_val; m_C(n2, n1) -= C_val; }
    }
    void stamp_capacitor(const std::string& n1, const std::string& n2, double C_val) {
        stamp_capacitor(m_node(n1), m_node(n2), C_val);
    }

    void stamp_voltage_source(const std::string& n_pos, const std::string& n_neg, double V) {
        int n1 = m_node(n_pos);
        if (n_neg != "gnd" && n_neg != "0") {
             throw std::runtime_error("This simple MNA only supports V-sources connected to ground.");
        }
        double G_large = 1e12;
        m_G_static[n1][n1] += G_large;
        m_b[n1] += V * G_large;
    }

    // --- Non-Linear Tube Model ---

    /**
     * @brief Adds the triode's contribution to the Jacobian and residual for Newton-Raphson.
     * @param J The Jacobian matrix (the system matrix 'A').
     * @param f The residual vector, which will have I_nonlinear subtracted from it.
     * @param x The current solution guess vector (node voltages).
     */
    void add_triode_nr_stamps(Matrix& J, Vector& f, const Vector& x) {
        double Vp = x[m_p_node], Vg = x[m_g_node], Vk = x[m_k_node];
        double Vpk = Vp - Vk, Vgk = Vg - Vk;

        // --- Plate Current (Koren Model) ---
        double Ip = 0.0, Gp = 0.0, Gm = 0.0;
        // ROBUSTNESS: check for E1 > 0 to prevent pow(negative, non-integer) -> NaN
        double E1 = Vpk / TUBE_MU + Vgk;
        if (E1 > 0) {
            double log_E1 = log(E1);
            Ip = exp(TUBE_EX * log_E1 + log(E1/TUBE_KG)); // More stable than pow

            double dIp_dE1 = Ip * TUBE_EX / E1;
            Gp = dIp_dE1 / TUBE_MU;
            Gm = dIp_dE1;
        }

        // --- Grid Current (Diode Model) ---
        double Ig = 0.0, Gg = 0.0;
        if (Vgk > 0) { // simplified threshold for stability
            Ig = TUBE_K_GRID * Vgk * Vgk * Vgk; // Smoother than a hard knee
            Gg = 3.0 * TUBE_K_GRID * Vgk * Vgk;
        }

        // --- Stamp contributions into Residual (f) and Jacobian (J) ---
        // f(x) = ... - I_nonlinear(x)
        // J(x) = ... + G_nonlinear(x)
        f[m_p_node] -= Ip;
        f[m_k_node] += Ip;
        f[m_g_node] -= Ig;
        f[m_k_node] += Ig;

        // Jacobian plate conductance stamps
        J[m_p_node][m_p_node] += Gp;
        J[m_p_node][m_k_node] -= Gp;
        J[m_p_node][m_g_node] += Gm;
        J[m_p_node][m_k_node] -= Gm;

        J[m_k_node][m_p_node] -= Gp;
        J[m_k_node][m_k_node] += Gp;
        // **FIXED HERE**: The sign was wrong in the original code. d(-Ip)/dVg = -Gm
        J[m_k_node][m_g_node] -= Gm;
        J[m_k_node][m_k_node] += Gm;

        // Jacobian grid conductance stamps
        J[m_g_node][m_g_node] += Gg;
        J[m_g_node][m_k_node] -= Gg;
        J[m_k_node][m_g_node] -= Gg;
        J[m_k_node][m_k_node] += Gg;
    }

    // --- Dense LU Solver and Helpers ---
    Vector matrix_vector_mult(const Matrix& A, const Vector& v) {
        int n = A.size();
        Vector result(n, 0.0);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                result[i] += A[i][j] * v[j];
            }
        }
        return result;
    }

    Vector lu_decompose(Matrix& A) {
        int n = A.size();
        Vector p(n);
        for (int i = 0; i < n; ++i) p[i] = i;

        for (int i = 0; i < n; ++i) {
            int max_j = i;
            for (int j = i + 1; j < n; ++j) {
                if (std::abs(A[j][i]) > std::abs(A[max_j][i])) max_j = j;
            }
            if (max_j != i) {
                std::swap(A[i], A[max_j]);
                std::swap(p[i], p[max_j]);
            }

            if (std::abs(A[i][i]) < 1e-15) { /* Matrix is singular */ }

            for (int j = i + 1; j < n; ++j) {
                A[j][i] /= A[i][i];
                for (int k = i + 1; k < n; ++k) {
                    A[j][k] -= A[j][i] * A[i][k];
                }
            }
        }
        return p;
    }

    Vector lu_solve(const Matrix& LU, const Vector& p_in, const Vector& b) {
        int n = LU.size();
        Vector x(n), y(n);

        // Permute b vector
        Vector pb(n);
        for(int i=0; i<n; ++i) pb[i] = b[static_cast<int>(p_in[i])];

        // Forward substitution: Ly = Pb
        for (int i = 0; i < n; ++i) {
            y[i] = pb[i];
            for (int j = 0; j < i; ++j) y[i] -= LU[i][j] * y[j];
        }

        // Backward substitution: Ux = y
        for (int i = n - 1; i >= 0; --i) {
            x[i] = y[i];
            for (int j = i + 1; j < n; ++j) x[i] -= LU[i][j] * x[j];
            x[i] /= LU[i][i];
        }
        return x;
    }
};
