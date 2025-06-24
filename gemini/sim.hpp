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
 * @brief A robust, single-header C++ class for realtime audio simulation of a tube amplifier circuit.
 *
 * This version corrects a critical bug in the time-varying source implementation that caused
 * the output to ramp up instead of processing audio.
 *
 * Key Improvements:
 * 1.  **Correct Time-Varying Source Model:** The input voltage source is now correctly
 *     implemented within the MNA framework. The right-hand-side vector is properly
 *     constructed at each time step to reflect V(t) * G_large, allowing the signal
 *     to pass into the simulation.
 * 2.  **Explicit Voltage Source Management:** A new struct manages voltage sources,
 *     making the code cleaner and less error-prone.
 * 3.  **Robust DC Solver:** The DC solver is updated to use the new source management,
 *     ensuring a correct quiescent point calculation.
 *
 * --- Simulated Circuit (Conceptual SPICE Netlist) ---
 * (Circuit remains the same as previous description)
 */
class RealtimeTubeSim {
public:
    using Matrix = std::vector<std::vector<double>>;
    using Vector = std::vector<double>;

    explicit RealtimeTubeSim(double sample_rate) {
        if (sample_rate <= 0) {
            throw std::invalid_argument("Sample rate must be positive.");
        }
        m_sample_rate = sample_rate;
        m_dt = 1.0 / sample_rate;

        build_system();
        m_x.assign(m_num_nodes, 0.0);
        m_x_prev.assign(m_num_nodes, 0.0);

        set_volume(0.5); set_bass(0.5); set_mid(0.5); set_treble(0.5); set_master(1.0);

        solve_dc();
    }

    void set_volume(double value) { m_params.volume = std::max(0.0, std::min(1.0, value)); }
    void set_bass(double value) { m_params.bass = std::max(0.0, std::min(1.0, value)); update_pots(); }
    void set_mid(double value) { m_params.mid = std::max(0.0, std::min(1.0, value)); update_pots(); }
    void set_treble(double value) { m_params.treble = std::max(0.0, std::min(1.0, value)); update_pots(); }
    void set_master(double value) { m_params.master = std::max(0.0, std::min(1.0, value)); update_pots(); }

    double process_sample(double input_sample) {
        // --- Assemble the time-dependent Right-Hand Side (RHS) vector ---
        // z = b_n + (C/dt) * x_{n-1}
        Vector z(m_num_nodes, 0.0);

        // Add voltage source contributions to z
        for (const auto& vs : m_v_sources) {
            double voltage = vs.dc_voltage;
            if (vs.is_time_varying) {
                voltage = input_sample * m_params.volume;
            }
            z[vs.node] += voltage * vs.g_large;
        }

        // Add capacitor history contributions to z
        for (int r = 0; r < m_num_nodes; ++r) {
            for (int c = 0; c < m_num_nodes; ++c) {
                if (m_C(r, c) != 0.0) {
                    z[r] += (m_C(r, c) / m_dt) * m_x_prev[c];
                }
            }
        }

        // --- Newton-Raphson Iteration ---
        m_x = m_x_prev; // Initial guess

        // Pre-calculate the linear part of the A matrix (G + C/dt)
        Matrix A_linear = m_G_dynamic;
        for (int r = 0; r < m_num_nodes; ++r) {
            for (int c = 0; c < m_num_nodes; ++c) {
                if (m_C(r, c) != 0.0) {
                    A_linear[r][c] += m_C(r, c) / m_dt;
                }
            }
        }

        for (int i = 0; i < m_nr_max_iter; ++i) {
            Matrix J = A_linear;
            Vector f = matrix_vector_mult(A_linear, m_x);
            for(size_t j = 0; j < f.size(); ++j) f[j] -= z[j];

            add_triode_nr_stamps(J, f, m_x);

            for(size_t j = 0; j < f.size(); ++j) f[j] = -f[j];

            auto p = lu_decompose(J);
            Vector dx = lu_solve(J, p, f);

            double damping_factor = 1.0;
            double max_step = 0.0;
            for(double val : dx) max_step = std::max(max_step, std::abs(val));
            if (max_step > m_nr_damping) {
                damping_factor = m_nr_damping / max_step;
            }

            double norm = 0.0;
            for (int j = 0; j < m_num_nodes; ++j) {
                double step = damping_factor * dx[j];
                m_x[j] += step;
                norm += step * step;
            }

            if (sqrt(norm) < m_nr_tolerance) break;
        }

        m_x_prev = m_x;
        return m_x[m_node_map["out_final"]];
    }

private:
    struct MatrixView { Matrix& mat; double& operator()(size_t r, size_t c) { return mat.at(r).at(c); } };
    struct VoltageSource { int node; double dc_voltage; double g_large; bool is_time_varying; };

    Matrix m_G_static, m_G_dynamic, m_C_mat;
    MatrixView m_C { m_C_mat };
    Vector m_x, m_x_prev;
    std::vector<VoltageSource> m_v_sources;

    int m_num_nodes = 0;
    std::map<std::string, int> m_node_map;

    double m_sample_rate, m_dt;
    const int m_nr_max_iter = 15;
    const double m_nr_tolerance = 1e-6;
    const double m_nr_damping = 1.0;

    struct Potentiometer { int n1, n2, n_wiper; double total_resistance; };
    Potentiometer m_pot_treble, m_pot_bass, m_pot_mid, m_pot_master;

    struct ControlParams { double volume, bass, mid, treble, master; } m_params;

    const double TUBE_MU = 100.0, TUBE_EX = 1.4, TUBE_KG = 1060.0;
    const double TUBE_KP = 600.0, TUBE_K_GRID = 2000.0;
    int m_p_node, m_g_node, m_k_node;

    void build_system() {
        auto add_node = [&](const std::string& name) { m_node_map[name] = m_num_nodes++; };
        add_node("in"); add_node("g"); add_node("p"); add_node("k");
        add_node("V_P"); add_node("ts1"); add_node("ts2"); add_node("ts3");
        add_node("out"); add_node("out_final");

        m_G_static.assign(m_num_nodes, Vector(m_num_nodes, 0.0));
        m_C_mat.assign(m_num_nodes, Vector(m_num_nodes, 0.0));

        stamp_resistor(m_G_static, "g", "gnd", 1e6);
        stamp_resistor(m_G_static, "V_P", "p", 100e3);
        stamp_resistor(m_G_static, "k", "gnd", 1.5e3);
        stamp_resistor(m_G_static, "out_final", "gnd", 1e6);

        stamp_capacitor("in", "g", 22e-9); stamp_capacitor("k", "gnd", 22e-6);
        stamp_capacitor("p", "ts1", 22e-9); stamp_capacitor("ts1", "ts3", 250e-12);
        stamp_capacitor("ts2", "ts3", 22e-9); stamp_capacitor("out", "out_final", 100e-9);

        stamp_voltage_source("V_P", "gnd", 300.0, false);
        stamp_voltage_source("in", "gnd", 0.0, true);

        m_pot_treble = {m_node("ts1"), m_node("ts2"), -1, 250e3};
        m_pot_bass = {m_node("ts2"), m_node("gnd"), -1, 1e6};
        m_pot_mid = {m_node("ts3"), m_node("gnd"), -1, 25e3};
        m_pot_master = {m_node("ts3"), m_node("gnd"), m_node("out"), 1e6};

        m_p_node = m_node("p"); m_g_node = m_node("g"); m_k_node = m_node("k");
    }

    void update_pots() {
        m_G_dynamic = m_G_static;
        auto taper = [](double val) { return val * val; };
        stamp_resistor(m_G_dynamic, m_pot_treble.n1, m_pot_treble.n2, std::max(1.0, m_pot_treble.total_resistance * taper(m_params.treble)));
        stamp_resistor(m_G_dynamic, m_pot_bass.n1, m_pot_bass.n2, std::max(1.0, m_pot_bass.total_resistance * taper(1.0 - m_params.bass)));
        stamp_resistor(m_G_dynamic, m_pot_mid.n1, m_pot_mid.n2, std::max(1.0, m_pot_mid.total_resistance * m_params.mid));
        double r1_val = m_pot_master.total_resistance * taper(m_params.master);
        double r2_val = m_pot_master.total_resistance - r1_val;
        stamp_resistor(m_G_dynamic, m_pot_master.n1, m_pot_master.n_wiper, std::max(1.0, r1_val));
        stamp_resistor(m_G_dynamic, m_pot_master.n_wiper, m_pot_master.n2, std::max(1.0, r2_val));
    }

    void solve_dc() {
        Vector b_dc(m_num_nodes, 0.0);
        for (const auto& vs : m_v_sources) {
            b_dc[vs.node] += vs.dc_voltage * vs.g_large;
        }

        m_x.assign(m_num_nodes, 0.0);
        for (int i = 0; i < 30; ++i) {
            Matrix J = m_G_dynamic;
            Vector f = matrix_vector_mult(m_G_dynamic, m_x);
            for(size_t j = 0; j < f.size(); ++j) f[j] -= b_dc[j];
            add_triode_nr_stamps(J, f, m_x);
            for(size_t j = 0; j < f.size(); ++j) f[j] = -f[j];

            auto p = lu_decompose(J);
            Vector dx = lu_solve(J, p, f);

            double norm = 0.0;
            for (int j = 0; j < m_num_nodes; ++j) { m_x[j] += dx[j]; norm += dx[j] * dx[j]; }
            if (sqrt(norm) < m_nr_tolerance) break;
        }
        m_x_prev = m_x;
    }

    int m_node(const std::string& name) { return (name == "gnd" || name == "0") ? -1 : m_node_map.at(name); }
    void stamp_resistor(Matrix& G, int n1, int n2, double R) {
        if (R <= 1e-9) R = 1e-9; double g = 1.0 / R;
        if (n1 != -1) G[n1][n1] += g; if (n2 != -1) G[n2][n2] += g;
        if (n1 != -1 && n2 != -1) { G[n1][n2] -= g; G[n2][n1] -= g; }
    }
    void stamp_resistor(Matrix& G, const std::string& n1, const std::string& n2, double R) { stamp_resistor(G, m_node(n1), m_node(n2), R); }
    void stamp_capacitor(const std::string& n1, const std::string& n2, double C) {
        int N1 = m_node(n1), N2 = m_node(n2);
        if (N1 != -1) m_C(N1, N1) += C; if (N2 != -1) m_C(N2, N2) += C;
        if (N1 != -1 && N2 != -1) { m_C(N1, N2) -= C; m_C(N2, N1) -= C; }
    }
    void stamp_voltage_source(const std::string& n_pos, const std::string& n_neg, double V, bool is_ac) {
        int n1 = m_node(n_pos);
        if (n_neg != "gnd" && n_neg != "0") throw std::runtime_error("V-sources must be grounded.");
        double G_large = 1e9;
        m_G_static[n1][n1] += G_large;
        m_v_sources.push_back({n1, V, G_large, is_ac});
    }

    void add_triode_nr_stamps(Matrix& J, Vector& f, const Vector& x) {
        double Vp = x[m_p_node], Vg = x[m_g_node], Vk = x[m_k_node];
        double Vpk = Vp - Vk, Vgk = Vg - Vk;

        double Ip = 0.0, Gp = 0.0, Gm = 0.0;
        double E1 = Vpk / TUBE_MU + Vgk;
        if (E1 > 1e-6) {
            double E1_ex = pow(E1, TUBE_EX);
            Ip = E1_ex / TUBE_KG;
            double dIp_dE1 = TUBE_EX * pow(E1, TUBE_EX - 1.0) / TUBE_KG;
            Gm = dIp_dE1;
            Gp = Gm / TUBE_MU;
        }

        double Ig = 0.0, Gg = 0.0;
        if (Vgk > 0) {
            Ig = TUBE_K_GRID * Vgk * Vgk * Vgk; Gg = 3.0 * TUBE_K_GRID * Vgk * Vgk;
        }

        f[m_p_node] -= Ip; f[m_k_node] += Ip;
        f[m_g_node] -= Ig; f[m_k_node] += Ig;

        J[m_p_node][m_p_node] += Gp; J[m_p_node][m_g_node] += Gm; J[m_p_node][m_k_node] -= (Gp + Gm);
        J[m_g_node][m_g_node] += Gg; J[m_g_node][m_k_node] -= Gg;
        J[m_k_node][m_p_node] -= Gp; J[m_k_node][m_g_node] -= (Gm + Gg); J[m_k_node][m_k_node] += (Gp + Gm + Gg);
    }

    Vector matrix_vector_mult(const Matrix& A, const Vector& v) {
        int n = A.size(); Vector result(n, 0.0);
        for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) result[i] += A[i][j] * v[j];
        return result;
    }
    Vector lu_decompose(Matrix& A) {
        int n = A.size(); Vector p(n); for (int i = 0; i < n; ++i) p[i] = i;
        for (int i = 0; i < n; ++i) {
            int max_j = i;
            for (int j = i + 1; j < n; ++j) if (std::abs(A[j][i]) > std::abs(A[max_j][i])) max_j = j;
            if (max_j != i) { std::swap(A[i], A[max_j]); std::swap(p[i], p[max_j]); }
            if (std::abs(A[i][i]) < 1e-15) continue;
            for (int j = i + 1; j < n; ++j) {
                A[j][i] /= A[i][i];
                for (int k = i + 1; k < n; ++k) A[j][k] -= A[j][i] * A[i][k];
            }
        } return p;
    }
    Vector lu_solve(const Matrix& LU, const Vector& p_in, const Vector& b) {
        int n = LU.size(); Vector x(n), y(n);
        for(int i=0; i<n; ++i) {
            y[i] = b[static_cast<int>(p_in[i])];
            for (int j = 0; j < i; ++j) y[i] -= LU[i][j] * y[j];
        }
        for (int i = n - 1; i >= 0; --i) {
            x[i] = y[i];
            for (int j = i + 1; j < n; ++j) x[i] -= LU[i][j] * x[j];
            if (std::abs(LU[i][i]) > 1e-15) x[i] /= LU[i][i]; else x[i] = 0;
        } return x;
    }
};
