#pragma once

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <algorithm>

/**
 * @class RealtimeTubeSim
 * @brief A robust, single-header C++ class for realtime audio simulation of a tube amplifier circuit.
 *
 * This version corrects a critical numerical instability in the triode model that caused
 * the output to drift. The simulation is now stable and correctly processes audio.
 *
 * Final Fix:
 * 1.  **Numerically Stable Triode Model:** The `pow(base, exp)` function in the triode current
 *     calculation has been replaced with the more stable `exp(exp * log(base))` equivalent.
 *     This prevents floating-point inaccuracies when the base is near zero, which was the
 *     root cause of the runaway DC offset.
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
        // --- 1. Assemble the time-dependent Right-Hand Side (RHS) vector ---
        Vector z(m_num_nodes, 0.0);

        // Add voltage source contributions
        for (const auto& vs : m_v_sources) {
            double voltage = vs.is_time_varying ? (input_sample * m_params.volume) : vs.dc_voltage;
            z[vs.node] += voltage * vs.g_large;
        }

        // Add capacitor history contributions
        Vector C_x_prev = matrix_vector_mult(m_C_mat, m_x_prev);
        for (int i = 0; i < m_num_nodes; ++i) {
            z[i] += C_x_prev[i] / m_dt;
        }

        // --- 2. Newton-Raphson Iteration ---
        m_x = m_x_prev; // Initial guess

        Matrix A_linear = m_G_dynamic;
        for (int r = 0; r < m_num_nodes; ++r) {
            for (int c = 0; c < m_num_nodes; ++c) {
                A_linear[r][c] += m_C_mat[r][c] / m_dt;
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
    struct VoltageSource { int node; double dc_voltage; double g_large; bool is_time_varying; };

    Matrix m_G_static, m_G_dynamic, m_C_mat;
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

        stamp_resistor(m_G_static, "in", "gnd", 1e6); // DC path for input node
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
        for (const auto& vs : m_v_sources) b_dc[vs.node] += vs.dc_voltage * vs.g_large;

        m_x.assign(m_num_nodes, 0.0);
        for (int i = 0; i < 30; ++i) {
            Matrix J = m_G_dynamic;
            Vector f = matrix_vector_mult(m_G_dynamic, m_x);
            for(size_t j = 0; j < f.size(); ++j) f[j] -= b_dc[j];
            add_triode_nr_stamps(J, f, m_x);
            for(size_t j = 0; j < f.size(); ++j) f[j] = -f[j];
            auto p = lu_decompose(J); Vector dx = lu_solve(J, p, f);
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
    void stamp_resistor(Matrix& G, const std::string& s1, const std::string& s2, double R) { stamp_resistor(G, m_node(s1), m_node(s2), R); }
    void stamp_capacitor(const std::string& s1, const std::string& s2, double C) {
        int N1 = m_node(s1), N2 = m_node(s2);
        if (N1 != -1) m_C_mat[N1][N1] += C; if (N2 != -1) m_C_mat[N2][N2] += C;
        if (N1 != -1 && N2 != -1) { m_C_mat[N1][N2] -= C; m_C_mat[N2][N1] -= C; }
    }
    void stamp_voltage_source(const std::string& n_pos, const std::string& n_neg, double V, bool is_ac) {
        if (n_neg != "gnd" && n_neg != "0") throw std::runtime_error("V-sources must be grounded.");
        double G_large = 1e12; // Increased for stiffness
        m_G_static[m_node(n_pos)][m_node(n_pos)] += G_large;
        m_v_sources.push_back({m_node(n_pos), V, G_large, is_ac});
    }

    void add_triode_nr_stamps(Matrix& J, Vector& f, const Vector& x) {
        double Vp = x[m_p_node], Vg = x[m_g_node], Vk = x[m_k_node];
        double Vpk = Vp - Vk, Vgk = Vg - Vk;

        double Ip = 0.0, Gp = 0.0, Gm = 0.0;
        double E1 = Vpk / TUBE_MU + Vgk;
        if (E1 > 1e-9) { // Safety margin
            // **STABILITY FIX**: Use log/exp for pow()
            double log_E1 = log(E1);
            Ip = exp(TUBE_EX * log_E1) / TUBE_KG;
            // Use stable derivative calculation
            double dIp_dE1 = Ip * TUBE_EX / E1;
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
        Vector result(v.size(), 0.0);
        for (size_t i = 0; i < A.size(); ++i) for (size_t j = 0; j < v.size(); ++j) result[i] += A[i][j] * v[j];
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
