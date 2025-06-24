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
 * @brief A flexible, component-based C++ realtime audio circuit simulator.
 *
 * This version has been refactored to support arbitrary circuit topologies with multiple
 * non-linear triode stages.
 *
 * Key Features:
 * - **Component API:** Add resistors, capacitors, voltage sources, and triodes to nodes programmatically.
 * - **Multi-Triode Support:** Simulate complex multi-stage amplifiers.
 * - **Advanced 12AX7 Model:** Implements a detailed SPICE model for high-fidelity simulation.
 * - **Numerical Jacobian:** Automatically and accurately calculates partial derivatives for the
 *   non-linear models using finite differences, ensuring stability and extensibility.
 *
 * Example Usage (to build the original preamp circuit):
 *
 *   RealtimeTubeSim sim(48000.0);
 *   // Define nodes
 *   sim.add_node("in"); sim.add_node("g"); sim.add_node("p"); sim.add_node("k");
 *   sim.add_node("V_P"); sim.add_node("ts1"); sim.add_node("ts2"); sim.add_node("ts3");
 *   sim.add_node("out"); sim.add_node("out_final");
 *
 *   // Add components
 *   sim.add_resistor("in", "gnd", 1e6);
 *   sim.add_resistor("g", "gnd", 1e6);
 *   // ... and so on for all R and C components ...
 *
 *   // Add the triode with specified model parameters
 *   sim.add_triode("p", "g", "k", 96.20, 1.437, 613.4, 740.3, 1672.0, 2000.0);
 *
 *   // Add sources
 *   sim.add_voltage_source("V_P", "gnd", 300.0, false);
 *   sim.add_voltage_source("in", "gnd", 0.0, true); // Input signal
 *
 *   sim.prepare_to_play(); // Finalizes matrices and solves DC point
 */
class RealtimeTubeSim {
public:
    using Matrix = std::vector<std::vector<double>>;
    using Vector = std::vector<double>;

    explicit RealtimeTubeSim(double sample_rate) {
        if (sample_rate <= 0) throw std::invalid_argument("Sample rate must be positive.");
        m_sample_rate = sample_rate;
        m_dt = 1.0 / sample_rate;
    }

    // --- 1. Circuit Building API ---

    int add_node(const std::string& name) {
        if (m_node_map.find(name) == m_node_map.end()) {
            m_node_map[name] = m_num_nodes++;
        }
        return m_node_map[name];
    }

    void add_resistor(const std::string& n1_str, const std::string& n2_str, double R) {
        m_resistors.push_back({m_node(n1_str), m_node(n2_str), R});
    }

    void add_capacitor(const std::string& n1_str, const std::string& n2_str, double C) {
        m_capacitors.push_back({m_node(n1_str), m_node(n2_str), C});
    }

    void add_voltage_source(const std::string& n_pos_str, const std::string& n_neg_str, double V, bool is_time_varying) {
        if (n_neg_str != "gnd" && n_neg_str != "0") throw std::runtime_error("V-sources must be grounded.");
        m_v_sources.push_back({m_node(n_pos_str), V, is_time_varying});
    }

    void add_triode(const std::string& p_str, const std::string& g_str, const std::string& k_str,
                    double mu, double ex, double kg1, double kp, double kvb, double rgi) {
        Triode t;
        t.p_node = m_node(p_str); t.g_node = m_node(g_str); t.k_node = m_node(k_str);
        t.mu = mu; t.ex = ex; t.kg1 = kg1; t.kp = kp; t.kvb = kvb; t.rgi = rgi;
        m_triodes.push_back(t);
    }

    // --- 2. Simulation Control ---

    void prepare_to_play() {
        if (m_is_prepared) return;
        m_G_static.assign(m_num_nodes, Vector(m_num_nodes, 0.0));
        m_C_mat.assign(m_num_nodes, Vector(m_num_nodes, 0.0));
        m_x.assign(m_num_nodes, 0.0);
        m_x_prev.assign(m_num_nodes, 0.0);

        for(const auto& r : m_resistors) stamp_resistor(m_G_static, r.n1, r.n2, r.R);
        for(const auto& c : m_capacitors) stamp_capacitor(c.n1, c.n2, c.C);
        for(const auto& vs : m_v_sources) {
            m_G_static[vs.node][vs.node] += G_LARGE;
        }
        // Stamp static part of triode models (e.g. RCP)
        for (const auto& t : m_triodes) {
            stamp_resistor(m_G_static, t.p_node, t.k_node, 1e9); // RCP = 1G
        }

        m_G_dynamic = m_G_static;
        solve_dc();
        m_is_prepared = true;
    }

    double process_sample(double input_sample) {
        if (!m_is_prepared) prepare_to_play();

        Vector b(m_num_nodes, 0.0);
        for (const auto& vs : m_v_sources) {
            double voltage = vs.is_time_varying ? input_sample : vs.dc_voltage;
            b[vs.node] += voltage * G_LARGE;
        }

        m_x = m_x_prev;

        for (int i=0; i<m_nr_max_iter; ++i) {
            Matrix J = m_G_dynamic;
            for(int r=0; r<m_num_nodes; ++r) for(int c=0; c<m_num_nodes; ++c) {
                J[r][c] += m_C_mat[r][c] / m_dt;
            }

            Vector Gx = matrix_vector_mult(m_G_dynamic, m_x);
            Vector x_diff = m_x;
            for(size_t k=0; k<m_x.size(); ++k) x_diff[k] -= m_x_prev[k];
            Vector C_x_diff_dt = matrix_vector_mult(m_C_mat, x_diff);
            for(size_t k=0; k<m_x.size(); ++k) C_x_diff_dt[k] /= m_dt;

            Vector f(m_num_nodes);
            for(size_t k=0; k<f.size(); ++k) f[k] = Gx[k] + C_x_diff_dt[k] - b[k];

            add_all_nonlinear_stamps(J, f, m_x);

            for(size_t j=0; j<f.size(); ++j) f[j] = -f[j];
            auto p = lu_decompose(J);
            Vector dx = lu_solve(J, p, f);

            double damping_factor = 1.0;
            double max_step = 0.0;
            for(double val : dx) max_step = std::max(max_step, std::abs(val));
            if (max_step > m_nr_damping) damping_factor = m_nr_damping / max_step;

            double norm = 0.0;
            for (int j=0; j<m_num_nodes; ++j) {
                double step = damping_factor * dx[j];
                m_x[j] += step; norm += step * step;
            }
            if (sqrt(norm) < m_nr_tolerance) break;
        }

        m_x_prev = m_x;
        // In a real plugin, you would add an output node name parameter.
        // For now, we assume the last added node is the output.
        return m_x.back();
    }

private:
    // --- Data Structures ---
    struct Resistor { int n1, n2; double R; };
    struct Capacitor { int n1, n2; double C; };
    struct VoltageSource { int node; double dc_voltage; bool is_time_varying; };
    struct Triode {
        int p_node, g_node, k_node;
        double mu, ex, kg1, kp, kvb, rgi;
    };

    std::vector<Resistor> m_resistors;
    std::vector<Capacitor> m_capacitors;
    std::vector<VoltageSource> m_v_sources;
    std::vector<Triode> m_triodes;

    Matrix m_G_static, m_G_dynamic, m_C_mat;
    Vector m_x, m_x_prev;
    int m_num_nodes = 0;
    std::map<std::string, int> m_node_map;
    bool m_is_prepared = false;

    // --- Simulation Parameters ---
    double m_sample_rate, m_dt;
    const int m_nr_max_iter = 15;
    const double m_nr_tolerance = 1e-6;
    const double m_nr_damping = 1.0;
    const double G_LARGE = 1e12;

    // --- Internal Methods ---
    int m_node(const std::string& name) {
        if (name == "gnd" || name == "0") return -1;
        if (m_node_map.find(name) == m_node_map.end()) {
            return add_node(name);
        }
        return m_node_map.at(name);
    }

    void solve_dc() {
        Vector b_dc(m_num_nodes, 0.0);
        for (const auto& vs : m_v_sources) b_dc[vs.node] += vs.dc_voltage * G_LARGE;
        m_x.assign(m_num_nodes, 0.0);
        for (int i=0; i<30; ++i) {
            Matrix J = m_G_dynamic;
            Vector f = matrix_vector_mult(m_G_dynamic, m_x);
            for(size_t j = 0; j < f.size(); ++j) f[j] -= b_dc[j];
            add_all_nonlinear_stamps(J, f, m_x);
            for(size_t j = 0; j < f.size(); ++j) f[j] = -f[j];
            auto p = lu_decompose(J); Vector dx = lu_solve(J, p, f);
            double norm = 0.0;
            for (int j=0; j<m_num_nodes; ++j) { m_x[j] += dx[j]; norm += dx[j] * dx[j]; }
            if (sqrt(norm) < m_nr_tolerance) break;
        }
        m_x_prev = m_x;
    }

    // --- Non-linear Model Implementation ---

    void add_all_nonlinear_stamps(Matrix& J, Vector& f, const Vector& x) {
        for (const auto& tube : m_triodes) {
            add_single_triode_stamps_numerical(tube, J, f, x);
        }
    }

    // New 12AX7 model based on SPICE subcircuit
    void get_triode_currents(const Triode& t, double Vp, double Vg, double Vk, double& Ip, double& Ig) {
        const double Vpk = Vp - Vk;
        const double Vgk = Vg - Vk;

        // Plate Current (Duncan Munro model)
        const double kvb_vpk_sq = t.kvb + Vpk * Vpk;
        if (kvb_vpk_sq <= 0) {
            Ip = 0;
        } else {
            const double inner_exp = t.kp * (1.0 / t.mu + Vgk / sqrt(kvb_vpk_sq));
            double E1 = 0.0;
            // Use log(1+exp(x)) "softplus" function for stability
            if (inner_exp < -100) E1 = 0.0; // exp(x) -> 0
            else if (inner_exp > 100) E1 = (Vpk / t.kp) * inner_exp; // exp(x) -> inf, log(exp(x)) -> x
            else E1 = (Vpk / t.kp) * log(1.0 + exp(inner_exp));

            // The 0.5*(pwr+pwrs) term is a smooth way of saying: if (E1 > 0) current = E1^ex else 0
            if (E1 > 0) {
                Ip = pow(E1, t.ex) / t.kg1;
            } else {
                Ip = 0;
            }
        }

        // Grid Current (approximated diode with series resistor RGI)
        if (Vgk > 0) {
            Ig = (Vgk * Vgk) / (t.rgi + Vgk); // Simple non-linear diode-like current
        } else {
            Ig = 0;
        }
    }

    void add_single_triode_stamps_numerical(const Triode& t, Matrix& J, Vector& f, const Vector& x) {
        double Vp = x[t.p_node], Vg = x[t.g_node], Vk = x[t.k_node];

        // 1. Calculate base currents at x
        double Ip0, Ig0;
        get_triode_currents(t, Vp, Vg, Vk, Ip0, Ig0);

        // 2. Add base currents to residual vector f
        f[t.p_node] += Ip0;
        f[t.k_node] -= Ip0;
        f[t.g_node] += Ig0;
        f[t.k_node] -= Ig0;

        // 3. Calculate partial derivatives using finite differences and add to Jacobian J
        const double delta = 1e-6; // Small voltage perturbation

        // d/dVp
        double Ip_p, Ig_p;
        get_triode_currents(t, Vp + delta, Vg, Vk, Ip_p, Ig_p);
        double dIp_dVp = (Ip_p - Ip0) / delta;
        double dIg_dVp = (Ig_p - Ig0) / delta;

        // d/dVg
        double Ip_g, Ig_g;
        get_triode_currents(t, Vp, Vg + delta, Vk, Ip_g, Ig_g);
        double dIp_dVg = (Ip_g - Ip0) / delta;
        double dIg_dVg = (Ig_g - Ig0) / delta;

        // d/dVk
        double Ip_k, Ig_k;
        get_triode_currents(t, Vp, Vg, Vk + delta, Ip_k, Ig_k);
        double dIp_dVk = (Ip_k - Ip0) / delta;
        double dIg_dVk = (Ig_k - Ig0) / delta;

        // Stamp Jacobian
        J[t.p_node][t.p_node] += dIp_dVp; J[t.p_node][t.g_node] += dIp_dVg; J[t.p_node][t.k_node] += dIp_dVk;
        J[t.g_node][t.p_node] += dIg_dVp; J[t.g_node][t.g_node] += dIg_dVg; J[t.g_node][t.k_node] += dIg_dVk;
        J[t.k_node][t.p_node] -= (dIp_dVp + dIg_dVp);
        J[t.k_node][t.g_node] -= (dIp_dVg + dIg_dVg);
        J[t.k_node][t.k_node] -= (dIp_dVk + dIg_dVk);
    }

    // --- MNA Stamping Helpers ---
    void stamp_resistor(Matrix& G, int n1, int n2, double R) {
        if (R <= 1e-9) R = 1e-9; double g = 1.0 / R;
        if (n1 != -1) G[n1][n1] += g; if (n2 != -1) G[n2][n2] += g;
        if (n1 != -1 && n2 != -1) { G[n1][n2] -= g; G[n2][n1] -= g; }
    }
    void stamp_capacitor(int n1, int n2, double C) {
        if (n1 != -1) m_C_mat[n1][n1] += C; if (n2 != -1) m_C_mat[n2][n2] += C;
        if (n1 != -1 && n2 != -1) { m_C_mat[n1][n2] -= C; m_C_mat[n2][n1] -= C; }
    }

    // --- Linear Algebra ---
    static RealtimeTubeSim::Vector matrix_vector_mult(const Matrix& A, const Vector& v) {
        Vector result(v.size(), 0.0);
        for (size_t i = 0; i < A.size(); ++i) for (size_t j = 0; j < v.size(); ++j) result[i] += A[i][j] * v[j];
        return result;
    }

    static RealtimeTubeSim::Vector lu_decompose(Matrix& A) {
        int n = A.size(); Vector p(n); for (int i = 0; i < n; ++i) p[i] = i;
        for (int i = 0; i < n; ++i) {
            int max_j = i;
            for (int j = i + 1; j < n; ++j) if (std::abs(A[j][i]) > std::abs(A[max_j][i])) max_j = j;
            if (max_j != i) { std::swap(A[i], A[max_j]); std::swap(p[i], p[max_j]); }
            if (std::abs(A[i][i]) < 1e-20) continue;
            for (int j = i + 1; j < n; ++j) {
                A[j][i] /= A[i][i];
                for (int k = i + 1; k < n; ++k) A[j][k] -= A[j][i] * A[i][k];
            }
        } return p;
    }

    static RealtimeTubeSim::Vector lu_solve(const Matrix& LU, const Vector& p_in, const Vector& b) {
        int n = LU.size(); Vector x(n), y(n);
        for(int i=0; i<n; ++i) {
            y[i] = b[static_cast<int>(p_in[i])];
            for (int j = 0; j < i; ++j) y[i] -= LU[i][j] * y[j];
        }
        for (int i = n - 1; i >= 0; --i) {
            x[i] = y[i];
            for (int j = i + 1; j < n; ++j) x[i] -= LU[i][j] * x[j];
            if (std::abs(LU[i][i]) > 1e-20) x[i] /= LU[i][i]; else x[i] = 0;
        } return x;
    }
};
