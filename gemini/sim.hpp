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
 * @brief A flexible, component-based C++ realtime audio circuit simulator with dynamic parameters.
 *
 * This version adds support for interactive, realtime control of potentiometers and variable resistors.
 *
 * Key Features:
 * - **Dynamic Parameters:** Add potentiometers and variable resistors by name, and control them
 *   in realtime with the `set_parameter()` method.
 * - **Optimized Performance:** A "dirty flag" system ensures that matrix recalculations only
 *   occur when a parameter changes, not on every sample.
 * - **Taper Support:** Models both Linear ('L') and Audio/Logarithmic ('A') tapers.
 * - **Component API:** Build arbitrary circuits by adding nodes, R, C, V, and Triode components.
 * - **Advanced 12AX7 Model:** Uses a detailed SPICE model with a numerically-computed Jacobian.
 *
 * Example Usage (to build the original preamp with tone stack):
 *
 *   RealtimeTubeSim sim(48000.0);
 *   // Define all nodes first
 *   sim.add_node("in"); sim.add_node("g"); sim.add_node("p"); sim.add_node("k");
 *   sim.add_node("V_P"); sim.add_node("ts1"); sim.add_node("ts2"); sim.add_node("ts3");
 *   sim.add_node("out"); sim.add_node("out_final");
 *
 *   // --- Static Components ---
 *   sim.add_resistor("g", "gnd", 1e6);        // Grid leak
 *   sim.add_resistor("V_P", "p", 100e3);      // Plate load
 *   sim.add_resistor("k", "gnd", 1.5e3);       // Cathode resistor
 *   sim.add_resistor("out_final", "gnd", 1e6);// Final load
 *   sim.add_capacitor("in", "g", 22e-9);      // Input cap
 *   sim.add_capacitor("k", "gnd", 22e-6);     // Cathode bypass
 *   sim.add_capacitor("p", "ts1", 22e-9);     // Coupling cap to tone stack
 *   sim.add_capacitor("ts1", "ts3", 250e-12);  // Treble cap
 *   sim.add_capacitor("ts2", "ts3", 22e-9);    // Mid cap
 *   sim.add_capacitor("out", "out_final", 100e-9); // Output cap
 *
 *   // --- Dynamic Components (Pots/Variable Resistors) ---
 *   sim.add_variable_resistor("treble", "ts1", "ts2", 250e3, 'L');
 *   sim.add_variable_resistor("bass", "ts2", "gnd", 1e6, 'A'); // Bass pot wired as rheostat
 *   sim.add_variable_resistor("mid", "ts3", "gnd", 25e3, 'L');
 *   sim.add_potentiometer("master", "ts3", "gnd", "out", 1e6, 'A');
 *
 *   // --- Active Components ---
 *   sim.add_triode("p", "g", "k", 96.2, 1.437, 613.4, 740.3, 1672.0, 2000.0);
 *
 *   // --- Sources ---
 *   sim.add_voltage_source("V_P", "gnd", 300.0, false);
 *   sim.add_voltage_source("in", "gnd", 0.0, true);
 *
 *   sim.prepare_to_play();
 *
 *   // In main loop:
 *   // sim.set_parameter("bass", 0.75);
 *   // output = sim.process_sample(input);
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

    [[nodiscard]] std::string name_of(int node) const {
        if (node == -1) return "0";
        return m_node_map_rev.find(node)->second;
    }

    // --- 1. Circuit Building API ---
    int add_node(const std::string& name) {
        if (m_node_map.find(name) == m_node_map.end()) {
            m_node_map[name] = m_num_nodes++;
            m_node_map_rev[m_num_nodes - 1] = name;
        }
        return m_node_map[name];
    }
    void add_resistor(const std::string& name, const std::string& n1, const std::string& n2, double R) { m_resistors.push_back({name, m_node(n1), m_node(n2), R}); }
    void add_capacitor(const std::string& name, const std::string& n1, const std::string& n2, double C) { m_capacitors.push_back({name, m_node(n1), m_node(n2), C}); }
    void add_voltage_source(const std::string& name, const std::string& n_pos, const std::string& n_neg, double V, bool is_ac) {
        if (n_neg != "gnd" && n_neg != "0") throw std::runtime_error("V-sources must be grounded.");
        m_v_sources.push_back({name, m_node(n_pos), V, is_ac});
    }
    void add_triode(const std::string& name, const std::string& p, const std::string& g, const std::string& k,
        double mu = 96.2,
        double ex = 1.437,
        double kg1 = 613.4,
        double kp = 740.3,
        double kvb = 1672.0,
        double rgi = 2000.0
    ) {
        m_triodes.push_back({name, m_node(p), m_node(g), m_node(k), mu, ex, kg1, kp, kvb, rgi});
    }
    void add_potentiometer(const std::string& name, const std::string& n1, const std::string& n2, const std::string& wiper, double R, char taper) {
        m_potentiometers[name] = {name, m_node(n1), m_node(n2), m_node(wiper), R, 0.5, taper};
    }
    void add_variable_resistor(const std::string& name, const std::string& n1, const std::string& n2, double R_max, char taper) {
        m_variable_resistors[name] = {name, m_node(n1), m_node(n2), R_max, 0.5, taper};
    }

    // --- 2. Simulation Control ---
    void set_parameter(const std::string& name, double value) {
        value = std::max(0.0, std::min(1.0, value));
        if (m_potentiometers.count(name)) {
            m_potentiometers.at(name).value = value;
            m_params_dirty = true;
        } else if (m_variable_resistors.count(name)) {
            m_variable_resistors.at(name).value = value;
            m_params_dirty = true;
        } else {
            // It's often better to fail silently in realtime audio
            // throw std::runtime_error("Parameter not found: " + name);
        }
    }

    void set_output_node(const std::string& name) {
        m_output_node = m_node(name);
    }

    void prepare_to_play() {
        if (m_is_prepared) return;
        m_G_static.assign(m_num_nodes, Vector(m_num_nodes, 0.0));
        m_C_mat.assign(m_num_nodes, Vector(m_num_nodes, 0.0));
        m_x.assign(m_num_nodes, 0.0);
        m_x_prev.assign(m_num_nodes, 0.0);
        m_x_prev2.assign(m_num_nodes, 0.0); // Init history

        for(const auto& r : m_resistors) stamp_resistor(m_G_static, r.n1, r.n2, r.R);
        for(const auto& c : m_capacitors) stamp_capacitor(c.n1, c.n2, c.C);
        for(const auto& vs : m_v_sources) m_G_static[vs.node][vs.node] += G_LARGE;
        for(const auto& t : m_triodes) stamp_resistor(m_G_static, t.p_node, t.k_node, 1e9);

        // select output_node if not already:
        if (m_output_node < 0) {
            m_output_node = m_x.size() - 1;
        }

        m_params_dirty = true; // Force initial calculation
        update_dynamic_components(); // Initial stamp of pots

        solve_dc();
        m_is_prepared = true;
    }

    double process_sample(double input_sample) {
        if (!m_is_prepared) prepare_to_play();
        if (m_params_dirty) update_dynamic_components();

        Vector b(m_num_nodes, 0.0);
        for (const auto& vs : m_v_sources) {
            double voltage = vs.is_time_varying ? input_sample : vs.dc_voltage;
            b[vs.node] += voltage * G_LARGE;
        }
#if 0
        // --- CONVERGENCE FIX 1: Extrapolate a better initial guess ---
        for (size_t i = 0; i < m_num_nodes; ++i) {
            m_x[i] = 2.0 * m_x_prev[i] - m_x_prev2[i];
        }

        for (int i = 0; i < m_nr_max_iter; ++i) {
            Matrix J = m_G_dynamic;
            for(int r=0; r<m_num_nodes; ++r) for(int c=0; c<m_num_nodes; ++c) J[r][c] += m_C_mat[r][c] / m_dt;

            Vector f = eval_f(m_x, b);
            add_all_nonlinear_stamps(J, f, m_x);

            // --- CONVERGENCE FIX 3: Adaptive Damping (Line Search) ---
            double initial_error_norm = vector_norm(f);

            for(size_t j=0; j<f.size(); ++j) f[j] = -f[j];
            auto p = lu_decompose(J);
            Vector dx = lu_solve(J, p, f);

            double step_damping = 1.0;
            Vector next_x = m_x;
            for(int k=0; k<10; ++k) { // Try up to 10 step reductions
                next_x = m_x;
                for(int n=0; n<m_num_nodes; ++n) next_x[n] += step_damping * dx[n];

                Vector next_f = eval_f(next_x, b);
                if (vector_norm(next_f) < initial_error_norm) {
                    m_x = next_x; // Good step, accept it
                    break;
                }
                step_damping /= 2.0; // Bad step, reduce size and try again
                if (k == 9) m_x = next_x; // Failsafe, accept last attempt
            }

            double norm = 0.0; for(double val : dx) norm += val * val;
            if (sqrt(norm * step_damping) < m_nr_tolerance) break;
        }

        m_x_prev2 = m_x_prev;
        m_x_prev = m_x;
#else
        m_x = m_x_prev;

        int iter;
        double norm_min = 1e6;
        for (iter=0; iter<m_nr_max_iter; ++iter) {
            Matrix J = m_G_dynamic;
            for(int r=0; r<m_num_nodes; ++r) for(int c=0; c<m_num_nodes; ++c) J[r][c] += m_C_mat[r][c] / m_dt;
            Vector Gx = matrix_vector_mult(m_G_dynamic, m_x);
            Vector x_diff = m_x; for(size_t k=0; k<m_x.size(); ++k) x_diff[k] -= m_x_prev[k];
            Vector C_x_diff_dt = matrix_vector_mult(m_C_mat, x_diff); for(size_t k=0; k<m_x.size(); ++k) C_x_diff_dt[k] /= m_dt;
            Vector f(m_num_nodes); for(size_t k=0; k<f.size(); ++k) f[k] = Gx[k] + C_x_diff_dt[k] - b[k];

            add_all_nonlinear_stamps(J, f, m_x);

            for(size_t j=0; j<f.size(); ++j) f[j] = -f[j];
            auto p = lu_decompose(J); Vector dx = lu_solve(J, p, f);
            double damp = 1.0, max_s = 0.0; for(double v : dx) max_s = std::max(max_s, std::abs(v));
            if (max_s > m_nr_damping) damp = m_nr_damping / max_s;
            double norm = 0.0; for (int j=0; j<m_num_nodes; ++j) { double s = damp*dx[j]; m_x[j]+=s; norm+=s*s; }
            if (sqrt(norm) < m_nr_tolerance) break;
            if (norm < norm_min) norm_min = norm;
        }
        if (iter >= m_nr_max_iter) {
            std::cerr << "too many iterations; norm=" << sqrt(norm_min) << std::endl;
        }
        m_x_prev = m_x;
#endif

        double output_sample = m_x.at(m_output_node);
        return output_sample;
    }

public:
    struct Resistor { std::string name; int n1, n2; double R; };
    struct Capacitor { std::string name; int n1, n2; double C; };
    struct VoltageSource { std::string name; int node; double dc_voltage; bool is_time_varying; };
    struct Triode { std::string name; int p_node, g_node, k_node; double mu, ex, kg1, kp, kvb, rgi; };
    struct Potentiometer { std::string name; int n1, n2, wiper; double R_total; double value; char taper; };
    struct VariableResistor { std::string name; int n1, n2; double R_max; double value; char taper; };

    std::vector<Resistor> m_resistors;
    std::vector<Capacitor> m_capacitors;
    std::vector<VoltageSource> m_v_sources;
    std::vector<Triode> m_triodes;
    std::map<std::string, Potentiometer> m_potentiometers;
    std::map<std::string, VariableResistor> m_variable_resistors;

    Matrix m_G_static, m_G_dynamic, m_C_mat;
    Vector m_x, m_x_prev, m_x_prev2;
    int m_num_nodes = 0;
    std::map<std::string, int> m_node_map;
    std::map<int, std::string> m_node_map_rev;
    bool m_is_prepared = false;
    bool m_params_dirty = true;

    int m_output_node = -1;

    double m_sample_rate, m_dt;
    const int m_nr_max_iter = 50, m_nr_damping = 1.0;
    const double m_nr_tolerance = 1e-6, G_LARGE = 1e12;

    void update_dynamic_components() {
        m_G_dynamic = m_G_static;
        for (const auto& pair : m_potentiometers) {
            const auto& p = pair.second;
            double w = std::clamp(p.value, 0.01, 0.99);
            if (p.taper == 'A' || p.taper == 'a') w = w * w; // Audio taper approximation
            double r1 = p.R_total * (1.0 - w);
            double r2 = p.R_total * (w);
            stamp_resistor(m_G_dynamic, p.n1, p.wiper, r1);
            stamp_resistor(m_G_dynamic, p.wiper, p.n2, r2);
        }
        for (const auto& pair : m_variable_resistors) {
            const auto& vr = pair.second;
            double w = std::clamp(vr.value, 0.01, 0.99);
            if (vr.taper == 'A' || vr.taper == 'a') w = w * w;
            stamp_resistor(m_G_dynamic, vr.n1, vr.n2, vr.R_max * (1.0 - w));
        }
        m_params_dirty = false;
    }

    void add_all_nonlinear_stamps(Matrix& J, Vector& f, const Vector& x) {
        for (const auto& tube : m_triodes) {
            add_single_triode_stamps_numerical(tube, J, f, x);
        }
    }

    int m_node(const std::string& name) {
        if (name == "gnd" || name == "0") return -1;
        if (m_node_map.find(name) == m_node_map.end()) return add_node(name);
        return m_node_map.at(name);
    }

    void solve_dc() {
        Vector b_dc(m_num_nodes, 0.0);
        for (const auto& vs : m_v_sources) b_dc[vs.node] += vs.dc_voltage * G_LARGE;
        m_x.assign(m_num_nodes, 0.0);
        for (int i=0; i<30; ++i) {
            Matrix J = m_G_dynamic; Vector f = matrix_vector_mult(m_G_dynamic, m_x);
            for(size_t j=0; j<f.size(); ++j) f[j] -= b_dc[j];
            add_all_nonlinear_stamps(J, f, m_x);
            for(size_t j=0; j<f.size(); ++j) f[j] = -f[j];
            auto p = lu_decompose(J); Vector dx = lu_solve(J, p, f);
            double norm=0.0; for (int j=0; j<m_num_nodes; ++j){ m_x[j]+=dx[j]; norm+=dx[j]*dx[j]; }
            if (sqrt(norm) < m_nr_tolerance) break;
        }
        m_x_prev = m_x;
    }

#if 1
    // A helper to evaluate the linear part of the system residual `f(x)`
    Vector eval_f(const Vector& x_eval, const Vector& b) {
        Vector Gx = matrix_vector_mult(m_G_dynamic, x_eval);
        Vector x_diff = x_eval; for(size_t k=0; k<x_eval.size(); ++k) x_diff[k] -= m_x_prev[k];
        Vector C_x_diff_dt = matrix_vector_mult(m_C_mat, x_diff); for(size_t k=0; k<x_eval.size(); ++k) C_x_diff_dt[k] /= m_dt;
        Vector f(m_num_nodes); for(size_t k=0; k<f.size(); ++k) f[k] = Gx[k] + C_x_diff_dt[k] - b[k];
        return f;
    }

    double vector_norm(const Vector& v) {
        double norm = 0.0; for (double val : v) norm += val * val; return sqrt(norm);
    }

    static void get_triode_currents(const Triode& t, double Vp, double Vg, double Vk, double& Ip, double& Ig) {
        const double Vpk = Vp - Vk, Vgk = Vg - Vk;
        const double kvb_vpk_sq = t.kvb + Vpk * Vpk;
        Ip = 0; Ig = 0;

        if (kvb_vpk_sq > 0) {
            const double inner_exp = t.kp * (1.0 / t.mu + Vgk / sqrt(kvb_vpk_sq));
            // The softplus function log(1+exp(x)) is naturally smooth.
            // Avoid artificial cutoffs. Check for infinity instead.
            double exp_val = exp(inner_exp);
            if (!std::isinf(exp_val)) {
                double E1 = (Vpk / t.kp) * log(1.0 + exp_val);
                if (E1 > 0) {
                    Ip = pow(E1, t.ex) / t.kg1;
                }
            } else {
                throw std::runtime_error("exp_val went InF!");
            }
        }
        if (Vgk > 0) {
            Ig = (Vgk * Vgk) / (t.rgi + Vgk);
        }
    }
#else
    static void get_triode_currents(const Triode& t, double Vp, double Vg, double Vk, double& Ip, double& Ig) {
        const double Vpk = Vp - Vk, Vgk = Vg - Vk;
        const double kvb_vpk_sq = t.kvb + Vpk * Vpk;
        Ip = 0;
        if (kvb_vpk_sq > 0) {
            const double inner_exp = t.kp * (1.0 / t.mu + Vgk / sqrt(kvb_vpk_sq));
            if (inner_exp < 100) { // log(1+exp(x)) is ~x for large x
                double E1 = (Vpk / t.kp) * log(1.0 + exp(inner_exp));
                if (E1 > 0) Ip = pow(E1, t.ex) / t.kg1;
            // } else {
            //     throw std::runtime_error("inner_exp >= 100");
            }
        }
        Ig = (Vgk > 0) ? (Vgk * Vgk) / (t.rgi + Vgk) : 0;
    }
#endif

    void add_single_triode_stamps_numerical(const Triode& t, Matrix& J, Vector& f, const Vector& x) {
        double Vp=x[t.p_node], Vg=x[t.g_node], Vk=x[t.k_node];
        double Ip0, Ig0; get_triode_currents(t, Vp, Vg, Vk, Ip0, Ig0);
        f[t.p_node] += Ip0; f[t.k_node] -= Ip0; f[t.g_node] += Ig0; f[t.k_node] -= Ig0;
        const double delta = 1e-6;
        double Ip_p, Ig_p; get_triode_currents(t, Vp+delta, Vg, Vk, Ip_p, Ig_p);
        double dIp_dVp=(Ip_p-Ip0)/delta, dIg_dVp=(Ig_p-Ig0)/delta;
        double Ip_g, Ig_g; get_triode_currents(t, Vp, Vg+delta, Vk, Ip_g, Ig_g);
        double dIp_dVg=(Ip_g-Ip0)/delta, dIg_dVg=(Ig_g-Ig0)/delta;
        double Ip_k, Ig_k; get_triode_currents(t, Vp, Vg, Vk+delta, Ip_k, Ig_k);
        double dIp_dVk=(Ip_k-Ip0)/delta, dIg_dVk=(Ig_k-Ig0)/delta;
        J[t.p_node][t.p_node]+=dIp_dVp; J[t.p_node][t.g_node]+=dIp_dVg; J[t.p_node][t.k_node]+=dIp_dVk;
        J[t.g_node][t.p_node]+=dIg_dVp; J[t.g_node][t.g_node]+=dIg_dVg; J[t.g_node][t.k_node]+=dIg_dVk;
        J[t.k_node][t.p_node]-=(dIp_dVp+dIg_dVp); J[t.k_node][t.g_node]-=(dIp_dVg+dIg_dVg); J[t.k_node][t.k_node]-=(dIp_dVk+dIg_dVk);
    }

    void stamp_resistor(Matrix& G, int n1, int n2, double R) {
        if (R <= 1e-9) R = 1e-9; double g = 1.0 / R;
        if (n1!=-1) G[n1][n1]+=g; if (n2!=-1) G[n2][n2]+=g;
        if (n1!=-1 && n2!=-1) { G[n1][n2]-=g; G[n2][n1]-=g; }
    }

    void stamp_capacitor(int n1, int n2, double C) {
        if (n1!=-1) m_C_mat[n1][n1]+=C; if (n2!=-1) m_C_mat[n2][n2]+=C;
        if (n1!=-1 && n2!=-1) { m_C_mat[n1][n2]-=C; m_C_mat[n2][n1]-=C; }
    }

    // Vector matrix_vector_mult(const Matrix& A, const Vector& v) { /* ... same ... */ }
    // Vector lu_decompose(Matrix& A) { /* ... same ... */ }
    // Vector lu_solve(const Matrix& LU, const Vector& p_in, const Vector& b) { /* ... same ... */ }

    // Re-paste the implementations to avoid snipping
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
            if (std::abs(A[i][i]) < 1e-20) continue;
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
            if (std::abs(LU[i][i]) > 1e-20) x[i] /= LU[i][i]; else x[i] = 0;
        } return x;
    }
};
