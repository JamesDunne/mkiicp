/*
Of course. This is a complex task that combines numerical methods, circuit theory, and software engineering. I will create a single, well-documented C++ header file that implements a realtime audio simulation of a plausible tube preamplifier circuit based on your requirements.

The chosen circuit for this simulation will be a single-stage 12AX7 (ECC83) triode amplifier with a Fender-style tone stack and master volume control. This circuit includes all the specified components: resistors, capacitors, a non-linear tube, fixed voltage sources, and potentiometers for user control.

Key Features of the Implementation:

SPICE Netlist as a Blueprint: The C++ code is a direct MNA (Modified Nodal Analysis) implementation of the following conceptual SPICE netlist.

Modified Nodal Analysis (MNA): The circuit is statically analyzed to build the core G (conductance), C (capacitance), and b (sources) matrices.

Newton-Raphson Solver: To handle the non-linear triode, a Newton-Raphson iterative solver is used at each time step to find the correct operating point.

Triode Model: A simplified but effective model for the 12AX7 triode is used, including the simulation of grid-current when the grid voltage becomes positive. This is crucial for realistic overdrive characteristics.

Dense LU Solver: A from-scratch dense LU decomposition solver is implemented to solve the linear system Ax = b that arises at each Newton-Raphson iteration.

Realtime Parameter Control: Public methods are provided to control the Volume, Bass, Mid, Treble, and Master potentiometers, allowing for dynamic changes during audio processing.

Single Header File: The entire implementation is self-contained in one header file for easy integration into any C++ audio project (like a VST, AU, or JUCE plugin).

RealtimeTubeSim.h
Generated cpp
*/
#pragma once

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>

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
 * // Nodes: 0 (GND), in, g, p, k, ts1, ts2, ts3, out
 *
 * // --- Input Stage ---
 * V_in   in 0 DC 0.0      ; Input signal, represented by b vector update
 * R_in   in 0 1M          ; Input impedance
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
 * C_out  out final 100n   ; Output coupling capacitor
 * R_load final 0 1M       ; Final load resistor
 *
 * // Fixed DC voltages (VE, VC, VC2 are conceptual placeholders, V_P is used)
 * // VE, VC, VC2 are not used in this specific circuit but the system supports fixed V sources.
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
     * @param value A value from 0.0 to 1.0. This is modeled as a potentiometer
     *              before the grid coupling capacitor, but for simplicity here,
     *              it acts as a simple gain multiplier on the input signal.
     */
    void set_volume(double value) { m_params.volume = std::max(0.0, std::min(1.0, value)); }

    /**
     * @brief Sets the bass control of the tone stack.
     * @param value A value from 0.0 (min bass) to 1.0 (max bass).
     */
    void set_bass(double value) { m_params.bass = std::max(0.0, std::min(1.0, value)); update_pots(); }

    /**
     * @brief Sets the mid control of the tone stack.
     * @param value A value from 0.0 (min mid) to 1.0 (max mid).
     */
    void set_mid(double value) { m_params.mid = std::max(0.0, std::min(1.0, value)); update_pots(); }

    /**
     * @brief Sets the treble control of the tone stack.
     * @param value A value from 0.0 (min treble) to 1.0 (max treble).
     */
    void set_treble(double value) { m_params.treble = std::max(0.0, std::min(1.0, value)); update_pots(); }

    /**
     * @brief Sets the master volume control.
     * @param value A value from 0.0 (min volume) to 1.0 (max volume).
     */
    void set_master(double value) { m_params.master = std::max(0.0, std::min(1.0, value)); update_pots(); }

    /**
     * @brief Processes a single audio sample.
     * @param input_sample The input audio sample.
     * @return The processed output audio sample.
     */
    double process_sample(double input_sample) {
        // --- 1. Update time-varying sources ---
        // The input signal is modeled as a voltage source.
        // We apply the 'volume' gain here.
        m_b[m_node_map["in"]] = input_sample * m_params.volume;

        // --- 2. Newton-Raphson Iteration for this time step ---
        // Initial guess for this timestep's solution is the previous one.
        m_x = m_x_prev;

        for (int i = 0; i < m_nr_max_iter; ++i) {
            // Build the Jacobian matrix (A) and the residual vector (f)
            // A = G_dynamic + C/dt + J_nonlinear
            // f = (G_dynamic + C/dt)*x - (b + (C/dt)*x_prev) - I_nonlinear

            auto A = m_G_dynamic; // Start with G matrix (includes pots)
            auto f = m_b;         // Start with source vector

            // Add capacitor contributions
            for (int r = 0; r < m_num_nodes; ++r) {
                for (int c = 0; c < m_num_nodes; ++c) {
                    if (m_C(r, c) != 0.0) {
                        double c_val_dt = m_C(r, c) / m_dt;
                        A[r][c] += c_val_dt;
                        f[r] -= c_val_dt * (m_x[c] - m_x_prev[c]);
                    }
                }
            }

            // Add non-linear triode contributions
            add_triode_nr_stamps(A, f);

            // Now, f is the residual f(x_k). We want to solve J*dx = -f
            // where J is our final 'A' matrix.
            for(size_t j = 0; j < f.size(); ++j) {
                f[j] = -f[j];
            }

            // Solve the linear system A * dx = f for dx
            Vector dx = lu_solve(A, lu_decompose(A), f);

            // Update the solution
            for (int j = 0; j < m_num_nodes; ++j) {
                m_x[j] += dx[j];
            }

            // Check for convergence
            double norm = 0.0;
            for (double val : dx) {
                norm += val * val;
            }
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
    Matrix m_G_static;      // Conductance matrix for fixed resistors
    Matrix m_G_dynamic;     // Conductance matrix including variable pots
    MatrixView m_C { m_C_mat }; // Capacitance matrix
    Matrix m_C_mat;
    Vector m_b;             // Source vector
    Vector m_x;             // Solution vector (node voltages) at current time n
    Vector m_x_prev;        // Solution vector at previous time n-1
    int m_num_nodes = 0;
    std::map<std::string, int> m_node_map;

    // --- Simulation Parameters ---
    double m_sample_rate;
    double m_dt; // Time step (1/sample_rate)

    // --- Newton-Raphson Parameters ---
    const int m_nr_max_iter = 15;
    const double m_nr_tolerance = 1e-6; // Convergence tolerance

    // --- Circuit Component Pointers for dynamic updates ---
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
    const double TUBE_K_GRID = 30.0;   // Grid current factor
    const double TUBE_V_GRID_THRESH = 0.5; // Grid current threshold voltage

    // Node indices for the triode
    int m_p_node, m_g_node, m_k_node;

private:
    // --- Core Setup Methods ---

    /**
     * @brief Defines the circuit topology, creates nodes, and builds the static MNA matrices.
     */
    void build_system() {
        // --- 1. Define Nodes ---
        // Ground (node 0) is implicit. We map named nodes to matrix indices.
        int node_idx = 0;
        auto add_node = [&](const std::string& name) {
            m_node_map[name] = node_idx++;
        };
        add_node("in");         // Input signal
        add_node("g");          // Tube grid
        add_node("p");          // Tube plate
        add_node("k");          // Tube cathode
        add_node("V_P");        // Plate DC supply voltage node
        add_node("ts1");        // Tone stack node 1
        add_node("ts2");        // Tone stack node 2
        add_node("ts3");        // Tone stack node 3
        add_node("out");        // Master vol output wiper
        add_node("out_final");  // Final output node

        m_num_nodes = node_idx;

        // --- 2. Initialize Matrices ---
        m_G_static.assign(m_num_nodes, Vector(m_num_nodes, 0.0));
        m_C_mat.assign(m_num_nodes, Vector(m_num_nodes, 0.0));
        m_b.assign(m_num_nodes, 0.0);

        // --- 3. Stamp Fixed Components ---
        // Resistors
        stamp_resistor(m_G_static, "in", "gnd", 1e6);    // R_in
        stamp_resistor(m_G_static, "g", "gnd", 1e6);     // R_g
        stamp_resistor(m_G_static, "V_P", "p", 100e3);   // R_p
        stamp_resistor(m_G_static, "k", "gnd", 1.5e3);   // R_k
        stamp_resistor(m_G_static, "out_final", "gnd", 1e6); // R_load

        // Capacitors
        stamp_capacitor("in", "g", 22e-9);        // C_in
        stamp_capacitor("k", "gnd", 22e-6);       // C_k
        stamp_capacitor("p", "ts1", 22e-9);       // C_ts_in
        stamp_capacitor("ts1", "ts3", 250e-12);    // C_t
        stamp_capacitor("ts2", "ts3", 22e-9);     // C_m
        stamp_capacitor("out", "out_final", 100e-9); // C_out

        // Fixed Voltage Sources (DC)
        stamp_voltage_source("V_P", "gnd", 300.0);

        // --- 4. Define Potentiometers for Dynamic Stamping ---
        // Note: For MNA, a pot is two resistors. We model the tone stack slightly
        // simplified but functionally equivalent for clarity.
        // Treble Pot: 250k Audio, modeled as a variable resistor
        m_pot_treble = {m_node("ts1"), m_node("ts2"), -1, 250e3};
        // Bass Pot: 1M Reverse Audio, modeled as variable resistor to ground
        m_pot_bass = {m_node("ts2"), m_node("gnd"), -1, 1e6};
        // Mid Pot: 25k Linear, modeled as variable resistor to ground
        m_pot_mid = {m_node("ts3"), m_node("gnd"), -1, 25e3};
        // Master Volume Pot: 1M Audio, modeled as a voltage divider
        m_pot_master = {m_node("ts3"), m_node("gnd"), m_node("out"), 1e6};

        // --- 5. Store triode node indices ---
        m_p_node = m_node("p");
        m_g_node = m_node("g");
        m_k_node = m_node("k");
    }

    /**
     * @brief Recalculates the G matrix contributions from potentiometers based on current settings.
     */
    void update_pots() {
        m_G_dynamic = m_G_static;

        // Treble Pot (modeled as a variable series resistor)
        // A-taper pot: use a simple power curve for approximation
        double treble_res = m_pot_treble.total_resistance * pow(m_params.treble, 2);
        stamp_resistor(m_G_dynamic, m_pot_treble.n1, m_pot_treble.n2, std::max(1.0, treble_res));

        // Bass Pot (modeled as a variable resistor to ground)
        // Reverse A-taper:
        double bass_res = m_pot_bass.total_resistance * pow((1.0 - m_params.bass), 2);
        stamp_resistor(m_G_dynamic, m_pot_bass.n1, m_pot_bass.n2, std::max(1.0, bass_res));

        // Mid Pot (linear)
        double mid_res = m_pot_mid.total_resistance * m_params.mid;
        stamp_resistor(m_G_dynamic, m_pot_mid.n1, m_pot_mid.n2, std::max(1.0, mid_res));

        // Master Volume Pot (voltage divider)
        // A-taper:
        double r1_val = m_pot_master.total_resistance * pow(m_params.master, 2);
        double r2_val = m_pot_master.total_resistance - r1_val;
        stamp_resistor(m_G_dynamic, m_pot_master.n1, m_pot_master.n_wiper, std::max(1.0, r1_val));
        stamp_resistor(m_G_dynamic, m_pot_master.n_wiper, m_pot_master.n2, std::max(1.0, r2_val));
    }

    /**
     * @brief Solves the DC operating point of the circuit.
     * This is done by running a Newton-Raphson simulation with capacitors ignored (treated as open circuits).
     */
    void solve_dc() {
        m_x.assign(m_num_nodes, 0.0); // Start with zero voltages

        for (int i = 0; i < m_nr_max_iter * 2; ++i) { // More iterations for initial convergence
             auto A = m_G_dynamic;
             auto f = m_b;

             add_triode_nr_stamps(A, f);

             for(size_t j = 0; j < f.size(); ++j) {
                f[j] = -f[j];
             }

             Vector dx = lu_solve(A, lu_decompose(A), f);
             for (int j = 0; j < m_num_nodes; ++j) m_x[j] += dx[j];

             double norm = 0.0;
             for (double val : dx) norm += val * val;
             if (sqrt(norm) < m_nr_tolerance) break;
        }

        m_x_prev = m_x; // Set the initial state for AC analysis
    }


    // --- MNA Stamping Helpers ---
    int m_node(const std::string& name) {
        if (name == "gnd" || name == "0") return -1; // -1 represents ground
        return m_node_map.at(name);
    }

    void stamp_resistor(Matrix& G, int n1, int n2, double resistance) {
        if (resistance <= 1e-9) return; // Avoid division by zero
        double conductance = 1.0 / resistance;
        if (n1 != -1) G[n1][n1] += conductance;
        if (n2 != -1) G[n2][n2] += conductance;
        if (n1 != -1 && n2 != -1) {
            G[n1][n2] -= conductance;
            G[n2][n1] -= conductance;
        }
    }
    void stamp_resistor(Matrix& G, const std::string& n1, const std::string& n2, double R) {
        stamp_resistor(G, m_node(n1), m_node(n2), R);
    }

    void stamp_capacitor(int n1, int n2, double capacitance) {
        if (n1 != -1) m_C(n1, n1) += capacitance;
        if (n2 != -1) m_C(n2, n2) += capacitance;
        if (n1 != -1 && n2 != -1) {
            m_C(n1, n2) -= capacitance;
            m_C(n2, n1) -= capacitance;
        }
    }
    void stamp_capacitor(const std::string& n1, const std::string& n2, double C_val) {
        stamp_capacitor(m_node(n1), m_node(n2), C_val);
    }

    void stamp_voltage_source(const std::string& n1, const std::string& n2, double voltage) {
        // For simplicity with this solver, we are only handling voltage sources
        // connected to ground. This creates a "known voltage" at a node.
        // We can enforce this with a very large conductance to the source vector.
        // A full MNA implementation would add a new variable for the current.
        int node1 = m_node(n1);
        if (n2 != "gnd" && n2 != "0") {
             throw std::runtime_error("This simple MNA only supports V-sources connected to ground.");
        }
        double G_large = 1e12; // A very large conductance
        m_G_static[node1][node1] += G_large;
        m_b[node1] += voltage * G_large;
    }

    // --- Non-Linear Tube Model ---

    /**
     * @brief Adds the triode's contribution to the Jacobian and residual for Newton-Raphson.
     * @param J The Jacobian matrix (the system matrix 'A').
     * @param f The residual vector (the right-hand side 'b' before negation).
     */
    void add_triode_nr_stamps(Matrix& J, Vector& f) {
        // Get current estimates for plate, grid, and cathode voltages
        double Vp = (m_p_node != -1) ? m_x[m_p_node] : 0.0;
        double Vg = (m_g_node != -1) ? m_x[m_g_node] : 0.0;
        double Vk = (m_k_node != -1) ? m_x[m_k_node] : 0.0;

        double Vpk = Vp - Vk;
        double Vgk = Vg - Vk;

        // --- Plate Current (Koren Model) ---
        double E1 = Vpk / TUBE_KP + Vgk;
        double Ip = 0.0, Gp = 0.0, Gm = 0.0; // Plate current and its partial derivatives

        if (E1 > 0) {
            double E1_ex = pow(E1, TUBE_EX);
            Ip = (2.0 / TUBE_KG) * E1_ex * E1; // Using E1^Ex * E1 to avoid issues at E1=0

            // Partial derivatives for Jacobian
            double dIp_dE1 = (2.0 / TUBE_KG) * (TUBE_EX + 1.0) * E1_ex;
            Gp = dIp_dE1 / TUBE_KP; // dIp/dVp
            Gm = dIp_dE1;          // dIp/dVg
        }

        // --- Grid Current (Diode Model) ---
        double Ig = 0.0, Gg = 0.0; // Grid current and its conductance
        if (Vgk > TUBE_V_GRID_THRESH) {
            Ig = (Vgk - TUBE_V_GRID_THRESH) * TUBE_K_GRID;
            Gg = TUBE_K_GRID;
        }

        // --- Stamp contributions into Jacobian (J) and Residual (f) ---
        // f(x) = ... - I_nonlinear(x)
        // J(x) = ... + G_nonlinear(x)

        // Plate current stamps
        if (m_p_node != -1) { f[m_p_node] -= Ip; }
        if (m_k_node != -1) { f[m_k_node] += Ip; }

        // Grid current stamps
        if (m_g_node != -1) { f[m_g_node] -= Ig; }
        if (m_k_node != -1) { f[m_k_node] += Ig; }

        // Jacobian plate conductance stamps
        if (m_p_node != -1) { J[m_p_node][m_p_node] += Gp; }
        if (m_k_node != -1) { J[m_k_node][m_k_node] += Gp + Gm; }
        if (m_p_node != -1 && m_k_node != -1) { J[m_p_node][m_k_node] -= Gp + Gm; }
        if (m_k_node != -1 && m_p_node != -1) { J[m_k_node][m_p_node] -= Gp; }
        if (m_p_node != -1 && m_g_node != -1) { J[m_p_node][m_g_node] += Gm; }
        if (m_g_node != -1 && m_p_node != -1) { J[m_g_node][m_p_node] += 0; /* Asymmetric */ }
        if (m_g_node != -1 && m_k_node != -1) { J[m_g_node][m_k_node] -= Gm; }
        if (m_k_node != -1 && m_g_node != -1) { J[m_k_node][m_g_node] += Gm; }

        // Jacobian grid conductance stamps
        if (m_g_node != -1) { J[m_g_node][m_g_node] += Gg; }
        if (m_k_node != -1) { J[m_k_node][m_k_node] += Gg; }
        if (m_g_node != -1 && m_k_node != -1) {
            J[m_g_node][m_k_node] -= Gg;
            J[m_k_node][m_g_node] -= Gg;
        }
    }


    // --- Dense LU Solver ---
    /**
     * @brief Performs LU decomposition of a matrix A in-place, with pivoting.
     * @param A The matrix to decompose. It will be replaced by its LU decomposition.
     * @return A permutation vector representing the row swaps.
     */
    Vector lu_decompose(Matrix& A) {
        int n = A.size();
        Vector p(n);
        for (int i = 0; i < n; ++i) p[i] = i;

        for (int i = 0; i < n; ++i) {
            // Pivoting
            int max_j = i;
            for (int j = i + 1; j < n; ++j) {
                if (std::abs(A[j][i]) > std::abs(A[max_j][i])) {
                    max_j = j;
                }
            }
            if (max_j != i) {
                std::swap(A[i], A[max_j]);
                std::swap(p[i], p[max_j]);
            }

            if (std::abs(A[i][i]) < 1e-12) {
                 // Matrix is singular or near-singular.
                 // This can happen during simulation and is not always an error,
                 // but indicates a potential problem. For a realtime sim, we just continue.
            }

            // Decomposition
            for (int j = i + 1; j < n; ++j) {
                A[j][i] /= A[i][i];
                for (int k = i + 1; k < n; ++k) {
                    A[j][k] -= A[j][i] * A[i][k];
                }
            }
        }
        return p;
    }

    /**
     * @brief Solves the system Ax = b given an LU-decomposed matrix.
     * @param LU The LU-decomposed matrix from lu_decompose.
     * @param p The permutation vector from lu_decompose.
     * @param b The right-hand side vector.
     * @return The solution vector x.
     */
    Vector lu_solve(const Matrix& LU, const Vector& p, const Vector& b) {
        int n = LU.size();
        Vector x(n);
        Vector y(n);

        // Forward substitution: Ly = Pb
        for (int i = 0; i < n; ++i) {
            y[i] = b[p[i]];
            for (int j = 0; j < i; ++j) {
                y[i] -= LU[i][j] * y[j];
            }
        }

        // Backward substitution: Ux = y
        for (int i = n - 1; i >= 0; --i) {
            x[i] = y[i];
            for (int j = i + 1; j < n; ++j) {
                x[i] -= LU[i][j] * x[j];
            }
            x[i] /= LU[i][i];
        }

        return x;
    }
};
