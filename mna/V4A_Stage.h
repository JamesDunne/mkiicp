#pragma once

#include "MNASolver.h"
#include "Triode.h"
#include <vector>

// This class correctly models Stage 6 (V4A) and its complex output coupling/mixing network.
// By solving them together, the AC and DC characteristics will match the simulation accurately.
//
// Input:  Voltage at node N030.
// Output: Voltage at the main mixer node, N002.
//
// SPICE components: XV4A, R27, R30, C29, C30, R31, C31, R32, C32, and VC power supply.
// The system has 6 unknown nodes and 2 voltage sources.
class V4A_and_Coupling : public MNASolver<6, 2> {
private:
    static constexpr int GND = -1;

    enum Var {
        V_N025 = 0, // V4A Plate
        V_N030,     // V4A Grid (Input)
        V_N034,     // V4A Cathode
        V_N010,     // Power Supply Node
        V_N026,     // Intermediate coupling node
        V_N002      // Mixer/Output node
        // V-Source Currents for Vin and Vsupply are at indices 6 and 7
    };

    // Component values
    double R27, R30, R31, R32;
    double C29, C30, C31, C32;

    // State variables for the capacitors
    std::vector<double> cap_z_state;

    const double VC = 410.0;

    void stampComponents(double inputVoltage) {
        resetMatrices();

        // --- Sources ---
        stampVoltageSource(V_N030, GND, 0, inputVoltage); // Input signal
        stampVoltageSource(V_N010, GND, 1, VC);         // Power supply

        // --- Static Components ---
        stampResistor(V_N025, V_N010, R27); // V4A Plate Load
        stampResistor(V_N034, GND, R30);    // V4A Cathode Resistor
        stampResistor(V_N002, V_N026, R31); // Mixer resistor
        stampResistor(V_N002, GND, R32);    // Mixer load resistor

        // --- Capacitors ---
        stampCapacitor(V_N034, GND, C29, cap_z_state[0]);       // V4A Cathode Bypass
        stampCapacitor(V_N026, V_N025, C30, cap_z_state[1]);   // V4A Output Coupling Cap
        stampCapacitor(V_N002, V_N026, C31, cap_z_state[2]);   // Mixer "bright" cap
        stampCapacitor(V_N002, GND, C32, cap_z_state[3]);       // Mixer load cap
    }

    void stampNonLinear(const std::array<double, NumUnknowns>& current_x) {
        // XV4A: Plate=N025, Grid=N030, Cathode=N034
        double v_p = current_x[V_N025];
        double v_g = current_x[V_N030];
        double v_c = current_x[V_N034];

        Triode::State ts = Triode::calculate(v_p - v_c, v_g - v_c);

        double i_p_lin = ts.ip - ts.g_p * (v_p - v_c) - ts.g_g * (v_g - v_c);
        stampCurrentSource(V_N025, V_N034, i_p_lin);
        stampConductance(V_N025, V_N034, ts.g_p);

        A[V_N025][V_N030] += ts.g_g;
        A[V_N025][V_N034] -= ts.g_g;
        A[V_N034][V_N030] -= ts.g_g;
        A[V_N034][V_N034] += ts.g_g;

        double i_g_lin = ts.ig - ts.g_ig * (v_g - v_c);
        stampCurrentSource(V_N030, V_N034, i_g_lin);
        stampConductance(V_N030, V_N034, ts.g_ig);
    }

public:
    V4A_and_Coupling() : cap_z_state(4, 0.0) {
        // V4A components
        R27 = 274e3;
        R30 = 3.3e3;
        C29 = 0.22e-6;
        // Coupling/Mixing components
        R31 = 220e3;
        C31 = 250e-12;
        R32 = 100e3;
        C32 = 500e-12;
        C30 = 0.047e-6;
    }

    double process(double in) {
        const int MAX_ITER = 15;
        const double CONVERGENCE_THRESH = 1e-6;

        std::array<double, NumUnknowns> current_x = x;

        for (int i = 0; i < MAX_ITER; ++i) {
            stampComponents(in);
            stampNonLinear(current_x);

            if (!lu_decompose()) { return 0.0; }

            std::array<double, NumUnknowns> next_x;
            lu_solve(next_x);

            double norm = 0.0;
            for(size_t j = 0; j < NumUnknowns; ++j) {
                double diff = next_x[j] - current_x[j];
                norm += diff * diff;
            }
            current_x = next_x;

            if (sqrt(norm) < CONVERGENCE_THRESH) { break; }
        }
        x = current_x;

        updateCapacitorState(x[V_N034], 0, C29, cap_z_state[0]);
        updateCapacitorState(x[V_N026], x[V_N025], C30, cap_z_state[1]);
        updateCapacitorState(x[V_N002], x[V_N026], C31, cap_z_state[2]);
        updateCapacitorState(x[V_N002], 0, C32, cap_z_state[3]);

        // The output of this combined stage is the voltage at the mixer node N002
        return x[V_N002];
    }
};