#pragma once

#include "MNASolver.h"
#include "Triode.h"
#include <vector>

// This class models the second tube stage (V1B) of the preamp.
// It is a standard common-cathode gain stage.
//
// Input:  Voltage at node N020 (from the Volume 1 pot wiper).
// Output: Voltage at node N001.
//
// SPICE components included: XV1B, R8, R7, C13, C7, R9.
// The system has 4 unknown nodes and 1 voltage source for the input.
class V1B_Stage : public MNASolver<4, 1> {
private:
    static constexpr int GND = -1;

    // A single enum for all unknowns in the solution vector 'x'.
    // The first 4 are node voltages (0-based indices).
    // The last 1 is the current through the input voltage source.
    enum Var {
        // Node Voltages (alphabetical by node name)
        V_N001 = 0, // Output node
        V_N009,     // Plate node
        V_N020,     // Grid / Input node
        V_N027,     // Cathode node
        // I_VIN at index 4 (NumNodes)
    };

    // Component values
    double R7, R8, R9;
    double C7, C13;

    // State variables for the capacitors
    std::vector<double> cap_z_state;

    // Power supply voltage connected to the plate load resistor
    const double VE = 405.0;

    void stampComponents(double inputVoltage) {
        resetMatrices();

        // --- Input Source ---
        // The input signal is the voltage at node N020.
        stampVoltageSource(V_N020, GND, 0, inputVoltage);

        // --- Static Components ---
        // Plate load resistor (R8) connected to power supply VE
        stampResistor(V_N009, GND, R8);
        b[V_N009] += VE / R8;

        // Cathode resistor (R7)
        stampResistor(V_N027, GND, R7);

        // Load resistor / next stage's grid-leak (R9)
        stampResistor(V_N001, GND, R9);

        // --- Capacitors ---
        // Cathode bypass capacitor (C13)
        stampCapacitor(V_N027, GND, C13, cap_z_state[0]);

        // Output coupling capacitor (C7)
        stampCapacitor(V_N001, V_N009, C7, cap_z_state[1]);
    }

    void stampLinear() override {
        stampResistorLinear(V_N009, GND, R8);
        b_linear[V_N009] += VE / R8; // Static contribution from power supply

        stampResistorLinear(V_N027, GND, R7);
        stampResistorLinear(V_N001, GND, R9);

        // Stamp the constant, linear part of the capacitors
        stampCapacitor_A(V_N027, GND, C13);
        stampCapacitor_A(V_N001, V_N009, C7);
    }

    void stampDynamic(double in) override {
        // Stamp the input signal source
        stampVoltageSource(V_N020, GND, 0, in);

        // Stamp only the dynamic history current of the capacitors
        stampCapacitor_b(V_N027, GND, cap_z_state[0]);    // C13
        stampCapacitor_b(V_N001, V_N009, cap_z_state[1]); // C7
    }

    void stampNonLinear(const std::array<double, NumUnknowns>& current_x) override {
        // XV1B: Plate=N009, Grid=N020, Cathode=N027
        double v_p = current_x[V_N009];
        double v_g = current_x[V_N020];
        double v_c = current_x[V_N027];

        Triode::State ts = Triode::calculate(v_p - v_c, v_g - v_c);

        // --- Plate Current ---
        double i_p_lin = ts.ip - ts.g_p * (v_p - v_c) - ts.g_g * (v_g - v_c);
        stampCurrentSource(V_N009, V_N027, i_p_lin);
        stampConductance(V_N009, V_N027, ts.g_p);

        A[V_N009][V_N020] += ts.g_g;
        A[V_N009][V_N027] -= ts.g_g;
        A[V_N027][V_N020] -= ts.g_g;
        A[V_N027][V_N027] += ts.g_g;

        // --- Grid Current ---
        double i_g_lin = ts.ig - ts.g_ig * (v_g - v_c);
        stampCurrentSource(V_N020, V_N027, i_g_lin);
        stampConductance(V_N020, V_N027, ts.g_ig);
    }

public:
    V1B_Stage() : cap_z_state(2, 0.0) {
        // Set fixed component values from SPICE netlist
        R7 = 1.5e3;
        R8 = 100e3;
        R9 = 100e3;
        C7 = 0.1e-6;
        C13 = 22e-6;
    }

    // No user-adjustable parameters in this stage

    double process(double in) {
        solveNonlinear(in);

        // Update capacitor states for the next time step
        updateCapacitorState(x[V_N027], 0, C13, cap_z_state[0]);
        updateCapacitorState(x[V_N001], x[V_N009], C7, cap_z_state[1]);

        // The output of this stage is the voltage at node N001
        return x[V_N001];
    }
};
