#pragma once

#include "MNASolver.h"
#include "Triode.h"
#include <vector>

// This class correctly models Stage 4 (V3B) and its subsequent load, Stage 5 (LeadPathCoupling).
// By solving them together, the AC and DC characteristics will match the simulation accurately.
//
// Input:  Voltage at node N029.
// Output: Voltage at node N030.
//
// SPICE components: XV3B, R26, R23, C23, C22, C25, R25, R24, C24, and VC power supply.
// The system has 6 unknown nodes and 2 voltage sources.
class V3B_and_Coupling : public MNASolver<6, 2> {
private:
    static constexpr int GND = -1;

    enum Var {
        V_N017 = 0, // V3B Plate
        V_N029,     // V3B Grid (Input)
        V_N035,     // V3B Cathode
        V_N010,     // Power Supply Node
        V_N018,     // Coupling Network intermediate node
        V_N030      // Coupling Network Output node
        // V-Source Currents at indices 6 and 7
    };

    // Component values
    double R23, R24, R25, R26;
    double C22, C23, C24, C25;

    // State variables for the capacitors
    std::vector<double> cap_z_state;

    const double VC = 410.0;

    void stampComponents(double inputVoltage) {
        resetMatrices();

        // --- Sources ---
        stampVoltageSource(V_N029, GND, 0, inputVoltage);
        stampVoltageSource(V_N010, GND, 1, VC);

        // --- Static Components ---
        stampResistor(V_N017, V_N010, R26); // V3B Plate Load
        stampResistor(V_N035, GND, R23);    // V3B Cathode Resistor
        stampResistor(V_N030, V_N018, R25); // Coupling Resistor
        stampResistor(V_N030, GND, R24);    // Coupling Load Resistor

        // --- Capacitors ---
        stampCapacitor(V_N035, GND, C23, cap_z_state[0]);       // V3B Cathode Bypass
        stampCapacitor(V_N035, V_N029, C22, cap_z_state[1]);   // V3B High-cut
        stampCapacitor(V_N018, V_N017, C25, cap_z_state[2]);   // Output Coupling Cap
        stampCapacitor(V_N030, GND, C24, cap_z_state[3]);       // Coupling Low-pass Cap
    }

    void stampLinear() override {
        // --- Sources ---
        stampVoltageSource_A(V_N029, GND, 0);
        stampVoltageSource_A(V_N010, GND, 1);

        // --- Static Components ---
        stampResistorLinear(V_N017, V_N010, R26); // V3B Plate Load
        stampResistorLinear(V_N035, GND, R23);    // V3B Cathode Resistor
        stampResistorLinear(V_N030, V_N018, R25); // Coupling Resistor
        stampResistorLinear(V_N030, GND, R24);    // Coupling Load Resistor

        // --- Capacitors ---
        stampCapacitor_A(V_N035, GND, C23);       // V3B Cathode Bypass
        stampCapacitor_A(V_N035, V_N029, C22);    // V3B High-cut
        stampCapacitor_A(V_N018, V_N017, C25);    // Output Coupling Cap
        stampCapacitor_A(V_N030, GND, C24);       // Coupling Low-pass Cap
    }

    void stampDynamic(double in) override {
        // --- Sources ---
        stampVoltageSource_b(0, in);
        stampVoltageSource_b(1, VC);

        // --- Capacitors ---
        stampCapacitor_b(V_N035, GND, cap_z_state[0]);       // C23 - V3B Cathode Bypass
        stampCapacitor_b(V_N035, V_N029, cap_z_state[1]);    // C22 - V3B High-cut
        stampCapacitor_b(V_N018, V_N017, cap_z_state[2]);    // C25 - Output Coupling Cap
        stampCapacitor_b(V_N030, GND, cap_z_state[3]);       // C24 - Coupling Low-pass Cap
    }

    void stampNonLinear(const std::array<double, NumUnknowns>& current_x) override {
        // XV3B: Plate=N017, Grid=N029, Cathode=N035
        double v_p = current_x[V_N017];
        double v_g = current_x[V_N029];
        double v_c = current_x[V_N035];

        Triode::State ts = Triode::calculate(v_p - v_c, v_g - v_c);

        double i_p_lin = ts.ip - ts.g_p * (v_p - v_c) - ts.g_g * (v_g - v_c);
        stampCurrentSource(V_N017, V_N035, i_p_lin);
        stampConductance(V_N017, V_N035, ts.g_p);

        A[V_N017][V_N029] += ts.g_g;
        A[V_N017][V_N035] -= ts.g_g;
        A[V_N035][V_N029] -= ts.g_g;
        A[V_N035][V_N035] += ts.g_g;

        double i_g_lin = ts.ig - ts.g_ig * (v_g - v_c);
        stampCurrentSource(V_N029, V_N035, i_g_lin);
        stampConductance(V_N029, V_N035, ts.g_ig);
    }

public:
    V3B_and_Coupling() : cap_z_state(4, 0.0) {
        // V3B components
        R23 = 1.5e3;
        R26 = 82e3;
        C22 = 120e-12;
        C23 = 2.2e-6;
        // Coupling components
        R24 = 68e3;
        R25 = 270e3;
        C24 = 1000e-12;
        C25 = 0.022e-6;
    }

    double process(double in) {
        // solveNonlinear(in);
        solveNonlinear_Simplified(in);

        updateCapacitorState(x[V_N035], 0, C23, cap_z_state[0]);
        updateCapacitorState(x[V_N035], x[V_N029], C22, cap_z_state[1]);
        updateCapacitorState(x[V_N018], x[V_N017], C25, cap_z_state[2]);
        updateCapacitorState(x[V_N030], 0, C24, cap_z_state[3]);

        // The output of the combined stage is the voltage at node N030
        return x[V_N030];
    }
};
