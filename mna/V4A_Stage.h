#pragma once

#include "MNASolver.h"
#include "Triode.h"
#include <vector>

// This class models Stage 6 (V4A) and its full output load, including the
// components that form the lead channel's contribution to the mixer.
// This is done to ensure V4A has the correct AC load for proper gain and frequency response.
//
// Input:  Voltage at node N030.
// Output: Voltage at the main mixer node, N002.
//
// SPICE components: XV4A, R27, R30, C29, C30, R31, C31, R32, C32, and VC power supply.
// The system has 6 unknown nodes and 2 voltage sources.
class V4A_and_MixerLoad : public MNASolver<6, 2> {
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

    std::vector<double> cap_z_state;
    const double VC = 410.0;

    void stampComponents(double inputVoltage) {
        resetMatrices();

        stampVoltageSource(V_N030, GND, 0, inputVoltage);
        stampVoltageSource(V_N010, GND, 1, VC);

        stampResistor(V_N025, V_N010, R27);
        stampResistor(V_N034, GND, R30);
        stampResistor(V_N002, V_N026, R31);
        stampResistor(V_N002, GND, R32);

        stampCapacitor(V_N034, GND, C29, cap_z_state[0]);
        stampCapacitor(V_N026, V_N025, C30, cap_z_state[1]);
        stampCapacitor(V_N002, V_N026, C31, cap_z_state[2]);
        stampCapacitor(V_N002, GND, C32, cap_z_state[3]);
    }

    void stampLinear() override {
        stampVoltageSource_A(V_N030, GND, 0); // inputVoltage
        stampVoltageSource_A(V_N010, GND, 1); // VC

        stampResistorLinear(V_N025, V_N010, R27);
        stampResistorLinear(V_N034, GND, R30);
        stampResistorLinear(V_N002, V_N026, R31);
        stampResistorLinear(V_N002, GND, R32);

        stampCapacitor_A(V_N034, GND, C29);
        stampCapacitor_A(V_N026, V_N025, C30);
        stampCapacitor_A(V_N002, V_N026, C31);
        stampCapacitor_A(V_N002, GND, C32);
    }

    void stampDynamic(double in) override {
        stampVoltageSource_b(0, in);
        stampVoltageSource_b(1, VC);

        stampCapacitor_b(V_N034, GND, cap_z_state[0]);    // C29
        stampCapacitor_b(V_N026, V_N025, cap_z_state[1]); // C30
        stampCapacitor_b(V_N002, V_N026, cap_z_state[2]); // C31
        stampCapacitor_b(V_N002, GND, cap_z_state[3]);    // C32
    }

    void stampNonLinear(const std::array<double, NumUnknowns>& current_x) override {
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
    V4A_and_MixerLoad() : cap_z_state(4, 0.0) {
        R27 = 274e3; R30 = 3.3e3; C29 = 0.22e-6;
        R31 = 220e3; C31 = 250e-12; R32 = 100e3;
        C32 = 500e-12; C30 = 0.047e-6;
    }

    double process(double in) {
        solveNonlinear(in);

        updateCapacitorState(x[V_N034], 0, C29, cap_z_state[0]);
        updateCapacitorState(x[V_N026], x[V_N025], C30, cap_z_state[1]);
        updateCapacitorState(x[V_N002], x[V_N026], C31, cap_z_state[2]);
        updateCapacitorState(x[V_N002], 0, C32, cap_z_state[3]);

        return x[V_N002];
    }
};