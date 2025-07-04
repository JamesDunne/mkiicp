#pragma once

#include "MNASolver.h"
#include "Triode.h"
#include <vector>

// This "super-stage" combines V3B with its full input/output coupling networks.
//
// THIS IS THE CORRECTED VERSION: It properly models the input signal from N001
// as a dedicated voltage source, creating a well-defined and solvable system.
//
// Replaces: Coupling1, V3B_Stage, and LeadPathCoupling.
//
// Input:  Voltage from the previous stage.
// Output: Voltage at node N030 (to V4A_Stage).
//
// The system now has 9 unknown nodes and 2 voltage sources.
class V3B_and_IO_Coupling : public MNASolver<9, 2> {
private:
    static constexpr int GND = -1;
    enum Var {
        V_N001 = 0, // Input signal node
        V_N011,     // After input cap C21
        V_N024,     // Input network (pot top)
        V_N029,     // V3B Grid
        V_N017,     // V3B Plate
        V_N035,     // V3B Cathode
        V_N010,     // Power Supply Node
        V_N018,     // Output network
        V_N030      // Final Output
        // V-sources at indices 9 and 10
    };

    double R21, R22, R23, R24, R25, R26;
    double RA_LEAD_DRIVE, RC_LEAD_DRIVE;
    double C21, C22, C23, C24, C25;

    std::vector<double> cap_z_state;
    const double VC = 410.0;

    void stampLinear() override {
        // --- Sources Topology ---
        stampVoltageSource_A(V_N001, GND, 0); // Input source
        stampVoltageSource_A(V_N010, GND, 1); // Power supply
        b_linear[NumNodes + 1] = VC;

        // --- Resistors ---
        stampResistorLinear(V_N011, V_N024, R21);
        stampResistorLinear(V_N024, V_N029, RA_LEAD_DRIVE);
        stampResistorLinear(V_N029, GND, RC_LEAD_DRIVE);
        stampResistorLinear(V_N029, GND, R22);
        stampResistorLinear(V_N017, V_N010, R26);
        stampResistorLinear(V_N035, GND, R23);
        stampResistorLinear(V_N030, V_N018, R25);
        stampResistorLinear(V_N030, GND, R24);

        // --- Capacitor Conductances ---
        stampCapacitor_A(V_N011, V_N001, C21); // Correctly stamped between two nodes
        stampCapacitor_A(V_N035, GND, C23);
        stampCapacitor_A(V_N035, V_N029, C22);
        stampCapacitor_A(V_N018, V_N017, C25);
        stampCapacitor_A(V_N030, GND, C24);
    }

    void stampDynamic(double in) override {
        // Set the dynamic input voltage value
        stampVoltageSource_b(0, in);

        // History currents
        stampCapacitor_b(V_N011, V_N001, cap_z_state[0]);
        stampCapacitor_b(V_N035, GND, cap_z_state[1]);
        stampCapacitor_b(V_N035, V_N029, cap_z_state[2]);
        stampCapacitor_b(V_N018, V_N017, cap_z_state[3]);
        stampCapacitor_b(V_N030, GND, cap_z_state[4]);
    }

    void stampNonLinear(const std::array<double, NumUnknowns>& current_x) override {
        double v_p=current_x[V_N017], v_g=current_x[V_N029], v_c=current_x[V_N035];
        Triode::State ts=Triode::calculate(v_p-v_c, v_g-v_c);
        stampConductance(V_N017, V_N035, ts.g_p);
        A[V_N017][V_N029] += ts.g_g; A[V_N017][V_N035] -= ts.g_g;
        A[V_N035][V_N029] -= ts.g_g; A[V_N035][V_N035] += ts.g_g;
        stampConductance(V_N029, V_N035, ts.g_ig);
        stampNonLinear_b_only(current_x);
    }

    void stampNonLinear_b_only(const std::array<double, NumUnknowns>& current_x) override {
        double v_p=current_x[V_N017], v_g=current_x[V_N029], v_c=current_x[V_N035];
        Triode::State ts=Triode::calculate(v_p-v_c, v_g-v_c);
        double i_p_lin=ts.ip-ts.g_p*(v_p-v_c)-ts.g_g*(v_g-v_c);
        stampCurrentSource(V_N017, V_N035, i_p_lin);
        double i_g_lin=ts.ig-ts.g_ig*(v_g-v_c);
        stampCurrentSource(V_N029, V_N035, i_g_lin);
    }

public:
    V3B_and_IO_Coupling() : cap_z_state(5, 0.0) {
        C21 = 0.02e-6; R21 = 680e3; R22 = 475e3;
        R23 = 1.5e3; R26 = 82e3; C22 = 120e-12; C23 = 2.2e-6;
        R24 = 68e3; R25 = 270e3; C24 = 1000e-12; C25 = 0.022e-6;
        setGain(0.75);
    }

    void setGain(double val) {
        double p = val * val;
        RA_LEAD_DRIVE = 1e6 * (1.0 - p);
        RC_LEAD_DRIVE = 1e6 * p;
        if (RA_LEAD_DRIVE < 1) RA_LEAD_DRIVE = 1;
        if (RC_LEAD_DRIVE < 1) RC_LEAD_DRIVE = 1;
        setDirty();
    }

    double process(double in) {
        solveNonlinear_Adaptive(in);

        // Update capacitor states
        updateCapacitorState(x[V_N011], x[V_N001], C21, cap_z_state[0]);
        updateCapacitorState(x[V_N035], 0, C23, cap_z_state[1]);
        updateCapacitorState(x[V_N035], x[V_N029], C22, cap_z_state[2]);
        updateCapacitorState(x[V_N018], x[V_N017], C25, cap_z_state[3]);
        updateCapacitorState(x[V_N030], 0, C24, cap_z_state[4]);

        return x[V_N030];
    }
};
