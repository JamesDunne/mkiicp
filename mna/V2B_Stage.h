#pragma once

#include "MNASolver.h"
#include "Triode.h"
#include <vector>

// This is the new, simplified V2B stage. It acts as a clean gain stage
// after the main Lead/Mixer block, with its own fully-modeled output load.
//
// Input:  Voltage from the Mixer output (node N002).
// Output: Voltage at the grid of the final stage (node N023).
//
// The system has 8 unknown nodes and 2 voltage sources.
class V2B_Stage : public MNASolver<8, 2> {
private:
    static constexpr int GND = -1;
    enum Var {
        V_N002 = 0, // V2B Grid (Input)
        V_N036,     // V2B Cathode
        V_N021,     // V2B Plate
        V_N006,     // Power Supply Node (VC2)
        V_P001,     // Output coupling intermediate node
        V_N022,     // Intermediate node in load network
        V_N032,     // New node in the final load
        V_N023      // New final output node (V2A Grid)
        // V-sources at indices 8 and 9
    };

    double R12, R13, R16, R46, R101, R102, R103, R105;
    double C9;
    double cap_z_state;
    const double VC2 = 410.0;

    void stampLinear() override {
        stampVoltageSource_A(V_N002, GND, 0);
        stampVoltageSource_A(V_N006, GND, 1);
        b_linear[NumNodes + 1] = VC2;
        stampResistorLinear(V_N021, V_N006, R13);
        stampResistorLinear(V_N036, GND, R16);
        stampResistorLinear(V_N022, V_P001, R105);
        stampResistorLinear(V_N022, GND, R46);
        stampResistorLinear(V_N022, V_N032, R102);
        stampResistorLinear(V_N032, GND, R101);
        stampResistorLinear(V_N023, V_N032, R12);
        stampResistorLinear(V_N023, GND, R103);
        stampCapacitor_A(V_P001, V_N021, C9);
    }

    void stampDynamic(double in) override {
        stampVoltageSource_b(0, in);
        stampCapacitor_b(V_P001, V_N021, cap_z_state);
    }

    void stampNonLinear(const std::array<double, NumUnknowns>& x) override {
        double v_p = x[V_N021]-x[V_N036];
        double v_g = x[V_N002]-x[V_N036];
        Triode::State ts = Triode::calculate(v_p, v_g);
        stampConductance(V_N021, V_N036, ts.g_p);
        A[V_N021][V_N002] += ts.g_g; A[V_N021][V_N036] -= ts.g_g;
        A[V_N036][V_N002] -= ts.g_g; A[V_N036][V_N036] += ts.g_g;
        stampConductance(V_N002, V_N036, ts.g_ig);
        double i_p_lin = ts.ip - ts.g_p*v_p - ts.g_g*v_g;
        stampCurrentSource(V_N021, V_N036, i_p_lin);
        double i_g_lin = ts.ig - ts.g_ig*v_g;
        stampCurrentSource(V_N002, V_N036, i_g_lin);
    }

public:
    V2B_Stage() : cap_z_state(0.0) {
        R12=2.2e3; R13=100e3; R16=1.5e3; R46=47e3;
        R101=4.7e3; R102=150e3; R103=47e3; R105=47e3;
        C9=0.047e-6;
    }

    double process(double in) {
        solveNonlinear_Adaptive(in);
        updateCapacitorState(x[V_P001], x[V_N021], C9, cap_z_state);
        return x[V_N023];
    }
};
