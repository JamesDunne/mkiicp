#pragma once

#include "MNASolver.h"
#include "Triode.h"
#include <vector>

// This "super-stage" models the entire lead channel and mixer network together.
// It combines the lead drive pot, V3B, V4A, and the mixer components into a
// single, highly-accurate MNA system. This solves impedance, loading, and
// divergence issues caused by splitting them apart.
//
// Input:  Voltage from V1B stage (which is also the rhythm channel input).
// Output: Voltage at the mixer output / V2B grid (node N002).
//
// Replaces: V3B_and_IO_Coupling and parts of the old MixerAndV2B.
// The system has 13 unknown nodes and 2 voltage sources.
class LeadAndMixerStage : public MNASolver<13, 2> {
private:
    static constexpr int GND = -1;
    enum Var {
        V_N001 = 0, // Input from V1B / Rhythm Channel
        V_N011,     // Lead path after C21
        V_N024,     // Lead pot top
        V_N029,     // V3B Grid
        V_N017,     // V3B Plate
        V_N035,     // V3B Cathode
        V_N010,     // Power Supply Node (VC)
        V_N018,     // Coupling between V3B and V4A
        V_N030,     // V4A Grid
        V_N025,     // V4A Plate
        V_N034,     // V4A Cathode
        V_N026,     // Coupling between V4A and Mixer
        V_N002      // Final Output / Mixer Node
        // V-sources at indices 13 and 14
    };

    // Component values
    double R10, R11, R21, R22, R23, R24, R25, R26, R27, R30, R31, R32;
    double RA_LEAD_DRIVE, RC_LEAD_DRIVE;
    double C10, C11, C21, C22, C23, C24, C25, C29, C30, C31, C32;

    std::vector<double> cap_z_state;
    const double VC = 410.0;

    void stampLinear() override {
        // --- Sources Topology ---
        stampVoltageSource_A(V_N001, GND, 0); // Input from V1B (v_idx=0)
        stampVoltageSource_A(V_N010, GND, 1); // VC Power Supply (v_idx=1)
        b_linear[NumNodes + 1] = VC;

        // --- Resistors ---
        stampResistorLinear(V_N002, V_N001, R10);  // Rhythm Path
        stampResistorLinear(V_N011, V_N024, R21);  // Lead Path
        stampResistorLinear(V_N024, V_N029, RA_LEAD_DRIVE);
        stampResistorLinear(V_N029, GND, RC_LEAD_DRIVE);
        stampResistorLinear(V_N029, GND, R22);
        stampResistorLinear(V_N017, V_N010, R26);
        stampResistorLinear(V_N035, GND, R23);
        stampResistorLinear(V_N030, V_N018, R25);
        stampResistorLinear(V_N030, GND, R24);
        stampResistorLinear(V_N025, V_N010, R27);
        stampResistorLinear(V_N034, GND, R30);
        stampResistorLinear(V_N002, V_N026, R31);
        stampResistorLinear(V_N002, GND, R11);
        stampResistorLinear(V_N002, GND, R32);

        // --- Capacitor Conductances ---
        stampCapacitor_A(V_N002, V_N001, C10);
        stampCapacitor_A(V_N011, V_N001, C21);
        stampCapacitor_A(V_N035, GND, C23);
        stampCapacitor_A(V_N035, V_N029, C22);
        stampCapacitor_A(V_N018, V_N017, C25);
        stampCapacitor_A(V_N030, GND, C24);
        stampCapacitor_A(V_N034, GND, C29);
        stampCapacitor_A(V_N026, V_N025, C30);
        stampCapacitor_A(V_N002, V_N026, C31);
        stampCapacitor_A(V_N002, GND, C11);
        stampCapacitor_A(V_N002, GND, C32);
    }

    void stampDynamic(double in) override {
        stampVoltageSource_b(0, in);
        stampCapacitor_b(V_N002, V_N001, cap_z_state[0]);
        stampCapacitor_b(V_N011, V_N001, cap_z_state[1]);
        stampCapacitor_b(V_N035, GND, cap_z_state[2]);
        stampCapacitor_b(V_N035, V_N029, cap_z_state[3]);
        stampCapacitor_b(V_N018, V_N017, cap_z_state[4]);
        stampCapacitor_b(V_N030, GND, cap_z_state[5]);
        stampCapacitor_b(V_N034, GND, cap_z_state[6]);
        stampCapacitor_b(V_N026, V_N025, cap_z_state[7]);
        stampCapacitor_b(V_N002, V_N026, cap_z_state[8]);
        stampCapacitor_b(V_N002, GND, cap_z_state[9]);
        stampCapacitor_b(V_N002, GND, cap_z_state[10]);
    }

    void stampNonLinear(const std::array<double, NumUnknowns>& x) override {
        // V3B
        Triode::State ts3b = Triode::calculate(x[V_N017]-x[V_N035], x[V_N029]-x[V_N035]);
        stampConductance(V_N017, V_N035, ts3b.g_p);
        A[V_N017][V_N029]+=ts3b.g_g; A[V_N017][V_N035]-=ts3b.g_g;
        A[V_N035][V_N029]-=ts3b.g_g; A[V_N035][V_N035]+=ts3b.g_g;
        stampConductance(V_N029, V_N035, ts3b.g_ig);
        // V4A
        Triode::State ts4a = Triode::calculate(x[V_N025]-x[V_N034], x[V_N030]-x[V_N034]);
        stampConductance(V_N025, V_N034, ts4a.g_p);
        A[V_N025][V_N030]+=ts4a.g_g; A[V_N025][V_N034]-=ts4a.g_g;
        A[V_N034][V_N030]-=ts4a.g_g; A[V_N034][V_N034]+=ts4a.g_g;
        stampConductance(V_N030, V_N034, ts4a.g_ig);
        stampNonLinear_b_only(x);
    }

    void stampNonLinear_b_only(const std::array<double, NumUnknowns>& x) override {
        // V3B
        Triode::State ts3b = Triode::calculate(x[V_N017]-x[V_N035], x[V_N029]-x[V_N035]);
        double i_p3b_lin=ts3b.ip-ts3b.g_p*(x[V_N017]-x[V_N035])-ts3b.g_g*(x[V_N029]-x[V_N035]);
        stampCurrentSource(V_N017, V_N035, i_p3b_lin);
        double i_g3b_lin=ts3b.ig-ts3b.g_ig*(x[V_N029]-x[V_N035]);
        stampCurrentSource(V_N029, V_N035, i_g3b_lin);
        // V4A
        Triode::State ts4a = Triode::calculate(x[V_N025]-x[V_N034], x[V_N030]-x[V_N034]);
        double i_p4a_lin=ts4a.ip-ts4a.g_p*(x[V_N025]-x[V_N034])-ts4a.g_g*(x[V_N030]-x[V_N034]);
        stampCurrentSource(V_N025, V_N034, i_p4a_lin);
        double i_g4a_lin=ts4a.ig-ts4a.g_ig*(x[V_N030]-x[V_N034]);
        stampCurrentSource(V_N030, V_N034, i_g4a_lin);
    }

public:
    LeadAndMixerStage() : cap_z_state(11, 0.0) {
        R10=3.3e6;  R11=680e3; R21=680e3; R22=475e3; R23=1.5e3;  R24=68e3;
        R25=270e3;  R26=82e3;  R27=274e3; R30=3.3e3; R31=220e3;  R32=100e3;
        C10=20e-12; C11=47e-12; C21=0.02e-6; C22=120e-12; C23=2.2e-6; C24=1000e-12;
        C25=0.022e-6; C29=0.22e-6; C30=0.047e-6; C31=250e-12; C32=500e-12;
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
        updateCapacitorState(x[V_N002], x[V_N001], C10, cap_z_state[0]);
        updateCapacitorState(x[V_N011], x[V_N001], C21, cap_z_state[1]);
        updateCapacitorState(x[V_N035], 0, C23, cap_z_state[2]);
        updateCapacitorState(x[V_N035], x[V_N029], C22, cap_z_state[3]);
        updateCapacitorState(x[V_N018], x[V_N017], C25, cap_z_state[4]);
        updateCapacitorState(x[V_N030], 0, C24, cap_z_state[5]);
        updateCapacitorState(x[V_N034], 0, C29, cap_z_state[6]);
        updateCapacitorState(x[V_N026], x[V_N025], C30, cap_z_state[7]);
        updateCapacitorState(x[V_N002], x[V_N026], C31, cap_z_state[8]);
        updateCapacitorState(x[V_N002], 0, C11, cap_z_state[9]);
        updateCapacitorState(x[V_N002], 0, C32, cap_z_state[10]);
        return x[V_N002];
    }
};
