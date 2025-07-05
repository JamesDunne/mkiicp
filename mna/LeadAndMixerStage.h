#pragma once

#include "MNASolver.h"
#include "Triode.h"
#include <vector>

// This is the definitive "super-stage" that models the V1B driver, the entire lead channel,
// and the mixer network together.
//
// THIS IS THE FINAL CORRECTED VERSION: It meticulously implements the V1B output coupling,
// which fixes the critical DC offset and gain-staging bugs. This should finally deliver
// the correct high-gain tone and response.
class LeadAndMixerStage : public MNASolver<17, 3> {
private:
    static constexpr int GND = -1;
    enum Var {
        V_N020 = 0, V_N009, V_N027, V_N003, V_N001, V_N011, V_N024, V_N029,
        V_N017, V_N035, V_N010, V_N018, V_N030, V_N025, V_N034, V_N026, V_N002
    };

    double R7,R8,R9,R10,R11,R21,R22,R23,R24,R25,R26,R27,R30,R31,R32;
    double RA_LEAD_DRIVE, RC_LEAD_DRIVE;
    double C7,C10,C11,C13,C21,C22,C23,C24,C25,C29,C30,C31,C32;

    std::vector<double> cap_z_state;
    const double VE = 405.0;
    const double VC = 410.0;

    void stampLinear() override {
        // --- Sources Topology ---
        stampVoltageSource_A(V_N020, GND, 0);
        stampVoltageSource_A(V_N003, GND, 1);
        stampVoltageSource_A(V_N010, GND, 2);
        b_linear[NumNodes + 1] = VE;
        b_linear[NumNodes + 2] = VC;

        // --- V1B Stage ---
        stampResistorLinear(V_N027, GND, R7);
        stampResistorLinear(V_N009, V_N003, R8);
        stampCapacitor_A(V_N027, GND, C13);

        // --- V1B Output Coupling (CRITICAL, NOW CORRECT) ---
        // C7 AC-couples the plate N009 to the junction node N001. This blocks DC.
        stampCapacitor_A(V_N001, V_N009, C7);
        // R9 provides the DC path to ground for N001 and everything connected to it.
        // This forces the DC voltage at N001 (and thus the lead grid) to be 0V.
        stampResistorLinear(V_N001, GND, R9);

        // --- Mixer (Rhythm Path from N001) ---
        stampResistorLinear(V_N002, V_N001, R10);
        stampCapacitor_A(V_N002, V_N001, C10);

        // --- Lead Path (Input from N001) ---
        stampCapacitor_A(V_N011, V_N001, C21);
        stampResistorLinear(V_N011, V_N024, R21);
        stampResistorLinear(V_N024, V_N029, RA_LEAD_DRIVE);
        stampResistorLinear(V_N029, GND, RC_LEAD_DRIVE);
        stampResistorLinear(V_N029, GND, R22);

        // --- V3B Stage ---
        stampResistorLinear(V_N017, V_N010, R26);
        stampResistorLinear(V_N035, GND, R23);
        stampCapacitor_A(V_N035, GND, C23);
        stampCapacitor_A(V_N035, V_N029, C22);

        // --- V3B->V4A Coupling ---
        stampCapacitor_A(V_N018, V_N017, C25);
        stampResistorLinear(V_N030, V_N018, R25);
        stampResistorLinear(V_N030, GND, R24);
        stampCapacitor_A(V_N030, GND, C24);

        // --- V4A Stage ---
        stampResistorLinear(V_N025, V_N010, R27);
        stampResistorLinear(V_N034, GND, R30);
        stampCapacitor_A(V_N034, GND, C29);

        // --- V4A->Mixer Coupling ---
        stampCapacitor_A(V_N026, V_N025, C30);
        stampResistorLinear(V_N002, V_N026, R31);
        stampCapacitor_A(V_N002, V_N026, C31);

        // --- Final Mixer Load on N002 ---
        stampResistorLinear(V_N002, GND, R11);
        stampResistorLinear(V_N002, GND, R32);
        stampCapacitor_A(V_N002, GND, C11);
        stampCapacitor_A(V_N002, GND, C32);
    }

    void stampDynamic(double in) override {
        stampVoltageSource_b(0, in);
        stampCapacitor_b(V_N027, GND, cap_z_state[0]);
        stampCapacitor_b(V_N001, V_N009, cap_z_state[1]);
        stampCapacitor_b(V_N002, V_N001, cap_z_state[2]);
        stampCapacitor_b(V_N011, V_N001, cap_z_state[3]);
        stampCapacitor_b(V_N035, GND, cap_z_state[4]);
        stampCapacitor_b(V_N035, V_N029, cap_z_state[5]);
        stampCapacitor_b(V_N018, V_N017, cap_z_state[6]);
        stampCapacitor_b(V_N030, GND, cap_z_state[7]);
        stampCapacitor_b(V_N034, GND, cap_z_state[8]);
        stampCapacitor_b(V_N026, V_N025, cap_z_state[9]);
        stampCapacitor_b(V_N002, V_N026, cap_z_state[10]);
        stampCapacitor_b(V_N002, GND, cap_z_state[11]);
        stampCapacitor_b(V_N002, GND, cap_z_state[12]);
    }

    void stampDynamic_b_only(double in, std::array<double, NumUnknowns>& b_vector) override {
        stampVoltageSource_b_vector(0, in, b_vector);
        stampCapacitor_b_vector(V_N027, GND, cap_z_state[0], b_vector);
        stampCapacitor_b_vector(V_N001, V_N009, cap_z_state[1], b_vector);
        stampCapacitor_b_vector(V_N002, V_N001, cap_z_state[2], b_vector);
        stampCapacitor_b_vector(V_N011, V_N001, cap_z_state[3], b_vector);
        stampCapacitor_b_vector(V_N035, GND, cap_z_state[4], b_vector);
        stampCapacitor_b_vector(V_N035, V_N029, cap_z_state[5], b_vector);
        stampCapacitor_b_vector(V_N018, V_N017, cap_z_state[6], b_vector);
        stampCapacitor_b_vector(V_N030, GND, cap_z_state[7], b_vector);
        stampCapacitor_b_vector(V_N034, GND, cap_z_state[8], b_vector);
        stampCapacitor_b_vector(V_N026, V_N025, cap_z_state[9], b_vector);
        stampCapacitor_b_vector(V_N002, V_N026, cap_z_state[10], b_vector);
        stampCapacitor_b_vector(V_N002, GND, cap_z_state[11], b_vector);
        stampCapacitor_b_vector(V_N002, GND, cap_z_state[12], b_vector);
    }

    void subtractSources(std::array<double, 20> &residual, double in) const override {
        // v_idx = 0: Input source
        residual[NumNodes + 0] -= in;
        // v_idx = 1: VE Power Supply
        residual[NumNodes + 1] -= VE;
        // v_idx = 2: VC Power Supply
        residual[NumNodes + 2] -= VC;
    }

private:
    void stampNonLinear(const std::array<double, NumUnknowns>& x) override {
        // Calculate the state of each triode ONCE.
        Triode::State ts1b = Triode::calculate(x[V_N009]-x[V_N027], x[V_N020]-x[V_N027]);
        Triode::State ts3b = Triode::calculate(x[V_N017]-x[V_N035], x[V_N029]-x[V_N035]);
        Triode::State ts4a = Triode::calculate(x[V_N025]-x[V_N034], x[V_N030]-x[V_N034]);

        // --- Stamp V1B ---
        // A Matrix (Conductances)
        stampConductance(V_N009, V_N027, ts1b.g_p);
        A[V_N009][V_N020]+=ts1b.g_g; A[V_N009][V_N027]-=ts1b.g_g;
        A[V_N027][V_N020]-=ts1b.g_g; A[V_N027][V_N027]+=ts1b.g_g;
        stampConductance(V_N020, V_N027, ts1b.g_ig);
        // b Vector (Linearized Current Sources)
        double i_p1b_lin=ts1b.ip-ts1b.g_p*(x[V_N009]-x[V_N027])-ts1b.g_g*(x[V_N020]-x[V_N027]);
        stampCurrentSource(V_N009, V_N027, i_p1b_lin);
        double i_g1b_lin=ts1b.ig-ts1b.g_ig*(x[V_N020]-x[V_N027]);
        stampCurrentSource(V_N020, V_N027, i_g1b_lin);

        // --- Stamp V3B ---
        // A Matrix
        stampConductance(V_N017, V_N035, ts3b.g_p);
        A[V_N017][V_N029]+=ts3b.g_g; A[V_N017][V_N035]-=ts3b.g_g;
        A[V_N035][V_N029]-=ts3b.g_g; A[V_N035][V_N035]+=ts3b.g_g;
        stampConductance(V_N029, V_N035, ts3b.g_ig);
        // b Vector
        double i_p3b_lin=ts3b.ip-ts3b.g_p*(x[V_N017]-x[V_N035])-ts3b.g_g*(x[V_N029]-x[V_N035]);
        stampCurrentSource(V_N017, V_N035, i_p3b_lin);
        double i_g3b_lin=ts3b.ig-ts3b.g_ig*(x[V_N029]-x[V_N035]);
        stampCurrentSource(V_N029, V_N035, i_g3b_lin);

        // --- Stamp V4A ---
        // A Matrix
        stampConductance(V_N025, V_N034, ts4a.g_p);
        A[V_N025][V_N030]+=ts4a.g_g; A[V_N025][V_N034]-=ts4a.g_g;
        A[V_N034][V_N030]-=ts4a.g_g; A[V_N034][V_N034]+=ts4a.g_g;
        stampConductance(V_N030, V_N034, ts4a.g_ig);
        // b Vector
        double i_p4a_lin=ts4a.ip-ts4a.g_p*(x[V_N025]-x[V_N034])-ts4a.g_g*(x[V_N030]-x[V_N034]);
        stampCurrentSource(V_N025, V_N034, i_p4a_lin);
        double i_g4a_lin=ts4a.ig-ts4a.g_ig*(x[V_N030]-x[V_N034]);
        stampCurrentSource(V_N030, V_N034, i_g4a_lin);
    }

        // Stamps ONLY the Jacobian (conductances)
    void stampNonLinear_A(const std::array<double, NumUnknowns>& x) override {
        // V1B
        Triode::State ts1b = Triode::calculate(x[V_N009]-x[V_N027], x[V_N020]-x[V_N027]);
        stampConductance(V_N009, V_N027, ts1b.g_p);
        A[V_N009][V_N020]+=ts1b.g_g; A[V_N009][V_N027]-=ts1b.g_g;
        A[V_N027][V_N020]-=ts1b.g_g; A[V_N027][V_N027]+=ts1b.g_g;
        stampConductance(V_N020, V_N027, ts1b.g_ig);
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
    }

    // Adds ONLY the pure non-linear currents
    void addNonlinear_b(std::array<double, NumUnknowns>& b, const std::array<double, NumUnknowns>& x) const override {
        // V1B
        Triode::State ts1b = Triode::calculate(x[V_N009]-x[V_N027], x[V_N020]-x[V_N027]);
        b[V_N009] += ts1b.ip; b[V_N027] -= ts1b.ip;
        b[V_N020] += ts1b.ig; b[V_N027] -= ts1b.ig;
        // V3B
        Triode::State ts3b = Triode::calculate(x[V_N017]-x[V_N035], x[V_N029]-x[V_N035]);
        b[V_N017] += ts3b.ip; b[V_N035] -= ts3b.ip;
        b[V_N029] += ts3b.ig; b[V_N035] -= ts3b.ig;
        // V4A
        Triode::State ts4a = Triode::calculate(x[V_N025]-x[V_N034], x[V_N030]-x[V_N034]);
        b[V_N025] += ts4a.ip; b[V_N034] -= ts4a.ip;
        b[V_N030] += ts4a.ig; b[V_N034] -= ts4a.ig;
    }

public:
    LeadAndMixerStage() : cap_z_state(13, 0.0) {
        R7=1.5e3; R8=100e3; R9=100e3; R10=3.3e6; R11=680e3; R21=680e3; R22=475e3;
        R23=1.5e3; R24=68e3; R25=270e3; R26=82e3; R27=274e3; R30=3.3e3; R31=220e3; R32=100e3;
        C7=0.1e-6; C10=20e-12; C11=47e-12; C13=22e-6; C21=0.02e-6; C22=120e-12; C23=2.2e-6;
        C24=1000e-12; C25=0.022e-6; C29=0.22e-6; C30=0.047e-6; C31=250e-12; C32=500e-12;
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
        // solveNonlinear(in);
        solveNonlinear_Adaptive(in);
        // solveNonlinear_Real(in);
        updateCapacitorState(x[V_N027], 0, C13, cap_z_state[0]);
        updateCapacitorState(x[V_N001], x[V_N009], C7, cap_z_state[1]);
        updateCapacitorState(x[V_N002], x[V_N001], C10, cap_z_state[2]);
        updateCapacitorState(x[V_N011], x[V_N001], C21, cap_z_state[3]);
        updateCapacitorState(x[V_N035], 0, C23, cap_z_state[4]);
        updateCapacitorState(x[V_N035], x[V_N029], C22, cap_z_state[5]);
        updateCapacitorState(x[V_N018], x[V_N017], C25, cap_z_state[6]);
        updateCapacitorState(x[V_N030], 0, C24, cap_z_state[7]);
        updateCapacitorState(x[V_N034], 0, C29, cap_z_state[8]);
        updateCapacitorState(x[V_N026], x[V_N025], C30, cap_z_state[9]);
        updateCapacitorState(x[V_N002], x[V_N026], C31, cap_z_state[10]);
        updateCapacitorState(x[V_N002], GND, C11, cap_z_state[11]);
        updateCapacitorState(x[V_N002], GND, C32, cap_z_state[12]);
        return x[V_N002];
    }
};
