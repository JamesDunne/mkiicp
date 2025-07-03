#pragma once

#include "MNASolver.h"
#include "Triode.h"
#include <vector>
#include <algorithm> // for std::max
#include <cmath>     // for std::abs, std::sqrt

// This class models the V2B gain stage and its surrounding network.
//
// THIS IS THE DEFINITIVE VERSION: It includes the full complex load on the V2B plate,
// ensuring the DC operating point and AC swing are physically correct.
//
// IT NOW INCLUDES A DAMPENED NEWTON-RAPHSON SOLVER to prevent divergence
// during large signal transients, fixing the "digital clipping" issue.
//
// Inputs: `in_lead`, `in_rhythm`.
// Output: Voltage at node N023 (the grid of the final tube stage, V2A).
//
// The system has 8 unknown nodes and 2 voltage sources.
class MixerAndV2B : public MNASolver<8, 2> {
private:
    class RhythmPathProcessor : public MNASolver<1, 0> {
    private:
        static constexpr int GND = -1;
        enum Var { V_OUT = 0 };
        double R10, C10;
        double cap_z;
    public:
        RhythmPathProcessor() : cap_z(0.0) { R10 = 3.3e6; C10 = 20e-12; }
        void setup(double sr) { MNASolver<1,0>::setup(sr); }
        double process(double in) {
            resetMatrices();
            A[V_OUT][V_OUT] += 1.0 / R10;
            b[V_OUT] += in / R10;
            stampCapacitor(V_OUT, -1, C10, cap_z);
            b[V_OUT] += (invT_2 * C10) * in;
            if (A[V_OUT][V_OUT] > 1e-12) { x[V_OUT] = b[V_OUT] / A[V_OUT][V_OUT]; } else { x[V_OUT] = 0.0; }
            updateCapacitorState(x[V_OUT], in, C10, cap_z);
            return x[V_OUT];
        }
    };
    RhythmPathProcessor rhythmProcessor;

    static constexpr int GND = -1;
    enum Var {
        V_N002 = 0, // V2B Grid
        V_N036,     // V2B Cathode
        V_N021,     // V2B Plate
        V_P001,     // Output coupling intermediate node
        V_N022,     // Intermediate node in load network
        V_N006,     // Power Supply Node
        V_N032,     // New node in the final load
        V_N023      // New final output node (V2A Grid)
        // V-source currents I_Vin and I_Vsupply are at indices 8 and 9
    };

    double R11, R12, R13, R16, R46, R101, R102, R103, R105;
    double C9, C11;
    std::vector<double> cap_z_state;
    const double VC2 = 410.0;

    void stampLinear() override {
        stampVoltageSource_A(V_N002, GND, 0); // in
        stampVoltageSource_A(V_N006, GND, 1); // VC2

        // Stamp all resistors (they are constant)
        stampResistorLinear(V_N021, V_N006, R13);
        stampResistorLinear(V_N036, GND, R16);
        stampResistorLinear(V_N002, GND, R11);
        stampResistorLinear(V_N022, V_P001, R105);
        stampResistorLinear(V_N022, GND, R46);
        stampResistorLinear(V_N022, V_N032, R102);
        stampResistorLinear(V_N032, GND, R101);
        stampResistorLinear(V_N023, V_N032, R12);
        stampResistorLinear(V_N023, GND, R103);

        // Stamp the constant, linear part of the capacitors
        stampCapacitor_A(V_N002, GND, C11);
        stampCapacitor_A(V_P001, V_N021, C9);
    }

    // 2. Stamps components that change each sample.
    void stampDynamic(double in) override {
        // Stamp the input signal source
        stampVoltageSource_b(0, in);
        stampVoltageSource_b(1, VC2);

        // Stamp only the dynamic history current of the capacitors
        stampCapacitor_b(V_N002, GND, cap_z_state[0]);    // C11
        stampCapacitor_b(V_P001, V_N021, cap_z_state[1]); // C9
    }

    void stampNonLinear(const std::array<double, NumUnknowns>& current_x) override {
        double v_p = current_x[V_N021];
        double v_g = current_x[V_N002];
        double v_c = current_x[V_N036];
        Triode::State ts = Triode::calculate(v_p - v_c, v_g - v_c);
        double i_p_lin = ts.ip - ts.g_p*(v_p-v_c) - ts.g_g*(v_g-v_c);
        stampCurrentSource(V_N021, V_N036, i_p_lin);
        stampConductance(V_N021, V_N036, ts.g_p);
        A[V_N021][V_N002] += ts.g_g; A[V_N021][V_N036] -= ts.g_g;
        A[V_N036][V_N002] -= ts.g_g; A[V_N036][V_N036] += ts.g_g;
        double i_g_lin = ts.ig - ts.g_ig*(v_g-v_c);
        stampCurrentSource(V_N002, V_N036, i_g_lin);
        stampConductance(V_N002, V_N036, ts.g_ig);
    }

public:
    MixerAndV2B() : cap_z_state(2, 0.0) {
        // V2B and mixer components
        R11 = 680e3; C11 = 47e-12; R13 = 100e3; R16 = 1.5e3;
        // Output coupling components
        C9 = 0.047e-6; R105 = 47e3; R46 = 47e3;
        // Newly added load components
        R102 = 150e3; R101 = 4.7e3; R12 = 2.2e3; R103 = 47e3;
    }

    void setup(double sr) {
        MNASolver<8, 2>::setup(sr);
        rhythmProcessor.setup(sr);
    }

    double process(double in_lead, double in_rhythm) {
        double rhythm_processed = rhythmProcessor.process(in_rhythm);
        double input_sum = in_lead + rhythm_processed;

        solveNonlinear(input_sum);

        updateCapacitorState(x[V_N002], 0, C11, cap_z_state[0]);
        updateCapacitorState(x[V_P001], x[V_N021], C9, cap_z_state[1]);

        return x[V_N023];
    }
};
