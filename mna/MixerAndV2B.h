#pragma once

#include "MNASolver.h"
#include "Triode.h"
#include <vector>

// This class models the V2B gain stage and its entire surrounding network.
//
// THIS IS THE DEFINITIVE VERSION: It includes the full complex load on the V2B plate,
// ensuring the DC operating point and AC swing are physically correct.
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
    // CORRECTED: Added nodes for the full output load network
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

    void stampComponents(double input_sum) {
        resetMatrices();

        stampVoltageSource(V_N002, GND, 0, input_sum);
        stampVoltageSource(V_N006, GND, 1, VC2);

        // --- Mixer & V2B Components ---
        stampResistor(V_N021, V_N006, R13); // Plate load
        stampResistor(V_N036, GND, R16);    // Cathode resistor
        stampResistor(V_N002, GND, R11);    // Grid-leak resistor
        stampCapacitor(V_N002, GND, C11, cap_z_state[0]);

        // --- Full Output Load Network ---
        stampCapacitor(V_P001, V_N021, C9, cap_z_state[1]); // Output coupling
        stampResistor(V_N022, V_P001, R105);
        stampResistor(V_N022, GND, R46);

        // CORRECTED: Add the newly included load resistors
        stampResistor(V_N022, V_N032, R102);
        stampResistor(V_N032, GND, R101);
        stampResistor(V_N023, V_N032, R12);
        stampResistor(V_N023, GND, R103);
    }

    void stampNonLinear(const std::array<double, NumUnknowns>& current_x) {
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

        const int MAX_ITER = 15;
        const double CONVERGENCE_THRESH = 1e-6;
        std::array<double, NumUnknowns> current_x = x;
        for (int i = 0; i < MAX_ITER; ++i) {
            stampComponents(input_sum);
            stampNonLinear(current_x);
            if (!lu_decompose()) { return 0.0; }
            std::array<double, NumUnknowns> next_x;
            lu_solve(next_x);
            double norm = 0.0;
            for(size_t j=0; j<NumUnknowns; ++j) { norm += (next_x[j] - current_x[j])*(next_x[j] - current_x[j]); }
            current_x = next_x;
            if (sqrt(norm) < CONVERGENCE_THRESH) { break; }
        }
        x = current_x;

        updateCapacitorState(x[V_N002], 0, C11, cap_z_state[0]);
        updateCapacitorState(x[V_P001], x[V_N021], C9, cap_z_state[1]);

        // CORRECTED: Return the voltage at the final node N023
        return x[V_N023];
    }
};
