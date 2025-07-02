#pragma once

#include "MNASolver.h"
#include "Triode.h"
#include <vector>

// This class models the main mixing stage and the V2B gain stage.
// It correctly combines the Rhythm (N001) and Lead (N026) signal paths
// into the mixer node N002, then amplifies the result with triode V2B.
// It also includes the full AC load on the V2B plate to ensure accuracy.
//
// Inputs: `in_rhythm` (from N001), `in_lead` (from N026).
// Output: Voltage at node N022.
//
// The system has 6 unknown nodes and 1 voltage source (power supply).
class MixerAndV2B : public MNASolver<6, 1> {
private:
    static constexpr int GND = -1;

    enum Var {
        V_N002 = 0, // Mixer Node / V2B Grid
        V_N036,     // V2B Cathode
        V_N021,     // V2B Plate
        V_N006,     // Power Supply Node
        V_P001,     // Output coupling intermediate node
        V_N022      // Output node
        // I_VSUPPLY at index 6
    };

    // Component values
    double R10, R11, R13, R16, R31, R32, R46, R105;
    double C9, C10, C11, C31, C32;
    
    std::vector<double> cap_z_state;

    const double VC2 = 410.0;

    void stampComponents(double in_rhythm, double in_lead) {
        resetMatrices();
        
        // --- Power Supply ---
        stampVoltageSource(V_N006, GND, 0, VC2);

        // --- Mixer Network ---
        // Stamping of input components requires special handling as the inputs
        // are known voltages, not unknown nodes in this MNA system.
        
        // Rhythm path mixing (R10, C10 from in_rhythm to N002)
        double G_R10 = 1.0 / R10;
        A[V_N002][V_N002] += G_R10;
        b[V_N002] += G_R10 * in_rhythm;
        stampCapacitor(V_N002, -1, C10, cap_z_state[0]); // Pass GND and handle input below
        b[V_N002] += (invT_2 * C10) * in_rhythm;

        // Lead path mixing (R31, C31 from in_lead to N002)
        double G_R31 = 1.0 / R31;
        A[V_N002][V_N002] += G_R31;
        b[V_N002] += G_R31 * in_lead;
        stampCapacitor(V_N002, -1, C31, cap_z_state[1]); // Pass GND and handle input below
        b[V_N002] += (invT_2 * C31) * in_lead;

        // Mixer load to ground
        stampResistor(V_N002, GND, R11);
        stampCapacitor(V_N002, GND, C11, cap_z_state[2]);
        stampResistor(V_N002, GND, R32);
        stampCapacitor(V_N002, GND, C32, cap_z_state[3]);

        // --- V2B Stage ---
        stampResistor(V_N021, V_N006, R13); // Plate load
        stampResistor(V_N036, GND, R16);    // Cathode resistor
        
        // --- V2B Output Load ---
        stampCapacitor(V_P001, V_N021, C9, cap_z_state[4]);
        stampResistor(V_N022, V_P001, R105);
        stampResistor(V_N022, GND, R46);
    }
    
    void stampNonLinear(const std::array<double, NumUnknowns>& current_x) {
        // XV2B: Plate=N021, Grid=N002, Cathode=N036
        double v_p = current_x[V_N021];
        double v_g = current_x[V_N002];
        double v_c = current_x[V_N036];
        Triode::State ts = Triode::calculate(v_p - v_c, v_g - v_c);
        
        double i_p_lin = ts.ip - ts.g_p * (v_p - v_c) - ts.g_g * (v_g - v_c);
        stampCurrentSource(V_N021, V_N036, i_p_lin);
        stampConductance(V_N021, V_N036, ts.g_p);
        
        A[V_N021][V_N002] += ts.g_g;
        A[V_N021][V_N036] -= ts.g_g;
        A[V_N036][V_N002] -= ts.g_g;
        A[V_N036][V_N036] += ts.g_g;

        double i_g_lin = ts.ig - ts.g_ig * (v_g - v_c);
        stampCurrentSource(V_N002, V_N036, i_g_lin);
        stampConductance(V_N002, V_N036, ts.g_ig);
    }

public:
    MixerAndV2B() : cap_z_state(5, 0.0) {
        R10 = 3.3e6;   C10 = 20e-12;
        R11 = 680e3;   C11 = 47e-12;
        R31 = 220e3;   C31 = 250e-12;
        R32 = 100e3;   C32 = 500e-12;
        R13 = 100e3;
        R16 = 1.5e3;
        C9 = 0.047e-6;
        R105 = 47e3;
        R46 = 47e3;
    }
    
    double process(double in_rhythm, double in_lead) {
        const int MAX_ITER = 15;
        const double CONVERGENCE_THRESH = 1e-6;
        
        std::array<double, NumUnknowns> current_x = x;

        for (int i = 0; i < MAX_ITER; ++i) {
            stampComponents(in_rhythm, in_lead);
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
        
        updateCapacitorState(x[V_N002], in_rhythm, C10, cap_z_state[0]);
        updateCapacitorState(x[V_N002], in_lead, C31, cap_z_state[1]);
        updateCapacitorState(x[V_N002], 0, C11, cap_z_state[2]);
        updateCapacitorState(x[V_N002], 0, C32, cap_z_state[3]);
        updateCapacitorState(x[V_P001], x[V_N021], C9, cap_z_state[4]);

        return x[V_N022];
    }
};
