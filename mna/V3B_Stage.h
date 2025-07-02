#pragma once

#include "MNASolver.h"
#include "Triode.h"
#include <vector>

// This class models the third tube stage (V3B), which is the first gain
// stage of the high-gain lead channel.
//
// Input:  Voltage at node N029 (from the Lead Drive pot).
// Output: Voltage at node N017 (the plate of the triode).
//
// SPICE components included: XV3B, R26, R23, C23, C22.
// The system has 3 unknown nodes and 1 voltage source for the input.
class V3B_Stage : public MNASolver<3, 1> {
private:
    static constexpr int GND = -1;

    // The unknowns are the voltages at the triode's terminals.
    enum Var {
        V_N017 = 0, // Plate node (output)
        V_N029,     // Grid node (input)
        V_N035      // Cathode node
        // I_VIN at index 3 (NumNodes)
    };

    // Component values
    double R23, R26;
    double C22, C23;
    
    // State variables for the capacitors
    std::vector<double> cap_z_state;

    // Power supply voltage for this stage
    const double VC = 410.0;

    void stampComponents(double inputVoltage) {
        resetMatrices();
        
        // --- Input Source ---
        // The input signal is the voltage at node N029.
        stampVoltageSource(V_N029, GND, 0, inputVoltage);

        // --- Static Components ---
        // Plate load resistor (R26) connected to power supply VC
        stampResistor(V_N017, GND, R26);
        b[V_N017] += VC / R26;

        // Cathode resistor (R23)
        stampResistor(V_N035, GND, R23);

        // --- Capacitors ---
        // Cathode bypass capacitor (C23)
        stampCapacitor(V_N035, GND, C23, cap_z_state[0]);
        
        // Grid-to-cathode high-cut capacitor (C22)
        stampCapacitor(V_N035, V_N029, C22, cap_z_state[1]);
    }
    
    void stampNonLinear(const std::array<double, NumUnknowns>& current_x) {
        // XV3B: Plate=N017, Grid=N029, Cathode=N035
        double v_p = current_x[V_N017];
        double v_g = current_x[V_N029];
        double v_c = current_x[V_N035];

        Triode::State ts = Triode::calculate(v_p - v_c, v_g - v_c);
        
        // --- Plate Current ---
        double i_p_lin = ts.ip - ts.g_p * (v_p - v_c) - ts.g_g * (v_g - v_c);
        stampCurrentSource(V_N017, V_N035, i_p_lin);
        stampConductance(V_N017, V_N035, ts.g_p);
        
        A[V_N017][V_N029] += ts.g_g;
        A[V_N017][V_N035] -= ts.g_g;
        A[V_N035][V_N029] -= ts.g_g;
        A[V_N035][V_N035] += ts.g_g;

        // --- Grid Current ---
        double i_g_lin = ts.ig - ts.g_ig * (v_g - v_c);
        stampCurrentSource(V_N029, V_N035, i_g_lin);
        stampConductance(V_N029, V_N035, ts.g_ig);
    }

public:
    V3B_Stage() : cap_z_state(2, 0.0) {
        // Set fixed component values from SPICE netlist
        R23 = 1.5e3;
        R26 = 82e3;
        C22 = 120e-12;
        C23 = 2.2e-6;
    }

    // No user-adjustable parameters in this stage.
    
    double process(double in) {
        const int MAX_ITER = 15;
        const double CONVERGENCE_THRESH = 1e-6;
        
        std::array<double, NumUnknowns> current_x = x;

        for (int i = 0; i < MAX_ITER; ++i) {
            stampComponents(in);
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
        
        // Update capacitor states for the next time step
        updateCapacitorState(x[V_N035], 0, C23, cap_z_state[0]);
        updateCapacitorState(x[V_N035], x[V_N029], C22, cap_z_state[1]);

        // The output of this stage is the voltage at the plate (node N017)
        return x[V_N017];
    }
};
