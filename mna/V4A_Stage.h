#pragma once

#include "MNASolver.h"
#include "Triode.h"
#include <vector>

// This class models the fourth tube stage (V4A), which is the second and final
// gain stage of the high-gain lead channel.
//
// Input:  Voltage at node N030 (from the first lead coupling network).
// Output: Voltage at node N025 (the plate of the triode).
//
// SPICE components included: XV4A, R27, R30, C29.
// The system has 3 unknown nodes and 1 voltage source for the input.
class V4A_Stage : public MNASolver<3, 1> {
private:
    static constexpr int GND = -1;

    // The unknowns are the voltages at the triode's terminals.
    enum Var {
        V_N025 = 0, // Plate node (output)
        V_N030,     // Grid node (input)
        V_N034      // Cathode node
        // I_VIN at index 3 (NumNodes)
    };

    // Component values
    double R27, R30;
    double C29;
    
    // State variable for the capacitor
    std::vector<double> cap_z_state;

    // Power supply voltage for this stage
    const double VC = 410.0;

    void stampComponents(double inputVoltage) {
        resetMatrices();
        
        // --- Input Source ---
        // The input signal is the voltage at node N030.
        stampVoltageSource(V_N030, GND, 0, inputVoltage);

        // --- Static Components ---
        // Plate load resistor (R27) connected to power supply VC
        stampResistor(V_N025, GND, R27);
        b[V_N025] += VC / R27;

        // Cathode resistor (R30)
        stampResistor(V_N034, GND, R30);

        // --- Capacitor ---
        // Cathode bypass capacitor (C29)
        stampCapacitor(V_N034, GND, C29, cap_z_state[0]);
    }
    
    void stampNonLinear(const std::array<double, NumUnknowns>& current_x) {
        // XV4A: Plate=N025, Grid=N030, Cathode=N034
        double v_p = current_x[V_N025];
        double v_g = current_x[V_N030];
        double v_c = current_x[V_N034];

        Triode::State ts = Triode::calculate(v_p - v_c, v_g - v_c);
        
        // --- Plate Current ---
        double i_p_lin = ts.ip - ts.g_p * (v_p - v_c) - ts.g_g * (v_g - v_c);
        stampCurrentSource(V_N025, V_N034, i_p_lin);
        stampConductance(V_N025, V_N034, ts.g_p);
        
        A[V_N025][V_N030] += ts.g_g;
        A[V_N025][V_N034] -= ts.g_g;
        A[V_N034][V_N030] -= ts.g_g;
        A[V_N034][V_N034] += ts.g_g;

        // --- Grid Current ---
        double i_g_lin = ts.ig - ts.g_ig * (v_g - v_c);
        stampCurrentSource(V_N030, V_N034, i_g_lin);
        stampConductance(V_N030, V_N034, ts.g_ig);
    }

public:
    V4A_Stage() : cap_z_state(1, 0.0) {
        // Set fixed component values from SPICE netlist
        R27 = 274e3;
        R30 = 3.3e3;
        C29 = 0.22e-6;
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
        
        // Update capacitor state for the next time step
        updateCapacitorState(x[V_N034], 0, C29, cap_z_state[0]);

        // The output of this stage is the voltage at the plate (node N025)
        return x[V_N025];
    }
};
