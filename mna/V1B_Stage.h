#pragma once

#include "MNASolver.h"
#include "Triode.h"
#include <vector>

// This class models the second tube stage (V1B) of the preamp.
// It is a standard common-cathode gain stage.
//
// Input:  Voltage at node N020 (from the Volume 1 pot wiper).
// Output: Voltage at node N001.
//
// SPICE components included: XV1B, R8, R7, C13, C7, R9.
// The system has 4 nodes and 1 voltage source (the input signal).
class V1B_Stage : public MNASolver<4, 1> {
private:
    // A single enum for all unknowns in the solution vector 'x'.
    // The first 4 are node voltages (0-based indices).
    // The last 1 is the current through the input voltage source.
    enum Var {
        // Node Voltages (alphabetical by node name)
        V_N001 = 0, // Output node
        V_N009,     // Plate node
        V_N020,     // Grid / Input node
        V_N027,     // Cathode node
        // V-Source Currents
        I_VIN // This will be at index 4 (NumNodes)
    };

    // Component values
    double R7, R8, R9;
    double C7, C13;
    
    // History for capacitors
    std::vector<double> cap_hist;

    // Power supply voltage connected to the plate load resistor
    const double VE = 405.0;

    void stampComponents(double inputVoltage) {
        resetMatrices();
        
        // --- Input Source ---
        // The input signal is the voltage at node N020.
        stampVoltageSource(V_N020 + 1, 0, 0, inputVoltage);

        // --- Static Components ---
        // Plate load resistor (R8) connected to power supply VE
        stampResistor(V_N009 + 1, 0, R8);
        b[V_N009] += VE / R8; // Model VE as a current source contribution

        // Cathode resistor (R7)
        stampResistor(V_N027 + 1, 0, R7);
        
        // Load resistor / next stage's grid-leak (R9)
        stampResistor(V_N001 + 1, 0, R9);

        // --- Capacitors ---
        // Cathode bypass capacitor (C13)
        stampCapacitor(V_N027 + 1, 0, C13, cap_hist[0]);
        
        // Output coupling capacitor (C7)
        stampCapacitor(V_N001 + 1, V_N009 + 1, C7, cap_hist[1]);
    }
    
    void stampNonLinear(const std::array<double, NumUnknowns>& current_x) {
        // XV1B: Plate=N009, Grid=N020, Cathode=N027
        double v_p = current_x[V_N009];
        double v_g = current_x[V_N020];
        double v_c = current_x[V_N027];

        Triode::State ts = Triode::calculate(v_p - v_c, v_g - v_c);
        
        // --- Plate Current ---
        double i_p_lin = ts.ip - ts.g_p * (v_p - v_c) - ts.g_g * (v_g - v_c);
        stampCurrentSource(V_N009 + 1, V_N027 + 1, i_p_lin);
        stampConductance(V_N009 + 1, V_N027 + 1, ts.g_p);
        
        // Transconductance (Vg influence on Ip)
        A[V_N009][V_N020] += ts.g_g;
        A[V_N009][V_N027] -= ts.g_g;
        A[V_N027][V_N020] -= ts.g_g;
        A[V_N027][V_N027] += ts.g_g;

        // --- Grid Current ---
        double i_g_lin = ts.ig - ts.g_ig * (v_g - v_c);
        stampCurrentSource(V_N020 + 1, V_N027 + 1, i_g_lin);
        stampConductance(V_N020 + 1, V_N027 + 1, ts.g_ig);
    }

public:
    V1B_Stage() : cap_hist(2, 0.0) {
        // Set fixed component values from SPICE netlist
        R7 = 1.5e3;
        R8 = 100e3;
        R9 = 100e3;
        C7 = 0.1e-6;
        C13 = 22e-6;
    }

    // No user-adjustable parameters in this stage
    
    double process(double in) {
        const int MAX_ITER = 15;
        const double CONVERGENCE_THRESH = 1e-6;
        
        std::array<double, NumUnknowns> current_x = x; // Start iteration with last solution

        for (int i = 0; i < MAX_ITER; ++i) {
            stampComponents(in);
            stampNonLinear(current_x);
            
            if (!lu_decompose()) break; // Matrix is singular, abort
            
            std::array<double, NumUnknowns> delta_x;
            lu_solve(delta_x);
            
            double norm = 0.0;
            for(size_t j = 0; j < NumUnknowns; ++j) {
                current_x[j] += delta_x[j];
                norm += delta_x[j] * delta_x[j];
            }
            if (sqrt(norm) < CONVERGENCE_THRESH) break; // Converged
        }
        x = current_x; // Store solution for next sample and history update
        
        // Update capacitor history using the new voltages
        updateCapacitorState(x[V_N027], 0,         C13, cap_hist[0]); // C13
        updateCapacitorState(x[V_N001], x[V_N009], C7,  cap_hist[1]); // C7

        // The output of this stage is the voltage at node N001
        return x[V_N001];
    }
};
