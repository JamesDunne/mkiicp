#pragma once

#include "MNASolver.h"
#include <vector>

// This class models the passive coupling network between the V3B and V4A stages
// in the lead channel. It consists of a CR high-pass filter (C25/R25+R24)
// and an RC low-pass filter (R25/C24).
//
// Input:  Voltage at node N017.
// Output: Voltage at node N030.
//
// SPICE components included: C25, R25, R24, C24.
// The system has 2 unknown node voltages and no internal voltage sources.
class LeadPathCoupling : public MNASolver<2, 0> {
private:
    static constexpr int GND = -1;

    // The unknowns are the voltages at the internal nodes of the network.
    enum Var {
        V_N018 = 0, // After coupling cap C25
        V_N030      // After resistor R25 (output node)
    };

    // Component values
    double R24, R25;
    double C24, C25;
    
    // State variables for the capacitors
    std::vector<double> cap_z_state;

    // This stage is linear, so we can solve it in one step.
    void solve(double in) {
        resetMatrices();
        
        // --- Stamp static components ---
        stampResistor(V_N030, V_N018, R25);
        stampResistor(V_N030, GND, R24);

        // --- Stamp capacitors ---
        // C24 is a simple shunt capacitor to ground
        stampCapacitor(V_N030, GND, C24, cap_z_state[0]);
        
        // C25 is the input coupling capacitor. Its input node is the known voltage 'in'.
        double Gc_in = invT_2 * C25;
        double Ic_hist_in = cap_z_state[1];
        
        // Stamp conductance Gc at the V_N018 node
        A[V_N018][V_N018] += Gc_in;
        
        // Stamp the equivalent history/input current source
        b[V_N018] += Ic_hist_in + Gc_in * in;
        
        // --- Solve the linear system ---
        if (lu_decompose()) {
            lu_solve(x);
        } else {
            x.fill(0.0); // Should not happen in a passive network
        }
    }

public:
    LeadPathCoupling() : cap_z_state(2, 0.0) {
        // Set fixed component values from SPICE netlist
        R24 = 68e3;
        R25 = 270e3;
        C24 = 1000e-12; // 1000pF
        C25 = 0.022e-6; // .022uF
    }

    // No user-adjustable parameters in this stage.
    
    double process(double in) {
        // Solve the linear system for the current input 'in'
        solve(in);
        
        // Update capacitor states for the next time step.
        updateCapacitorState(x[V_N030], 0, C24, cap_z_state[0]);
        updateCapacitorState(x[V_N018], in, C25, cap_z_state[1]);

        // The output of this stage is the voltage at node N030
        return x[V_N030];
    }
};
