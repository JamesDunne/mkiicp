#pragma once

#include "MNASolver.h"
#include <vector>

// This class models the final coupling network of the lead channel before
// the signal is mixed back into the main path at node N002.
//
// Input:  Voltage at node N025.
// Output: Voltage at node N026.
//
// SPICE components included: C30, R31, C31.
//
// APPROXIMATION NOTE: R31 and C31 connect N026 to the mixer node N002.
// To create an independent processing block, we are treating N002 as ground
// for this stage. This makes R31/C31 a load for this stage. A more
// comprehensive simulation would combine this network into the main mixer MNA system.
class LeadPathCoupling2 : public MNASolver<1, 0> {
private:
    static constexpr int GND = -1;

    // The only unknown is the voltage at the output node.
    enum Var {
        V_N026 = 0
    };

    // Component values
    double R31;
    double C30, C31;
    
    // State variables for the capacitors
    std::vector<double> cap_z_state;

    // This stage is linear, so we can solve it in one step.
    void solve(double in) {
        resetMatrices();
        
        // --- Stamp the load components (connected to ground per approximation) ---
        stampResistor(V_N026, GND, R31);
        stampCapacitor(V_N026, GND, C31, cap_z_state[0]);
        
        // --- Stamp the input coupling capacitor C30 ---
        // Its input node is the known voltage 'in'.
        double Gc_in = invT_2 * C30;
        double Ic_hist_in = cap_z_state[1];
        
        // Stamp conductance Gc at the V_N026 node
        A[V_N026][V_N026] += Gc_in;
        
        // Stamp the equivalent history/input current source
        b[V_N026] += Ic_hist_in + Gc_in * in;
        
        // --- Solve the 1x1 linear system ---
        if (A[V_N026][V_N026] > 1e-12) {
             x[V_N026] = b[V_N026] / A[V_N026][V_N026];
        } else {
             x.fill(0.0);
        }
    }

public:
    LeadPathCoupling2() : cap_z_state(2, 0.0) {
        // Set fixed component values from SPICE netlist
        R31 = 220e3;
        C30 = 0.047e-6; // .047uF
        C31 = 250e-12;  // 250pF
    }

    // No user-adjustable parameters in this stage.
    
    double process(double in) {
        // Solve the linear system for the current input 'in'
        solve(in);
        
        // Update capacitor states for the next time step.
        updateCapacitorState(x[V_N026], 0, C31, cap_z_state[0]);
        updateCapacitorState(x[V_N026], in, C30, cap_z_state[1]);

        // The output of this stage is the voltage at node N026
        return x[V_N026];
    }
};
