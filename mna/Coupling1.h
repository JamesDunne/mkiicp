#pragma once

#include "MNASolver.h"
#include <vector>

// This class models the passive coupling network that feeds the lead channel.
// It includes the coupling capacitor C21 and the "Lead Drive" potentiometer.
//
// Input:  Voltage at node N001 (from V1B_Stage).
// Output: Voltage at node N029 (the wiper of the Lead Drive pot).
//
// SPICE components included: C21, R21, RA_LEAD_DRIVE, RC_LEAD_DRIVE, R22.
// The system has 3 unknown node voltages and no internal voltage sources.
class Coupling1 : public MNASolver<3, 0> {
private:
    static constexpr int GND = -1;

    // The unknowns are the voltages at the internal nodes of the network.
    enum Var {
        V_N011 = 0, // After coupling cap C21
        V_N024,     // Top of the Lead Drive pot
        V_N029      // Wiper of the Lead Drive pot (output)
    };

    // Component values
    double R21, R22;
    double C21;
    
    // Potentiometer resistances
    double RA_LEAD_DRIVE, RC_LEAD_DRIVE;
    
    // State variable for the capacitor
    std::vector<double> cap_z_state;

    // This stage is linear, so we can solve it in one step.
    void solve(double in) {
        resetMatrices();
        
        // --- Stamp static components ---
        stampResistor(V_N011, V_N024, R21);
        stampResistor(V_N029, GND, R22);

        // --- Stamp potentiometer ---
        stampResistor(V_N024, V_N029, RA_LEAD_DRIVE);
        stampResistor(V_N029, GND, RC_LEAD_DRIVE);
        
        // --- Stamp capacitor C21 (connected between input 'in' and node V_N011) ---
        // This is a special case since one side is connected to a known voltage 'in'
        // instead of another unknown node in our system.
        double Gc = invT_2 * C21;
        double Ic_hist = cap_z_state[0];
        
        // Stamp conductance Gc at the V_N011 node
        A[V_N011][V_N011] += Gc;
        
        // Stamp the equivalent current source for the capacitor.
        // I_eq = I_history + Gc * v_in
        b[V_N011] += Ic_hist + Gc * in;
        
        // --- Solve the linear system ---
        if (lu_decompose()) {
            lu_solve(x);
        } else {
            x.fill(0.0); // Should not happen in a passive network
        }
    }

public:
    Coupling1() : cap_z_state(1, 0.0) {
        // Set fixed component values from SPICE netlist
        C21 = 0.02e-6;
        R21 = 680e3;
        R22 = 475e3;

        // Set default gain parameter based on .param
        setGain(0.75);
    }

    // Set the Lead Drive parameter (0-1)
    void setGain(double val) {
        // From SPICE: {1Meg*(1-(gain*gain))} and {1Meg*(gain*gain)}
        double p = val * val;
        RA_LEAD_DRIVE = 1e6 * (1.0 - p);
        RC_LEAD_DRIVE = 1e6 * p;
    }

    double process(double in) {
        // Solve the linear system for the current input 'in'
        solve(in);
        
        // Update capacitor state for the next time step.
        // The voltage across the cap is v_n11 - v_in.
        updateCapacitorState(x[V_N011], in, C21, cap_z_state[0]);

        // The output of this stage is the voltage at the pot wiper (node N029)
        return x[V_N029];
    }
};
