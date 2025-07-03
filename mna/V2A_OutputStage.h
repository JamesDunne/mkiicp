#pragma once

#include "MNASolver.h"
#include "Triode.h"
#include <vector>
#include <algorithm>
#include <cmath>

// This class models the final tube stage (V2A), which is configured as a
// cathode follower buffer and includes the Master Volume control.
//
// Input:  Voltage at node N023.
// Output: Voltage at the master volume wiper, node N014.
//
// Includes a dampened Newton-Raphson solver for robustness.
class V2A_OutputStage : public MNASolver<6, 2> {
private:
    static constexpr int GND = -1;
    enum Var {
        V_N012 = 0, // V2A Plate
        V_N023,     // V2A Grid (Input)
        V_N031,     // V2A Cathode
        V_N006,     // Power Supply Node
        V_N013,     // Intermediate coupling node
        V_N014      // Final Output node (Master Volume wiper)
        // V-source currents I_Vin and I_Vsupply at indices 6 and 7
    };

    // Component values
    double R19, R104, R106;
    double C12, C15, C16;
    double RA_MASTER, RC_MASTER; // Potentiometer resistances

    std::vector<double> cap_z_state;
    const double VC2 = 410.0;

    void stampComponents(double inputVoltage) {
        resetMatrices();
        
        stampVoltageSource(V_N023, GND, 0, inputVoltage);
        stampVoltageSource(V_N006, GND, 1, VC2);

        // --- V2A Stage Components ---
        stampResistor(V_N012, V_N006, R19); // Plate resistor
        stampResistor(V_N031, GND, R104);   // Cathode resistor
        stampCapacitor(V_N031, GND, C15, cap_z_state[0]);
        stampCapacitor(V_N031, GND, C16, cap_z_state[1]);

        // --- Output & Master Volume Network ---
        // SPICE shows RA_MASTER from N014 to 0. This implements a standard
        // voltage divider pot for the master volume.
        stampCapacitor(V_N013, V_N012, C12, cap_z_state[2]);
        stampResistor(V_N013, V_N014, R106);
        stampResistor(V_N014, GND, RA_MASTER); // Correctly represents the pot
    }
    
    void stampNonLinear(const std::array<double, NumUnknowns>& current_x) {
        // XV2A: Plate=N012, Grid=N023, Cathode=N031
        double v_p = current_x[V_N012];
        double v_g = current_x[V_N023];
        double v_c = current_x[V_N031];
        Triode::State ts = Triode::calculate(v_p - v_c, v_g - v_c);
        
        double i_p_lin = ts.ip - ts.g_p*(v_p-v_c) - ts.g_g*(v_g-v_c);
        stampCurrentSource(V_N012, V_N031, i_p_lin);
        stampConductance(V_N012, V_N031, ts.g_p);
        
        A[V_N012][V_N023] += ts.g_g; A[V_N012][V_N031] -= ts.g_g;
        A[V_N031][V_N023] -= ts.g_g; A[V_N031][V_N031] += ts.g_g;

        double i_g_lin = ts.ig - ts.g_ig*(v_g-v_c);
        stampCurrentSource(V_N023, V_N031, i_g_lin);
        stampConductance(V_N023, V_N031, ts.g_ig);
    }

public:
    V2A_OutputStage() : cap_z_state(3, 0.0) {
        R19 = 120e3;
        R104 = 1e3;
        R106 = 15e3;
        C12 = 0.047e-6;
        C15 = 15e-6;
        C16 = 0.47e-6;
        
        setMaster(0.5); // Default value from .param
    }
    
    void setMaster(double val) {
        // The netlist shows a 1M pot on N014. We'll model it as a simple variable resistor to ground.
        // RA_MASTER N014 0 {1Meg*(master*master)}
        double p = val * val;
        RA_MASTER = 1e6 * p;
        if (RA_MASTER < 1e-3) RA_MASTER = 1e-3; // Prevent division by zero
    }

    double process(double in) {
        const int MAX_ITER = 25;
        const double CONVERGENCE_THRESH = 1e-6;
        const double DAMPING_LIMIT = 1.0; 

        std::array<double, NumUnknowns> current_x = x;
        for (int i = 0; i < MAX_ITER; ++i) {
            stampComponents(in);
            stampNonLinear(current_x);
            
            if (!lu_decompose()) { return 0.0; } 
            
            std::array<double, NumUnknowns> next_x;
            lu_solve(next_x);

            double max_delta = 0.0;
            for (size_t j = 0; j < NumUnknowns; ++j) {
                max_delta = std::max(max_delta, std::abs(next_x[j] - current_x[j]));
            }

            if (max_delta < CONVERGENCE_THRESH) {
                x = next_x;
                break;
            }

            if (max_delta > DAMPING_LIMIT) {
                double scale = DAMPING_LIMIT / max_delta;
                for (size_t j = 0; j < NumUnknowns; ++j) {
                    current_x[j] += scale * (next_x[j] - current_x[j]);
                }
            } else {
                current_x = next_x;
            }
            x = current_x;
        }
        
        updateCapacitorState(x[V_N031], 0, C15, cap_z_state[0]);
        updateCapacitorState(x[V_N031], 0, C16, cap_z_state[1]);
        updateCapacitorState(x[V_N013], x[V_N012], C12, cap_z_state[2]);

        // Final output is the voltage at the master volume wiper
        return x[V_N014];
    }
};
