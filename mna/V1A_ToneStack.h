#pragma once
#include "MNASolver.h"
#include "Triode.h"
#include <vector>

// This stage combines the first tube stage (V1A) and the TMB tone stack.
// It models all components connected to nodes N019, N004, N033, N015, N005, N008, N007, N016, N028, N020.
// The system has 10 nodes and 1 voltage source (the input signal).
class V1A_ToneStack : public MNASolver<10, 1> {
private:
    // A single enum for all unknowns in the solution vector 'x'.
    // The first 10 are node voltages (0-based indices).
    // The last 1 is the current through the input voltage source.
    enum Var {
        // Node Voltages
        V_N004 = 0, V_N033, V_N015, V_N005, V_N008, V_N007, V_N016, V_N028, V_N020, V_N019,
        // V-Source Currents
        I_VIN // This will be at index 10 (NumNodes)
    };

    // Component values
    double R1, R2, R4, R5, R5A;
    double C1, C2, C3, C4, C5, C6, C13B;

    // Potentiometer resistances
    double RA_VOLUME1, RC_VOLUME1;
    double RA_TREBLE, RC_TREBLE;
    double RA_BASS, RA_MID;

    // History for capacitors
    std::vector<double> cap_hist;

    // Power supply
    const double VE = 405.0;

    void stampComponents(double inputVoltage) {
        resetMatrices();

        // --- Static Components ---
        // Power supply resistor to V1A Plate
        stampResistor(V_N004 + 1, 0, R4);
        b[V_N004] += VE / R4;

        // V1A Cathode Resistor
        stampResistor(V_N033 + 1, 0, R2);

        // Input Grid-Leak Resistor
        stampResistor(V_N019 + 1, 0, R1);

        // Tone Stack & Volume Pot Resistors
        stampResistor(V_N004 + 1, V_N015 + 1, R5);
        stampResistor(V_N008 + 1, V_N007 + 1, R5A);
        stampResistor(V_N020 + 1, V_N008 + 1, RA_VOLUME1);
        stampResistor(0, V_N020 + 1, RC_VOLUME1);
        stampResistor(V_N005 + 1, V_N007 + 1, RA_TREBLE);
        stampResistor(V_N007 + 1, V_N016 + 1, RC_TREBLE);
        stampResistor(V_N016 + 1, V_N028 + 1, RA_BASS);
        stampResistor(V_N028 + 1, 0, RA_MID);

        // --- Capacitors (and their history) ---
        stampCapacitor(V_N033 + 1, 0, C1, cap_hist[0]);
        stampCapacitor(V_N033 + 1, 0, C2, cap_hist[1]);
        stampCapacitor(V_N028 + 1, V_N015 + 1, C3, cap_hist[2]);
        stampCapacitor(V_N016 + 1, V_N015 + 1, C4, cap_hist[3]);
        stampCapacitor(V_N005 + 1, V_N004 + 1, C5, cap_hist[4]);
        stampCapacitor(V_N005 + 1, V_N004 + 1, C6, cap_hist[5]);
        stampCapacitor(V_N020 + 1, V_N008 + 1, C13B, cap_hist[6]);

        // --- Input Source ---
        // V1 in SPICE: Connects node N019 to ground.
        stampVoltageSource(V_N019 + 1, 0, 0, inputVoltage);
    }

    void stampNonLinear(const std::array<double, NumUnknowns>& current_x) {
        // V1A (XV1A): Plate=N004, Grid=N019, Cathode=N033
        double v_p = current_x[V_N004];
        double v_g = current_x[V_N019];
        double v_c = current_x[V_N033];

        Triode::State ts = Triode::calculate(v_p - v_c, v_g - v_c);

        // --- Plate Current ---
        // Linearized current source I_lin = I(v0) - g_p*v_p - g_g*v_g
        double i_p_lin = ts.ip - ts.g_p * (v_p - v_c) - ts.g_g * (v_g - v_c);
        stampCurrentSource(V_N004 + 1, V_N033 + 1, i_p_lin);

        // Plate-cathode conductance
        stampConductance(V_N004 + 1, V_N033 + 1, ts.g_p);

        // Transconductance (V_g influence on I_p)
        // KCL at plate (V_N004): +g_g*(v_g - v_c) -> A[p][g]+=g_g, A[p][c]-=g_g
        A[V_N004][V_N019] += ts.g_g;
        A[V_N004][V_N033] -= ts.g_g;
        // KCL at cathode (V_N033): -g_g*(v_g - v_c) -> A[c][g]-=g_g, A[c][c]+=g_g
        A[V_N033][V_N019] -= ts.g_g;
        A[V_N033][V_N033] += ts.g_g;

        // --- Grid Current ---
        // Linearized current source from grid to cathode
        double i_g_lin = ts.ig - ts.g_ig * (v_g - v_c);
        stampCurrentSource(V_N019 + 1, V_N033 + 1, i_g_lin);
        stampConductance(V_N019 + 1, V_N033 + 1, ts.g_ig);
    }

public:
    V1A_ToneStack() : cap_hist(7, 0.0) {
        // Set fixed component values from SPICE netlist
        R1 = 1e6; R2 = 1.5e3; R4 = 150e3; R5 = 100e3; R5A = 100e3;
        C1 = 0.47e-6; C2 = 22e-6; C3 = 0.047e-6; C4 = 0.1e-6;
        C5 = 250e-12; C6 = 750e-12; C13B = 180e-12;

        // Set default pot values based on .param
        setVolume1(0.75);
        setTreble(0.8);
        setBass(0.25);
        setMid(0.5);
    }

    void setVolume1(double val) { double v = val * val; RA_VOLUME1 = 1e6 * (1.0-v); RC_VOLUME1 = 1e6 * v; }
    void setTreble(double val)  { double v = val * val; RA_TREBLE = 250e3 * (1.0-v); RC_TREBLE = 250e3 * v; }
    void setBass(double val)    { RA_BASS = 250e3 * (val * val); }
    void setMid(double val)     { RA_MID = 10e3 * (val * val); }

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
        x = current_x; // Store solution for next sample and for history update

        // Update capacitor history using the new voltages
        updateCapacitorHistory(x[V_N033], 0,         cap_hist[0]); // C1
        updateCapacitorHistory(x[V_N033], 0,         cap_hist[1]); // C2
        updateCapacitorHistory(x[V_N028], x[V_N015], cap_hist[2]); // C3
        updateCapacitorHistory(x[V_N016], x[V_N015], cap_hist[3]); // C4
        updateCapacitorHistory(x[V_N005], x[V_N004], cap_hist[4]); // C5
        updateCapacitorHistory(x[V_N005], x[V_N004], cap_hist[5]); // C6
        updateCapacitorHistory(x[V_N020], x[V_N008], cap_hist[6]); // C13B

        return x[V_N020]; // Output is the voltage at the volume pot wiper (node N020)
    }
};
