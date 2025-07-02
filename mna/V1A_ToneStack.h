#pragma once
#include "MNASolver.h"
#include "Triode.h"
#include <vector>

class V1A_ToneStack : public MNASolver<10, 1> {
private:
    static constexpr int GND = -1; // Define Ground node for clarity

    enum Var {
        V_N004 = 0, V_N033, V_N015, V_N005, V_N008, V_N007, V_N016, V_N028, V_N020, V_N019,
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
        stampResistor(V_N004, GND, R4);
        b[V_N004] += VE / R4;
        stampResistor(V_N033, GND, R2);
        stampResistor(V_N019, GND, R1);
        stampResistor(V_N004, V_N015, R5);
        stampResistor(V_N008, V_N007, R5A);
        stampResistor(V_N020, V_N008, RA_VOLUME1);
        stampResistor(GND, V_N020, RC_VOLUME1);
        stampResistor(V_N005, V_N007, RA_TREBLE);
        stampResistor(V_N007, V_N016, RC_TREBLE);
        stampResistor(V_N016, V_N028, RA_BASS);
        stampResistor(V_N028, GND, RA_MID);

        // --- Capacitors ---
        stampCapacitor(V_N033, GND, C1, cap_hist[0]);
        stampCapacitor(V_N033, GND, C2, cap_hist[1]);
        stampCapacitor(V_N028, V_N015, C3, cap_hist[2]);
        stampCapacitor(V_N016, V_N015, C4, cap_hist[3]);
        stampCapacitor(V_N005, V_N004, C5, cap_hist[4]);
        stampCapacitor(V_N005, V_N004, C6, cap_hist[5]);
        stampCapacitor(V_N020, V_N008, C13B, cap_hist[6]);

        // --- Input Source ---
        stampVoltageSource(V_N019, GND, 0, inputVoltage);
    }

    void stampNonLinear(const std::array<double, NumUnknowns>& current_x) {
        double v_p = current_x[V_N004];
        double v_g = current_x[V_N019];
        double v_c = current_x[V_N033];

        Triode::State ts = Triode::calculate(v_p - v_c, v_g - v_c);

        double i_p_lin = ts.ip - ts.g_p * (v_p - v_c) - ts.g_g * (v_g - v_c);
        stampCurrentSource(V_N004, V_N033, i_p_lin);
        stampConductance(V_N004, V_N033, ts.g_p);

        A[V_N004][V_N019] += ts.g_g;
        A[V_N004][V_N033] -= ts.g_g;
        A[V_N033][V_N019] -= ts.g_g;
        A[V_N033][V_N033] += ts.g_g;

        double i_g_lin = ts.ig - ts.g_ig * (v_g - v_c);
        stampCurrentSource(V_N019, V_N033, i_g_lin);
        stampConductance(V_N019, V_N033, ts.g_ig);
    }

public:
    V1A_ToneStack() : cap_hist(7, 0.0) {
        R1 = 1e6; R2 = 1.5e3; R4 = 150e3; R5 = 100e3; R5A = 100e3;
        C1 = 0.47e-6; C2 = 22e-6; C3 = 0.047e-6; C4 = 0.1e-6;
        C5 = 250e-12; C6 = 750e-12; C13B = 180e-12;

        setVolume1(0.75);
        setTreble(0.8);
        setBass(0.25);
        setMid(0.5);
    }

    void setVolume1(double val) { double v = val * val; RA_VOLUME1 = 1e6 * (1.0 - v); RC_VOLUME1 = 1e6 * v; }
    void setTreble(double val)  { double v = val * val; RA_TREBLE = 250e3 * (1.0-v); RC_TREBLE = 250e3 * v; }
    void setBass(double val)    { RA_BASS = 250e3 * (val * val); }
    void setMid(double val)     { RA_MID = 10e3 * (val * val); }

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

        updateCapacitorHistory(x[V_N033], 0,         cap_hist[0]);
        updateCapacitorHistory(x[V_N033], 0,         cap_hist[1]);
        updateCapacitorHistory(x[V_N028], x[V_N015], cap_hist[2]);
        updateCapacitorHistory(x[V_N016], x[V_N015], cap_hist[3]);
        updateCapacitorHistory(x[V_N005], x[V_N004], cap_hist[4]);
        updateCapacitorHistory(x[V_N005], x[V_N004], cap_hist[5]);
        updateCapacitorHistory(x[V_N020], x[V_N008], cap_hist[6]);

        return x[V_N020];
    }
};