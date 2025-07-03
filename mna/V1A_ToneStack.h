#pragma once
#include "MNASolver.h"
#include "Triode.h"
#include <vector>

class V1A_ToneStack : public MNASolver<10, 1> {
private:
    static constexpr int GND = -1;

    enum Var {
        V_N004 = 0, V_N033, V_N015, V_N005, V_N008, V_N007, V_N016, V_N028, V_N020, V_N019,
        I_VIN
    };

    // Component values
    double R1, R2, R4, R5, R5A;
    double C1, C2, C3, C4, C5, C6, C13B;

    // Potentiometers
    double RA_VOLUME1, RC_VOLUME1, RA_TREBLE, RC_TREBLE, RA_BASS, RA_MID;

    // CORRECTED: State variables for capacitors (z-domain history)
    std::array<double, 7> cap_z_state;

    const double VE = 405.0;

    void stampLinear() override {
        // Stamp the input signal source
        stampVoltageSource_A(V_N019, GND, 0);

        // Stamp all resistors, as they only change when a pot moves.
        // The power supply contribution is also linear.
        stampResistorLinear(V_N004, GND, R4);
        b_linear[V_N004] += VE / R4;
        stampResistorLinear(V_N033, GND, R2);
        stampResistorLinear(V_N019, GND, R1);
        stampResistorLinear(V_N004, V_N015, R5);
        stampResistorLinear(V_N008, V_N007, R5A);
        stampResistorLinear(V_N020, V_N008, RA_VOLUME1);
        stampResistorLinear(GND, V_N020, RC_VOLUME1);
        stampResistorLinear(V_N005, V_N007, RA_TREBLE);
        stampResistorLinear(V_N007, V_N016, RC_TREBLE);
        stampResistorLinear(V_N016, V_N028, RA_BASS);
        stampResistorLinear(V_N028, GND, RA_MID);

        // Stamp the constant, linear part of the capacitors
        stampCapacitor_A(V_N033, GND, C1);
        stampCapacitor_A(V_N033, GND, C2);
        stampCapacitor_A(V_N028, V_N015, C3);
        stampCapacitor_A(V_N016, V_N015, C4);
        stampCapacitor_A(V_N005, V_N004, C5);
        stampCapacitor_A(V_N005, V_N004, C6);
        stampCapacitor_A(V_N020, V_N008, C13B);
    }

    void stampDynamic(double in) override {
        // Stamp the input signal source
        stampVoltageSource_b(0, in);

        // Stamp only the dynamic history current of the capacitors
        stampCapacitor_b(V_N033, GND, cap_z_state[0]);      // C1
        stampCapacitor_b(V_N033, GND, cap_z_state[1]);      // C2
        stampCapacitor_b(V_N028, V_N015, cap_z_state[2]);   // C3
        stampCapacitor_b(V_N016, V_N015, cap_z_state[3]);   // C4
        stampCapacitor_b(V_N005, V_N004, cap_z_state[4]);   // C5
        stampCapacitor_b(V_N005, V_N004, cap_z_state[5]);   // C6
        stampCapacitor_b(V_N020, V_N008, cap_z_state[6]);   // C13B
    }

    void stampComponents(double inputVoltage) {
        resetMatrices();

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

        // --- Use new capacitor stamping method ---
        stampCapacitor(V_N033, GND, C1, cap_z_state[0]);
        stampCapacitor(V_N033, GND, C2, cap_z_state[1]);
        stampCapacitor(V_N028, V_N015, C3, cap_z_state[2]);
        stampCapacitor(V_N016, V_N015, C4, cap_z_state[3]);
        stampCapacitor(V_N005, V_N004, C5, cap_z_state[4]);
        stampCapacitor(V_N005, V_N004, C6, cap_z_state[5]);
        stampCapacitor(V_N020, V_N008, C13B, cap_z_state[6]);

        stampVoltageSource(V_N019, GND, 0, inputVoltage);
    }

    // stampNonLinear is unchanged
    void stampNonLinear(const std::array<double, NumUnknowns>& current_x) override {
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
    V1A_ToneStack() : cap_z_state() { // Changed name from cap_hist
        R1 = 1e6; R2 = 1.5e3; R4 = 150e3; R5 = 100e3; R5A = 100e3;
        C1 = 0.47e-6; C2 = 22e-6; C3 = 0.047e-6; C4 = 0.1e-6;
        C5 = 250e-12; C6 = 750e-12; C13B = 180e-12;

        setVolume1(0.75); setTreble(0.8); setBass(0.25); setMid(0.5);
    }

    void setVolume1(double val) { double v = val * val; RA_VOLUME1 = 1e6 * (1.0 - v); RC_VOLUME1 = 1e6 * v; setDirty(); }
    void setTreble(double val)  { double v = val * val; RA_TREBLE = 250e3 * (1.0-v); RC_TREBLE = 250e3 * v; setDirty(); }
    void setBass(double val)    { RA_BASS = 250e3 * (val * val); setDirty(); }
    void setMid(double val)     { RA_MID = 10e3 * (val * val); setDirty(); }

    double process(double in) {
        solveNonlinear(in);

        updateCapacitorState(x[V_N033], 0, C1, cap_z_state[0]);
        updateCapacitorState(x[V_N033], 0, C2, cap_z_state[1]);
        updateCapacitorState(x[V_N028], x[V_N015], C3, cap_z_state[2]);
        updateCapacitorState(x[V_N016], x[V_N015], C4, cap_z_state[3]);
        updateCapacitorState(x[V_N005], x[V_N004], C5, cap_z_state[4]);
        updateCapacitorState(x[V_N005], x[V_N004], C6, cap_z_state[5]);
        updateCapacitorState(x[V_N020], x[V_N008], C13B, cap_z_state[6]);

        return x[V_N020];
    }
};
