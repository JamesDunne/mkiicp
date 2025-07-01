#pragma once
#include "IIRBiquad.h"

class TubeStage {
public:
    TubeStage();
    void prepare(double sampleRate, double R_k, double R_L, double V_supply);
    void reset();
    double process(double in);

private:
    void solveDC(); // Solves for the quiescent operating point

    // --- Physical Parameters ---
    double R_k, R_L, V_supply;

    // --- SPICE Model Constants ---
    const double Mu = 96.20, Ex = 1.437, KG1 = 613.4, KP = 740.3, KVB = 1672.0;

    // --- Solved DC Operating Point ---
    double Vp_dc; // Quiescent Plate Voltage
    double Vk_dc; // Quiescent Cathode Voltage
    double Ip_q;  // Quiescent Plate Current

    // --- Real-time State ---
    double last_Ip; // Plate current from the previous sample

public:
    // Filters are now public to be configured by the main Preamp class
    IIRBiquad inputFilter, outputFilter, cathodeBypass, interStageLPF;
};
