#pragma once
#include "IIRBiquad.h"

class TubeStage {
public:
    TubeStage();
    // New signature: Now takes bypass capacitance C_k. 0.0 means unbypassed.
    void prepare(double sampleRate, double R_k, double R_L, double V_supply, double C_k);
    void reset();
    double process(double in);

private:
    void solveDC();

    // Physical Parameters
    double R_k, R_L, V_supply, C_k;
    bool isBypassed;

    // SPICE Model Constants
    const double Mu = 96.20, Ex = 1.437, KG1 = 613.4, KP = 740.3, KVB = 1672.0;

    // Solved DC Operating Point
    double Vp_dc;
    double Ip_q;

    // Filters managed by the stage itself
    IIRBiquad cathodeBypassFilter;

public:
    // Public filters for inter-stage coupling, configured externally
    IIRBiquad inputFilter, outputFilter, interStageLPF;
};
