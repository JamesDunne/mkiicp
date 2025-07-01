#pragma once
#include "IIRBiquad.h"

class TubeStage {
public:
    TubeStage();
    // Final signature: Takes sampleRate and all physical components.
    void prepare(double sampleRate, double R_k, double R_L, double V_supply, double C_k);
    void reset();
    // Process no longer needs R_load, as the solver handles it implicitly.
    double process(double in);
    double getVpDC() const { return Vp_dc; }

private:
    void solveDC();

    // Physical Parameters
    double R_k, R_L, V_supply;
    bool isBypassed;

    // SPICE Model Constants
    const double Mu = 96.20;

    // Calculated DC Operating Point
    double Vp_dc;
    double Ip_q; // Quiescent (DC) Plate Current

    // Internal stateful filter for the cathode voltage
    IIRBiquad cathodeFilter;
};
