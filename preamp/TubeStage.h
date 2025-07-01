#pragma once
#include "IIRBiquad.h"

class TubeStage {
public:
    TubeStage();
    void prepare(double sampleRate, double R_k, double R_L, double V_supply, double C_k);
    void reset();
    double process(double in);
    double getVpDC() const { return Vp_dc; }

private:
    void solveDC();
    double getIp(double V_gk, double V_pk) const;

    // Physical Parameters
    double R_k, R_L, V_supply;
    bool isBypassed;

    // SPICE Model Constants
    const double Mu = 96.20;

    // Calculated Operating Point & Thevenin Parameters
    double Vp_dc;
    double r_p;        // Dynamic plate resistance
    double R_out;      // Thevenin equivalent output resistance
    double gain;       // Open-circuit AC gain

    double Ip_q;

    IIRBiquad cathodeFilter;
};
