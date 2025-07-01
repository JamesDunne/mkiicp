#pragma once
#include "IIRBiquad.h"

class TubeStage {
public:
    TubeStage();
    void prepare(double sampleRate, double R_k, double R_L, double V_supply, double C_k);
    void reset();
    double process(double in, double R_load);

    // Accessor to get the calculated DC plate voltage
    double getVpDC() const { return Vp_dc; }

private:
    void calculateOperatingPoint();

    // Physical Parameters
    double R_k, R_L, V_supply;
    bool isBypassed;

    // SPICE Model Constants
    const double Mu = 96.20, Ex = 1.437, KG1 = 613.4;

    // Calculated Operating Point & Thevenin Parameters
    double Vp_dc;
    double r_p;
    double R_out;
    double gain;

    // Filters
    IIRBiquad cathodeBypassFilter;

public:
    // Public filters for external coupling networks
    IIRBiquad inputFilter, outputFilter, interStageLPF;
};
