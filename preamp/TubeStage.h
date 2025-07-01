#pragma once

class TubeStage {
public:
    TubeStage();
    void prepare(double sampleRate, double R_k, double R_L, double V_supply, double C_k);
    void reset();
    double process(double in, double R_load);
    double getVpDC() const { return Vp_dc; }

private:
    void calculateOperatingPoint();

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
};
