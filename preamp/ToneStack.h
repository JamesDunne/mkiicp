#pragma once
#include "IIRBiquad.h"

class ToneStack {
public:
    ToneStack();
    void prepare(double sampleRate, double v1a_Vp_dc);
    void reset();
    void setParams(double treble, double mid, double bass, double volume);
    double process(double in);

private:
    void calculateCoefficients();

    double sampleRate;
    double p_treble, p_mid, p_bass, p_vol;

    IIRBiquad filter;
    IIRBiquad dcBlocker;
    double vol_gain;
};
