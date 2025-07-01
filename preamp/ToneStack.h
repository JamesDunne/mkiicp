//
// Created by Jim Dunne on 6/30/25.
//

#pragma once

#include "IIRBiquad.h" // We'll go back to the trusty biquad

// ... other class definitions ...

class ToneStack {
public:
    ToneStack();
    void prepare(double sampleRate);
    void reset();
    void setParams(double treble, double mid, double bass, double volume);
    double process(double in);

private:
    void calculateCoefficients();

    double sampleRate;
    double p_treble, p_mid, p_bass, p_vol;

    IIRBiquad bassMidFilter;   // Models the main body of the tone
    IIRBiquad trebleFilter;    // Models the bright cap and treble control
    double volumeGain;
};
