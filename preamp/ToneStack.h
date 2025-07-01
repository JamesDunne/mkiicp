//
// Created by Jim Dunne on 6/30/25.
//

#pragma once

#include "IIRFilter.h" // Use the new filter header

class ToneStack {
public:
    ToneStack() = default;
    void prepare(double sampleRate);
    void reset();
    void setParams(double treble, double mid, double bass, double volume);
    double process(double in);

private:
    void calculateCoefficients();

    double sampleRate = 44100.0;
    // Store pot positions (0-1) directly
    double p_treble = 0.8, p_mid = 0.5, p_bass = 0.25, p_vol = 0.75;

    DirectFormIIR<3> filter;
    double volumeGain = 1.0;
};
