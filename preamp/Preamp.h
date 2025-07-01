#pragma once

#include "IIRBiquad.h"
#include "ToneStack.h"
#include <iostream>
#include <ostream>

double tubeSaturate12AX7(double x, double drive, double bias, double Vp);

// --- Individual Stage Classes ---

class V1AStage {
public:
    V1AStage() = default;
    void prepare(double sampleRate);
    void reset();
    double process(double in);
private:
    IIRBiquad cathodeBypassFilter;
    IIRBiquad outputCouplingFilter;
    IIRBiquad interStageLPF;
    double gain = 50.0;
};

class V1BStage {
public:
    V1BStage() = default;
    void prepare(double sampleRate);
    void reset();
    double process(double in);
private:
    IIRBiquad cathodeBypassFilter;
    IIRBiquad outputCouplingFilter;
    IIRBiquad interStageLPF;
    double gain = 50.0;
};

class V3BV4AStage {
public:
    V3BV4AStage() = default;
    void prepare(double sampleRate);
    void reset();
    void setGain(double gain);
    double process(double in);
private:
    IIRBiquad fizzFilter;
    IIRBiquad v3b_inputFilter;
    IIRBiquad v3b_cathodeBypass;
    IIRBiquad interStageHPF;

    IIRBiquad v4a_cathodeBypass;
    IIRBiquad v4a_outputCoupling;

    double drive = 0.75;
};

class V2AStage {
public:
    V2AStage() = default;
    void prepare(double sampleRate);
    void reset();
    void setMaster(double master);
    double process(double in);
private:
    IIRBiquad inputCouplingFilter;
    IIRBiquad cathodeBypassFilter;
    IIRBiquad outputCouplingFilter;
    double gain = 50.0;
    double masterVol = 0.5;
};

class V2BStage {
public:
    V2BStage() = default;
    void prepare(double sampleRate);
    void reset();
    double process(double in);
private:
    IIRBiquad outputCouplingFilter;
    double gain = 50.0;
};

// --- Main Preamp Class ---

class Preamp {
public:
    Preamp() = default;

    void prepare(double sampleRate) {
        v1a.prepare(sampleRate);
        toneStack.prepare(sampleRate);
        v1b.prepare(sampleRate);
        v3b_v4a.prepare(sampleRate);
        v2a.prepare(sampleRate);
        v2b.prepare(sampleRate);
        reset(); // Ensure a clean state on preparation
    }

    void reset() {
        m_min = 1.0; m_max = -1.0;

        v1a.reset();
        toneStack.reset();
        v1b.reset();
        v3b_v4a.reset();
        v2a.reset();
        v2b.reset();
    }

    void setParameters(double treble, double mid, double bass, double vol1, double gain, double master) {
        toneStack.setParams(treble, mid, bass, vol1);
        v3b_v4a.setGain(gain);
        v2a.setMaster(master);
    }

    double processSample(double in);

    void measureMinMax(double sample) {
        if (sample < m_min) m_min = sample;
        if (sample > m_max) m_max = sample;
    }

    void printMinMax() {
        std::cout << m_min << " " << m_max << std::endl;
    }

private:
    V1AStage v1a;
    ToneStack toneStack;
    V1BStage v1b;
    V3BV4AStage v3b_v4a;
    V2AStage v2a;
    V2BStage v2b;

    double m_min,m_max;
};
