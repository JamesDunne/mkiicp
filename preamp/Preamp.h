#pragma once
#include <iomanip>
#include <iostream>
#include <string>

#include "IIRBiquad.h"
#ifdef USE_TONESTACK_MNA
#include "ToneStackMNA.h"
using ToneStack = ToneStackFilter;
#else
#include "ToneStack.h"
#endif
#include "TubeStage.h"

class MinMax {
public:
    MinMax() : m_min(1.0e6), m_max(-1.0e6) {}

    void measureMinMax(double sample) {
        if (sample < m_min) m_min = sample;
        if (sample > m_max) m_max = sample;
    }

    void printMinMax(const std::string& name) {
        std::cout << name << ": "
            << std::setprecision(6) << m_min << " "
            << std::setprecision(6) << m_max << std::endl;
    }

    void reset() {
        m_min = 1.0e6;
        m_max = -1.0e6;
    }

private:
    double m_min,m_max;
};

// --- Main Preamp Class ---

class Preamp {
public:
    Preamp();
    void prepare(double sampleRate);
    void setParameters(double treble, double mid, double bass, double vol1, double gain, double master);
    double processSample(double in);

private:
    double lead_drive, master_vol;

    TubeStage v1a, v1b, v3b, v4a, v2a, v2b;
    ToneStack toneStack;

    // --- All inter-stage filters are managed here ---
    IIRBiquad v1a_output_lpf;
    IIRBiquad v1b_to_rhythm_lpf;
    IIRBiquad ts_to_v3b_hpf;
    IIRBiquad v3b_to_v4a_hpf;
    IIRBiquad v3b_to_v4a_lpf;
    IIRBiquad v4a_to_mixer_hpf;
    IIRBiquad mixer_output_lpf;
    IIRBiquad v2b_to_v2a_hpf;
    IIRBiquad v2b_to_v2a_lpf;
    IIRBiquad v2b_output_hpf;

public:
    MinMax mm_v1a, mm_toneStack, mm_v1b, mm_v3b_in, mm_v3b_out, mm_v4a_in, mm_v4a_out, mm_v2b_in, mm_v2b_out, mm_v2a_in, mm_v2a_out, mm_output;
};
