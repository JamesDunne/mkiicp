#pragma once

#include <iomanip>
#include <iostream>
#include <ostream>

#include "ToneStack.h"
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
public:
    double lead_drive, master_vol;

    // --- Simulation Components ---
    TubeStage v1a, v1b, v3b, v4a, v2a, v2b;
    ToneStack toneStack;

    // --- Inter-stage Filter Networks ---
    // Each filter is named by the stage it follows.
    IIRBiquad v1a_pre_ts_lpf;  // LPF before tone stack
    IIRBiquad v1a_bypass_hpf;  // Simulates V1A's cathode bypass gain boost

    IIRBiquad v1b_bypass_hpf;  // Simulates V1B's bypass
    IIRBiquad v1b_to_rhythm_lpf; // R10/C10 filter

    IIRBiquad ts_to_v3b_hpf;   // C7/R9/C21/R21 network
    IIRBiquad v3b_bypass_hpf;  // Simulates V3B's bypass
    IIRBiquad v3b_to_v4a_hpf;
    IIRBiquad v3b_to_v4a_lpf;

    IIRBiquad v4a_bypass_hpf;
    IIRBiquad v4a_to_mixer_hpf;

    IIRBiquad mixer_output_lpf; // R11/C11/etc. filter
    IIRBiquad v2a_bypass_hpf;

    IIRBiquad v2b_output_hpf;

public:
    MinMax mm_v1a, mm_toneStack, mm_v1b, mm_v3b, mm_v4a, mm_v2a, mm_output;
};
