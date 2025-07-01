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
    double lead_drive;
    double master_vol;

    TubeStage v1a, v1b, v3b, v4a, v2a, v2b;
    ToneStack toneStack;
public:
    MinMax mm_v1a, mm_toneStack, mm_v1b, mm_v3b, mm_v4a, mm_v2a, mm_output;
};
