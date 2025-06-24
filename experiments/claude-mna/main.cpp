
#include <iostream>

#include "mkiicp.hpp"

int main(int argc, char *argv[]) {
    const double sampleRate = 48000.0;
    MKIICSimulator amp { sampleRate };

    // Set tone controls
    amp.setVolume1(0.75);
    amp.setTreble(0.8);
    amp.setMid(0.33);
    amp.setBass(0.05);
    amp.setGain(0.5);
    amp.setMaster(0.5);

    double startFreq = 20.0f;        // Start frequency in Hz
    double endFreq = 20000.0f;       // End frequency in Hz
    double duration = 5.0f;          // Sweep duration in seconds

    double k = (endFreq / startFreq);
    // float L = duration / std::log(k);

    double phase = 0.0;
    for (long i = 0; i < 48000 * 5; ++i) {
        double t = static_cast<double>(i) / sampleRate;
        double instantaneousFreq = startFreq * std::pow(k, t / duration);

        // Calculate phase increment for this sample
        double phaseIncrement = 2.0f * M_PI * instantaneousFreq / sampleRate;

        // Add to accumulated phase and generate sample
        phase += phaseIncrement;
        double x = std::sin(phase);

        // process sample with first gain stage and tone stack:
        x = amp.processSample(x);

        std::cout
            << t << ","
            << x
            << std::endl;
    }

    return 0;
}