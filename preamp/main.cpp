#include <iostream>
#include <vector>
#include "Preamp.h"

// Dummy audio data for demonstration
const int SAMPLE_RATE = 48000;

void processWavFile(
    const std::string& inputFilename,
    const std::string& outputFilename,
    const std::function<double(double)>& processor
);

int main(int argc, char* argv[]) {
    if (argc <= 2) {
        std::cerr << "Usage: [infile] [outfile]" << std::endl;
        return 1;
    }

    // 1. Create and prepare the preamp instance
    Preamp myPreamp;
    myPreamp.prepare(SAMPLE_RATE);

    // 2. Set the control parameters (values from 0.0 to 1.0)
    // These correspond to the .param values in the SPICE file
    double treble = 0.8;
    double mid = 0.5;
    double bass = 0.25;
    double vol1 = 0.75;
    double gain = 0.75; // Lead Drive
    double master = 0.5;

    myPreamp.setParameters(treble, mid, bass, vol1, gain, master);

    processWavFile(
        argv[1],
        argv[2],
        [&](double sample) -> double {
            return myPreamp.processSample(sample);
        }
    );

    myPreamp.printMinMax();

    return 0;
}
