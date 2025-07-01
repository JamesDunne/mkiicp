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

#if 0
    for (int i = 0; i < 1024; i++) {
        double x = (i - 512) / 512.0;
        double y = tubeSaturate12AX7(x, 1.0, 1.5, 210.0);
        std::cout << x << "\t" << y << std::endl;
    }
#endif

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

    std::cout << "v1a: ";
    myPreamp.mm_v1a.printMinMax();
    std::cout << "tone: ";
    myPreamp.mm_toneStack.printMinMax();
    std::cout << "out: ";
    myPreamp.mm_output.printMinMax();

    return 0;
}
