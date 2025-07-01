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

    // let it settle:
    for (int i = 0; i < 10000; i++) {
        myPreamp.processSample(0);
    }
    myPreamp.mm_v1a.reset();
    myPreamp.mm_toneStack.reset();
    myPreamp.mm_v1b.reset();
    myPreamp.mm_v3b_in.reset();
    myPreamp.mm_v3b_out.reset();
    myPreamp.mm_v4a_in.reset();
    myPreamp.mm_v4a_out.reset();
    myPreamp.mm_v2b_in.reset();
    myPreamp.mm_v2b_out.reset();
    myPreamp.mm_v2a_in.reset();
    myPreamp.mm_v2a_out.reset();
    myPreamp.mm_output.reset();
#if 0
    for (int i = 0; i < 1024; i++) {
        double x = (i - 512) / 512.0;
        double y = myPreamp.processSample(x);
        std::cout << x << "\t" << y << std::endl;
    }
#else

    processWavFile(
        argv[1],
        argv[2],
        [&](double sample) -> double {
            return myPreamp.processSample(sample);
        }
    );

    myPreamp.mm_v1a.printMinMax(" v1a");
    myPreamp.mm_toneStack.printMinMax("tone");
    myPreamp.mm_v1b.printMinMax(" v1b");
    myPreamp.mm_v3b_in.printMinMax("v3bI");
    myPreamp.mm_v3b_out.printMinMax("v3bO");
    myPreamp.mm_v4a_in.printMinMax("v4aI");
    myPreamp.mm_v4a_out.printMinMax("v4aO");
    myPreamp.mm_v2b_in.printMinMax("v2bI");
    myPreamp.mm_v2b_out.printMinMax("v2bO");
    myPreamp.mm_v2a_in.printMinMax("v2aI");
    myPreamp.mm_v2a_out.printMinMax("v2aO");
    myPreamp.mm_output.printMinMax(" out");
#endif

    return 0;
}
