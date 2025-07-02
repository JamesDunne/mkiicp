#include <iostream>

#include "WaveFile.h"
#include "MarkIICPreamp.h"

int main(int argc, const char *argv[]) {
    if (argc <= 2) {
        std::cout << "Usage: " << argv[0] << "[infile.wave] [outfile.wav]" << std::endl;
        return 1;
    }

    MarkIICPreamp preamp;
    preamp.setup(48000.0);

    std::cout << "Hello, World!" << std::endl;

    double x = 0.0;
    for (int i = 0; i < 100000; i++) {
        x = preamp.processSample(0.0);
    }

#if 0
    std::cout << x << std::endl;
    for (int i = 0; i < 10; i++) {
        x = preamp.processSample(0.0);
        std::cout << x << std::endl;
    }
#endif

    double min = 1.0, max = -1.0;
    processWavFile(
        argv[1],
        argv[2],
        [&](double sample) -> double {
            double out = preamp.processSample(sample);
            if (out > max) max = out;
            if (out < min) min = out;
            return out / 410.0;
        }
    );

    std::cout << min << " " << max << std::endl;

    return 0;
}