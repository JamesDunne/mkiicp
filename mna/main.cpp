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
        // std::cout << x << std::endl;
    }

    converged = diverged = failed = 0;

#if 0
    std::cout << x << std::endl;
    for (int i = 0; i < 10; i++) {
        x = preamp.processSample(0.0);
        std::cout << x << std::endl;
    }
#endif

    auto start = std::chrono::high_resolution_clock::now();

    double min = 1.0, max = -1.0;
    processWavFile(
        argv[1],
        argv[2],
        [&](double sample, long time) -> double {
            double out = preamp.processSample(sample);
            if (out > max) max = out;
            if (out < min) min = out;
            // std::cout << std::fixed << std::setprecision(6) << out << std::endl;
            return out / 410.0;
        }
    );

    auto end = std::chrono::high_resolution_clock::now();

    std::cout << "time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << std::endl;

    std::cout << min << " " << max << std::endl;

    Triode::printStats();
    std::cout << converged << "," << diverged << "," << failed << std::endl;
    std::cout << ((double)diverged / (double)(diverged + converged + failed) * 100.0) << "% diverged" << std::endl;
    std::cout << ((double)failed / (double)(diverged + converged + failed) * 100.0) << "% failed" << std::endl;

    return 0;
}