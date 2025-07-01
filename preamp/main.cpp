#include <iostream>
#include <vector>
#include "Preamp.h"

// Dummy audio data for demonstration
const int SAMPLE_RATE = 48000;
const int DURATION_SECONDS = 5;
const int NUM_SAMPLES = SAMPLE_RATE * DURATION_SECONDS;

int main() {
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

    // 3. Create a dummy input signal (e.g., an impulse)
    std::vector<double> inputSignal(NUM_SAMPLES, 0.0);
    inputSignal[10] = 1.0; // Impulse to check the response

    // 4. Process the audio buffer
    std::vector<double> outputSignal(NUM_SAMPLES);

    std::cout << "Processing " << DURATION_SECONDS << " seconds of audio..." << std::endl;

    for (int i = 0; i < NUM_SAMPLES; ++i) {
        outputSignal[i] = myPreamp.processSample(inputSignal[i]);
    }

    std::cout << "Processing complete." << std::endl;
    std::cout << "You can now write the 'outputSignal' vector to a .wav file to hear the result." << std::endl;

    // Example: Print first few output samples
    for (int i = 10; i < 30; ++i) {
        std::cout << "Sample " << i << ": " << outputSignal[i] << std::endl;
    }

    return 0;
}
