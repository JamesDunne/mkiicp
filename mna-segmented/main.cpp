#include "MkiicPreamp.h"
#include <iostream>
#include <vector>

// Dummy audio buffer for demonstration
std::vector<float> audio_buffer(48000, 0.0f);

int main() {
    const double sampleRate = 48000.0;

    // 1. Create an instance of the preamp
    MkiicPreamp preamp(sampleRate);
    std::cout << "MKIIC+ Preamp simulation initialized." << std::endl;

    // 2. Set parameters for a high-gain lead sound
    std::cout << "Setting parameters for a lead tone..." << std::endl;
    preamp.setToneParams(0.7, 0.2, 0.6, 0.8); // Treble, Mid, Bass, Vol1
    preamp.setLeadDrive(0.9);                  // High gain
    preamp.setMasterVolume(0.6);               // Master volume

    // 3. Create a test signal (e.g., an impulse)
    audio_buffer[10] = 1.0f;

    // 4. Process the audio
    std::cout << "Processing audio block..." << std::endl;
    preamp.processBlock(audio_buffer.data(), audio_buffer.size());

    std::cout << "Processing complete." << std::endl;
    // At this point, 'audio_buffer' contains the processed sound.
    // You could write it to a WAV file to listen to the impulse response.

    // --- Real-time control example ---
    // Imagine this is happening in a GUI thread while audio runs:
    std::cout << "\nChanging to a scooped rhythm tone in 'real-time'..." << std::endl;
    preamp.setToneParams(0.6, 0.05, 0.6, 0.9);
    preamp.setLeadDrive(0.2); // Low gain
    preamp.setMasterVolume(0.8);

    // Process another block of audio with the new settings
    // ...

    return 0;
}
