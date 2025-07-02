#pragma once

#include <memory>
#include <vector>

// Forward declare all the circuit block classes
class InputAndToneStage;
class TubeStageV1B;
class Interstage_C7_R9;
class LeadDriveNetwork;
class TubeStageV3B;
class LeadChannelFinisher;
class MainMixer;
class TubeStageV2B;
class FinalStageAndMaster;

class MkiicPreamp {
public:
    MkiicPreamp(double sampleRate);
    ~MkiicPreamp(); // Needed for unique_ptr with forward-declared types

    // Sets the audio sample rate for all internal blocks
    void setSampleRate(double newSampleRate);

    // Main processing function for a single audio sample
    double process(double inputSample);
    
    // Processes a block of audio samples
    void processBlock(float* buffer, int numSamples);

    // --- Parameter Control ---

    // Tone stack and volume 1 controls
    void setToneParams(double treble, double mid, double bass, double vol1);

    // Lead channel gain
    void setLeadDrive(double gain);

    // Final master volume
    void setMasterVolume(double master);

private:
    // PIMPL (Pointer to Implementation) idiom is used here with unique_ptr
    // to hide the implementation details and dependencies from the header.

    // --- Circuit Segments ---
    std::unique_ptr<InputAndToneStage> block_inputAndTone;
    std::unique_ptr<TubeStageV1B> block_v1b;
    std::unique_ptr<Interstage_C7_R9> block_coupling1;
    
    // Lead Channel Path
    std::unique_ptr<LeadDriveNetwork> block_leadDrive;
    std::unique_ptr<TubeStageV3B> block_v3b;
    std::unique_ptr<LeadChannelFinisher> block_leadFinisher;

    // Mixer and Final Stages
    std::unique_ptr<MainMixer> block_mixer;
    std::unique_ptr<TubeStageV2B> block_v2b;
    std::unique_ptr<FinalStageAndMaster> block_finalStage;

    // Final output gain scaling from the SPICE model
    const double finalGain;
};
