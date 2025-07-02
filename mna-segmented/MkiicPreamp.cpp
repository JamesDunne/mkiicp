#include "MkiicPreamp.h"

// Include all the circuit block implementation headers
#include "InputAndToneStage.h"
#include "TubeStageV1B.h"
#include "Interstage_C7_R9.h"
#include "LeadDriveNetwork.h"
#include "TubeStageV3B.h"
#include "LeadChannelFinisher.h"
#include "MainMixer.h"
#include "TubeStageV2B.h"
#include "FinalStageAndMaster.h"


MkiicPreamp::MkiicPreamp(double sampleRate) 
    : finalGain(1.0 / 1400.0) // From E1 element in SPICE
{
    // Instantiate all the circuit blocks
    block_inputAndTone = std::make_unique<InputAndToneStage>(sampleRate);
    block_v1b = std::make_unique<TubeStageV1B>(sampleRate);
    block_coupling1 = std::make_unique<Interstage_C7_R9>(sampleRate);
    
    block_leadDrive = std::make_unique<LeadDriveNetwork>(sampleRate);
    block_v3b = std::make_unique<TubeStageV3B>(sampleRate);
    block_leadFinisher = std::make_unique<LeadChannelFinisher>(sampleRate);

    block_mixer = std::make_unique<MainMixer>(sampleRate);
    block_v2b = std::make_unique<TubeStageV2B>(sampleRate);
    block_finalStage = std::make_unique<FinalStageAndMaster>(sampleRate);

    // Set default initial parameters
    setToneParams(0.5, 0.5, 0.5, 0.75);
    setLeadDrive(0.5);
    setMasterVolume(0.5);
}

// The destructor needs to be defined in the .cpp file where the
// complete types of the unique_ptrs are known.
MkiicPreamp::~MkiicPreamp() = default;

void MkiicPreamp::setSampleRate(double newSampleRate) {
    // Propagate sample rate change to all blocks
    block_inputAndTone->setSampleRate(newSampleRate);
    block_v1b->setSampleRate(newSampleRate);
    block_coupling1->setSampleRate(newSampleRate);
    block_leadDrive->setSampleRate(newSampleRate);
    block_v3b->setSampleRate(newSampleRate);
    block_leadFinisher->setSampleRate(newSampleRate);
    block_mixer->setSampleRate(newSampleRate);
    block_v2b->setSampleRate(newSampleRate);
    block_finalStage->setSampleRate(newSampleRate);
}

void MkiicPreamp::setToneParams(double treble, double mid, double bass, double vol1) {
    block_inputAndTone->setParams(treble, mid, bass, vol1);
}

void MkiicPreamp::setLeadDrive(double gain) {
    block_leadDrive->setGain(gain);
}

void MkiicPreamp::setMasterVolume(double master) {
    block_finalStage->setMaster(master);
}

double MkiicPreamp::process(double inputSample) {
    // --- Define the Signal Flow ---
    // This method connects the output of one block to the input of the next.

    // 1. Input Stage (V1A) and Tone Stack
    double v_N020 = block_inputAndTone->process(inputSample);

    // 2. V1B Recovery Stage
    double v_N009 = block_v1b->process(v_N020);

    // 3. Coupling to mixer and lead channel
    double v_N001 = block_coupling1->process(v_N009);

    // 4. Process the parallel Lead Channel path
    double v_N029 = block_leadDrive->process(v_N001);
    double v_N017 = block_v3b->process(v_N029);
    double v_N026 = block_leadFinisher->process(v_N017);

    // 5. Mix the rhythm path (v_N001) and lead path (v_N026)
    double v_N002 = block_mixer->process(v_N001, v_N026);

    // 6. Post-mixer gain stage (V2B)
    double v_N021 = block_v2b->process(v_N002);

    // 7. Final stage (V2A) and master volume
    double v_N014 = block_finalStage->process(v_N021);

    // 8. Apply final output scaling
    return v_N014 * finalGain;
}

void MkiicPreamp::processBlock(float* buffer, int numSamples) {
    for (int i = 0; i < numSamples; ++i) {
        // Cast to double for internal processing to maintain precision
        double input = static_cast<double>(buffer[i]);
        double output = process(input);
        // Cast back to float for the output buffer
        buffer[i] = static_cast<float>(output);
    }
}
