#pragma once

#include "V1A_ToneStack.h"
#include "LeadAndMixerStage.h"
#include "V2B_Stage.h"
#include "V2A_OutputStage.h"

class MarkIICPreamp {
private:
    // Instantiate all 9 processing stages
    // This is a conceptual layout. The actual segments would be implemented
    // following the `V1A_ToneStack` example.
    V1A_ToneStack v1aToneStack;
    PostToneDriver leadAndMixer;
    V2B_Stage v2b;
    V2A_OutputStage v2aOutput;

public:
    MarkIICPreamp() {
        Triode::initializeLUT();
    }

    void setup(double sampleRate) {
        v1aToneStack.setup(sampleRate);
        leadAndMixer.setup(sampleRate);
        v2b.setup(sampleRate);
        v2aOutput.setup(sampleRate);
    }

    // Parameter setting methods
    void setVolume1(double val) { v1aToneStack.setVolume1(val); }
    void setTreble(double val) { v1aToneStack.setTreble(val); }
    void setBass(double val) { v1aToneStack.setBass(val); }
    void setMid(double val) { v1aToneStack.setMid(val); }
    void setGain(double val) { leadAndMixer.setGain(val); }
    void setMaster(double val) { v2aOutput.setMaster(val); }

    double processSample(double in) {
        double out_v1a = v1aToneStack.process(in);
        // return out_v1a;

        double out_lead = leadAndMixer.process(out_v1a);
        // return out_lead;

        double out_v2b = v2b.process(out_lead);
        // return out_v2b;

        // Final V2A output stage:
        double out_master = v2aOutput.process(out_v2b);
        return out_master;
    }
};
