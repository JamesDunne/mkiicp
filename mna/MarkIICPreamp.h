#pragma once

#include "Coupling1.h"
#include "MixerAndV2B.h"
#include "V1A_ToneStack.h"
#include "V1B_Stage.h"
#include "V2A_OutputStage.h"
#include "V3B_Stage.h"
#include "V4A_Stage.h"

class MarkIICPreamp {
private:
    // Instantiate all 9 processing stages
    // This is a conceptual layout. The actual segments would be implemented
    // following the `V1A_ToneStack` example.
    V1A_ToneStack v1aToneStack;
    V1B_Stage v1b;
    V3B_and_IO_Coupling v3b;
    V4A_and_MixerLoad v4a;
    MixerAndV2B mixerv2b;
    V2A_OutputStage v2aOutput;

public:
    MarkIICPreamp() {
        Triode::initializeLUT();
    }

    void setup(double sampleRate) {
        v1aToneStack.setup(sampleRate);
        v1b.setup(sampleRate);
        v3b.setup(sampleRate);
        v4a.setup(sampleRate);
        mixerv2b.setup(sampleRate);
        v2aOutput.setup(sampleRate);
    }

    // Parameter setting methods
    void setVolume1(double val) { v1aToneStack.setVolume1(val); }
    void setTreble(double val) { v1aToneStack.setTreble(val); }
    void setBass(double val) { v1aToneStack.setBass(val); }
    void setMid(double val) { v1aToneStack.setMid(val); }
    void setGain(double val) { v3b.setGain(val); }
    void setMaster(double val) { v2aOutput.setMaster(val); }

    double processSample(double in) {
        double out_v1a = v1aToneStack.process(in);
        // return out_v1a;

        double out_v1b = v1b.process(out_v1a);
        // return out_v1b;

        double out_v3b = v3b.process(out_v1b);
        // return out_v3b;

        double out_v4a = v4a.process(out_v3b);
        // return out_v4a;

        // Mixer stage:
        double out_v2b = mixerv2b.process(out_v4a, out_v1b);
        // return out_v2b;

        // Final V2A output stage:
        double out_master = v2aOutput.process(out_v2b);
        return out_master;
    }
};
