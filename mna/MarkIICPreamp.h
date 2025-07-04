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
    Coupling1 v1bCoupling;
    V3B_and_Coupling v3b;
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
        v1bCoupling.setup(sampleRate);
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
    void setGain(double val) { v1bCoupling.setGain(val); }
    void setMaster(double val) { v2aOutput.setMaster(val); }

    double processSample(double in) {
        double out_s1 = v1aToneStack.process(in);
        // return out_s1;

        double out_s2 = v1b.process(out_s1);
        // return out_s2;

        double out_s3 = v1bCoupling.process(out_s2); // out is N001
        // return out_s3;

        double lead_path_in = out_s3;
        double out_s4 = v3b.process(lead_path_in);
        return out_s4;

        double out_s5 = v4a.process(out_s4);
        // return out_s5;

        // Mixer stage:
        double out_s8 = mixerv2b.process(out_s5, out_s3);
        // return out_s8;

        // Final V2A output stage:
        double out_s9 = v2aOutput.process(out_s8);
        return out_s9;
    }
};
