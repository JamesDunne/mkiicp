#pragma once

#include "Coupling1.h"
// #include "LeadPathCoupling.h"
#include "LeadPathCoupling2.h"
#include "V1A_ToneStack.h"
#include "V1B_Stage.h"
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
    V4A_and_Coupling v4a;
    LeadPathCoupling2 v4aCoupling;
    // MixerAndV2B stage8;
    // V2A_OutputStage stage9;

public:
    MarkIICPreamp() = default;

    void setup(double sampleRate) {
        v1aToneStack.setup(sampleRate);
        v1b.setup(sampleRate);
        v1bCoupling.setup(sampleRate);
        v3b.setup(sampleRate);
        v4a.setup(sampleRate);
        v4aCoupling.setup(sampleRate);
    }

    // Parameter setting methods
    void setVolume1(double val) { v1aToneStack.setVolume1(val); }
    void setTreble(double val) { v1aToneStack.setTreble(val); }
    void setBass(double val) { v1aToneStack.setBass(val); }
    void setMid(double val) { v1aToneStack.setMid(val); }
    void setGain(double val) { v1bCoupling.setGain(val); }
    // void setMaster(double val) { stage9.setMaster(val); }

    double processSample(double in) {
        double out_s1 = v1aToneStack.process(in);
        // return out_s1;

        double out_s2 = v1b.process(out_s1);
        // return out_s2;

        double out_s3 = v1bCoupling.process(out_s2); // out is N001
        // return out_s3;

        double lead_path_in = out_s3;
        double out_s4 = v3b.process(lead_path_in);
        // return out_s4;

        double out_s5 = v4a.process(out_s4);
        return out_s5;
        double out_s7 = v4aCoupling.process(out_s5); // out is N025
        // return out_s7;

        // double out_s8 = stage8.process(out_s3, out_s7); // Mixer stage
        // double out_s9 = stage9.process(out_s8);

        // The final E1 block is a simple gain stage
        // return out_s9 / 1400.0;
    }
};
