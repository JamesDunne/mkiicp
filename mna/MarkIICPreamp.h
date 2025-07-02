#pragma once

#include "Coupling1.h"
#include "LeadPathCoupling.h"
#include "V1A_ToneStack.h"
#include "V1B_Stage.h"
#include "V3B_Stage.h"

class MarkIICPreamp {
private:
    // Instantiate all 9 processing stages
    // This is a conceptual layout. The actual segments would be implemented
    // following the `V1A_ToneStack` example.
    V1A_ToneStack stage1;
    V1B_Stage stage2;
    Coupling1 stage3;
    V3B_Stage stage4;
    LeadPathCoupling stage5;
    // V4A_Stage stage6;
    // LeadPathCoupling2 stage7;
    // MixerAndV2B stage8;
    // V2A_OutputStage stage9;

public:
    MarkIICPreamp() = default;

    void setup(double sampleRate) {
        stage1.setup(sampleRate);
        stage2.setup(sampleRate);
        stage3.setup(sampleRate);
        stage4.setup(sampleRate);
        stage5.setup(sampleRate);
    }

    // Parameter setting methods
    void setVolume1(double val) { stage1.setVolume1(val); }
    void setTreble(double val) { stage1.setTreble(val); }
    void setBass(double val) { stage1.setBass(val); }
    void setMid(double val) { stage1.setMid(val); }
    void setGain(double val) { stage3.setGain(val); }
    // void setMaster(double val) { stage9.setMaster(val); }

    double processSample(double in) {
        double out_s1 = stage1.process(in);
        // return out_s1;

        double out_s2 = stage2.process(out_s1);
        // return out_s2;

        double out_s3 = stage3.process(out_s2); // out is N001
        // return out_s3;

        double lead_path_in = out_s3;
        double out_s4 = stage4.process(lead_path_in);
        // return out_s4;

        double out_s5 = stage5.process(out_s4);
        return out_s5;
        // double out_s6 = stage6.process(out_s5);
        // double out_s7 = stage7.process(out_s6); // out is N025

        // double out_s8 = stage8.process(out_s3, out_s7); // Mixer stage
        // double out_s9 = stage9.process(out_s8);

        // The final E1 block is a simple gain stage
        // return out_s9 / 1400.0;
    }
};
