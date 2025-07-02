#pragma once
#include "MNASolver.h"

constexpr int V1A_TONESTACK_NODES = 10;

// Node mapping for clarity inside the implementation
namespace V1A_Nodes {
    constexpr int N004 = 0, N033 = 1, E1_V1A_UNUSED = 2, N020 = 3, N015 = 4, 
                  N005 = 5, N007 = 6, N008 = 7, N016 = 8, N028 = 9;
}

class InputAndToneStage : public MNASolver<V1A_TONESTACK_NODES> {
public:
    InputAndToneStage(double sampleRate);

    // Provide a public, user-friendly process method
    double processStage(double inputSample) {
        return process(inputSample);
    }
    
    void setParams(double treble, double mid, double bass, double vol1);

protected:
    // Implement the virtual functions required by MNASolver
    void stampLinear(Matrix<V1A_TONESTACK_NODES>& G, Matrix<V1A_TONESTACK_NODES>& C) override;
    
    void stampInputs(std::vector<std::pair<Vector<V1A_TONESTACK_NODES>, Vector<V1A_TONESTACK_NODES>>>& inputs) override {
        // This is a non-linear circuit, so this function is not used by the solver.
        // The input is handled directly in stampNonlinear.
    }

    void stampNonlinear(Matrix<V1A_TONESTACK_NODES>& J, Vector<V1A_TONESTACK_NODES>& b, const Vector<V1A_TONESTACK_NODES>& x_k, const std::vector<double>& u_n) override;

    int getOutputNodeIndex() const override {
        return V1A_Nodes::N020;
    }

private:
    // User-controlled parameters
    double p_treble, p_mid, p_bass, p_vol1;

    // Koren Tube Model Parameters for 12AX7
    const double MU = 96.20, EX = 1.437, KG1 = 613.4, KP = 740.3, KVB = 1672.0;
};
