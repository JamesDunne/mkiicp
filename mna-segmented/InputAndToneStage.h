#pragma once
#include "MNASolver.h"

// The number of unknown nodes for this combined circuit segment.
constexpr int V1A_TONESTACK_NODES = 9;

// Namespace to provide clear, readable names for the node indices within the implementation.
namespace V1A_Nodes {
    constexpr int N004 = 0; // Plate of V1A, input to the tone stack
    constexpr int N033 = 1; // Cathode of V1A
    constexpr int N020 = 2; // Output of the segment (wiper of vol1)
    constexpr int N015 = 3; // Tone stack internal node
    constexpr int N005 = 4; // Tone stack internal node
    constexpr int N007 = 5; // Tone stack internal node (treble pot wiper)
    constexpr int N008 = 6; // Tone stack internal node (volume pot input)
    constexpr int N016 = 7; // Tone stack internal node (bass/treble junction)
    constexpr int N028 = 8; // Tone stack internal node (bass/mid junction)
}

/**
 * @class InputAndToneStage
 * @brief Simulates the first tube stage (V1A) and the entire tone stack as a single, interacting, non-linear system.
 *
 * This class inherits from MNASolver and implements the virtual "stamping" functions
 * to describe the circuit's topology. It provides a public interface to set the
 * tone parameters and process audio.
 */
class InputAndToneStage : public MNASolver<V1A_TONESTACK_NODES> {
public:
    /**
     * @brief Constructs the InputAndToneStage.
     * @param sampleRate The audio sample rate.
     */
    InputAndToneStage(double sampleRate);

    /**
     * @brief Processes a single audio sample.
     * @param inputSample The input voltage at the grid of V1A.
     * @return The output voltage at the wiper of the volume pot (N020).
     */
    double process(double inputSample) {
        return MNASolver::process(inputSample);
    }

    /**
     * @brief Sets the tone and volume parameters.
     * @param treble Control value from 0.0 to 1.0.
     * @param mid Control value from 0.0 to 1.0.
     * @param bass Control value from 0.0 to 1.0.
     * @param vol1 Control value from 0.0 to 1.0.
     */
    void setParams(double treble, double mid, double bass, double vol1);

protected:
    // --- Implementation of MNASolver's virtual functions ---

    /**
     * @brief Stamps all linear components (resistors, capacitors) onto the G and C matrices.
     * This is called by the solver only when parameters have changed.
     */
    void stampLinear(Matrix<V1A_TONESTACK_NODES>& G, Matrix<V1A_TONESTACK_NODES>& C) override;

    /**
     * @brief Stamps the non-linear tube model (V1A) onto the Jacobian and RHS vector.
     * This is called by the solver every single audio sample.
     */
    void stampNonlinear(Matrix<V1A_TONESTACK_NODES>& J, Vector<V1A_TONESTACK_NODES>& b, const Vector<V1A_TONESTACK_NODES>& x_k, const std::vector<double>& u_n) override;

    /**
     * @brief Not used for non-linear circuits where the input is handled in stampNonlinear.
     * Must be defined as it's a pure virtual function in the base class.
     */
    void stampInputs(std::vector<std::pair<Vector<V1A_TONESTACK_NODES>, Vector<V1A_TONESTACK_NODES>>>& i) override {}

    /**
     * @brief Returns the index of the node that serves as the output of this circuit block.
     * @return The integer index for node N020.
     */
    int getOutputNodeIndex() const override {
        return V1A_Nodes::N020;
    }

private:
    // Member variables to hold the current state of the user-controlled parameters.
    double p_treble, p_mid, p_bass, p_vol1;
};
