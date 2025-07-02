#pragma once
#include "MNASolver.h"

// Namespace for clear, readable node indices within the implementation.
namespace V1B_Nodes {
    constexpr int Plate = 0;   // Corresponds to circuit node N009
    constexpr int Cathode = 1; // Corresponds to circuit node N027
}

/**
 * @class TubeStageV1B
 * @brief Simulates the V1B common-cathode amplifier stage.
 *
 * This stage provides gain recovery after the signal loss from the passive tone stack.
 * It's a non-linear stage whose behavior is modeled with a per-sample Newton-Raphson solve.
 */
class TubeStageV1B : public MNASolver<2> {
public:
    /**
     * @brief Constructs the TubeStageV1B.
     * @param sampleRate The audio sample rate.
     */
    TubeStageV1B(double sampleRate);

    /**
     * @brief Processes a single audio sample.
     * @param v_grid_N020 The input voltage at the grid (from the tone stack output).
     * @return The output voltage at the plate (N009).
     */
    double process(double v_grid_N020) {
        return MNASolver::process(v_grid_N020);
    }

protected:
    // --- Implementation of MNASolver's virtual functions ---

    /**
     * @brief Stamps the linear resistors and capacitors for this stage.
     */
    void stampLinear(Matrix<2>& G, Matrix<2>& C) override;

    /**
     * @brief Stamps the non-linear tube model.
     */
    void stampNonlinear(Matrix<2>& J, Vector<2>& b, const Vector<2>& x, const std::vector<double>& u) override;

    /**
     * @brief Not used for this non-linear circuit.
     */
    void stampInputs(std::vector<std::pair<Vector<2>, Vector<2>>>& i) override {}

    /**
     * @brief Returns the index of the output node.
     * @return The integer index for the Plate node (0).
     */
    int getOutputNodeIndex() const override {
        return V1B_Nodes::Plate;
    }
};
