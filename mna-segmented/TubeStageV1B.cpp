#include "TubeStageV1B.h"
#include "TubeModel.h" // The shared helper for non-linear tube calculations

// Use the namespace for node indices for improved readability.
using namespace V1B_Nodes;

TubeStageV1B::TubeStageV1B(double sampleRate)
    : MNASolver<2>(sampleRate, true) // true = this is a non-linear circuit
{
    // Flag for the initial calculation of the system matrices.
    flagForRecalculation();
}

void TubeStageV1B::stampLinear(Matrix<2>& G, Matrix<2>& C) {
    // Helper lambda for converting resistance to conductance.
    auto g = [](double r) { return (r > 1e-12) ? 1.0 / r : 1e12; };

    // --- Define Component Values for the V1B Stage ---
    const double R8 = 100e3;  // Plate load resistor to B+ (AC ground)
    const double R7 = 1.5e3;   // Cathode resistor to ground
    const double C13 = 22e-6; // Cathode bypass capacitor to ground

    // Stamp the components into the G (conductance) and C (capacitance) matrices.

    // R8 is connected from the Plate (node 0) to AC ground.
    G[Plate][Plate] += g(R8);

    // R7 is connected from the Cathode (node 1) to ground.
    G[Cathode][Cathode] += g(R7);

    // C13 is connected from the Cathode (node 1) to ground.
    C[Cathode][Cathode] += C13;
}

void TubeStageV1B::stampNonlinear(Matrix<2>& J, Vector<2>& b, const Vector<2>& x, const std::vector<double>& u) {
    // This function adds the contribution of the non-linear V1B tube.
    // The grid voltage is the external input to this block.
    double v_grid = u[0];

    // Delegate the complex calculations to the TubeModel helper function.
    // We provide the Jacobian (J), RHS vector (b), previous state (x),
    // current grid voltage (v_grid), and the node indices for this stage's plate and cathode.
    TubeModel::stamp(J, b, x, v_grid, Plate, Cathode);
}
