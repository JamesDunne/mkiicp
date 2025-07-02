#include "TubeStageV3B.h"
#include "TubeModel.h"

namespace V3B_Nodes { constexpr int Plate = 0, Cathode = 1; }

TubeStageV3B::TubeStageV3B(double sr) : MNASolver<2>(sr, true) {
    flagForRecalculation();
}

void TubeStageV3B::stampLinear(Matrix<2>& G, Matrix<2>& C) {
    using namespace V3B_Nodes;
    const double R26 = 82e3;  // Plate load
    const double R23 = 1.5e3; // Cathode resistor
    const double C23 = 2.2e-6; // Cathode bypass
    const double C22 = 120e-12; // Grid-to-cathode coupling cap

    G[Plate][Plate] += 1.0 / R26;
    G[Cathode][Cathode] += 1.0 / R23;
    C[Cathode][Cathode] += C23;

    // C22 connects Cathode to Grid (which is an input).
    // This adds a capacitive load on the cathode and also an input term.
    // In a full MNA, this would be: C22 * d/dt (Vc - Vg).
    // The Vc part is stamped here. The Vg part is handled in stampNonlinear
    // as it relates to the input.
    C[Cathode][Cathode] += C22;
}

void TubeStageV3B::stampNonlinear(Matrix<2>& J, Vector<2>& b, const Vector<2>& x, const std::vector<double>& u) {
    using namespace V3B_Nodes;
    // Call the shared tube model helper
    TubeModel::stamp(J, b, x, u[0], Plate, Cathode);

    // Add contribution from C22 (grid-cathode cap) related to the input grid voltage.
    // Current into cathode from C22 is C * d/dt(Vg - Vc).
    // The Vc part is in stampLinear. The Vg part is an input source I = C * d/dt(Vg).
    // Using Bilinear Transform I_n = C * g_s * (Vg_n - Vg_n-1).
    // We can approximate this by adding it to the RHS.
    // This is a simplification; a more rigorous model would pass u_n-1.
    // However, given it's part of the Newton-Raphson iteration, we add the derivative.
    double g_s = 2.0 * 48000.0; // Needs to be passed properly from base
    const double C22 = 120e-12;
    J[Cathode][Cathode] += g_s * C22;
    b[Cathode] += g_s * C22 * u[0];
}
