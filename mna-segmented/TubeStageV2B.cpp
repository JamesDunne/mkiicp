#include "TubeStageV2B.h"
#include "TubeModel.h"

namespace V2B_Nodes { constexpr int Plate = 0, Cathode = 1; }

TubeStageV2B::TubeStageV2B(double sr) : MNASolver<2>(sr, true) {
    flagForRecalculation();
}

void TubeStageV2B::stampLinear(Matrix<2>& G, Matrix<2>& C) {
    using namespace V2B_Nodes;
    const double R13 = 100e3; // Plate load
    const double R16 = 1.5e3; // Cathode resistor
    G[Plate][Plate] += 1.0 / R13;
    G[Cathode][Cathode] += 1.0 / R16;
}

void TubeStageV2B::stampNonlinear(Matrix<2>& J, Vector<2>& b, const Vector<2>& x, const std::vector<double>& u) {
    using namespace V2B_Nodes;
    TubeModel::stamp(J, b, x, u[0], Plate, Cathode);
}
