#include "TubeStageV1B.h"
#include "TubeParams.h"
#include <cmath>

namespace V1B_Nodes { constexpr int Plate = 0, Cathode = 1; }

TubeStageV1B::TubeStageV1B(double sampleRate) : MNASolver<2>(sampleRate, true) {
    flagForRecalculation();
}

void TubeStageV1B::stampLinear(Matrix<2>& G, Matrix<2>& C) {
    using namespace V1B_Nodes;
    const double R8 = 100e3;  // Plate load
    const double R7 = 1.5e3;   // Cathode resistor
    const double C13 = 22e-6; // Cathode bypass

    G[Plate][Plate] += 1.0 / R8;
    G[Cathode][Cathode] += 1.0 / R7;
    C[Cathode][Cathode] += C13;
}

void TubeStageV1B::stampNonlinear(Matrix<2>& J, Vector<2>& b, const Vector<2>& x, const std::vector<double>& u) {
    using namespace V1B_Nodes;
    double v_p = x[Plate], v_c = x[Cathode], v_g = u[0];
    
    // Koren tube model evaluation (this logic will be repeated in all tube stages)
    double v_gp = v_g - v_p;
    double v_gc = v_g - v_c;
    double e1_arg = Tube::KP * (1.0/Tube::MU + v_gc / sqrt(Tube::KVB + v_gp*v_gp));
    double e1 = log1p(exp(e1_arg)) / Tube::KP;
    double g1 = pow(e1, Tube::EX) / Tube::KG1;
    
    // Partial derivatives (simplified for brevity, a full implementation is needed)
    double dG1_dVp = 0.0, dG1_dVc = 0.0, dG1_dVg = 0.0;
    // ... calculate actual derivatives here ...

    // Stamp current source and its derivatives
    b[Plate] -= g1; b[Cathode] += g1;
    b[Plate] -= dG1_dVg * v_g; b[Cathode] += dG1_dVg * v_g;
    
    J[Plate][Plate]   += -dG1_dVp; J[Plate][Cathode]   += -dG1_dVc;
    J[Cathode][Plate] += dG1_dVp;  J[Cathode][Cathode] += dG1_dVc;
}
