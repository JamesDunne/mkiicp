#pragma once
#include "MNASolver.h" // For Matrix/Vector types
#include "TubeParams.h"
#include <cmath>

namespace TubeModel {
    // Helper to evaluate the Koren tube model and its partial derivatives
    // It will stamp the results directly onto the Jacobian and RHS vector
    template<int N>
    void stamp(Matrix<N>& J, Vector<N>& b,
               const Vector<N>& x_k, double v_g_k,
               int plate_idx, int cathode_idx) {
        
        double v_p_k = x_k[plate_idx];
        double v_c_k = x_k[cathode_idx];

        // --- Evaluate Model at x_k ---
        double v_gp = v_g_k - v_p_k;
        double v_gc = v_g_k - v_c_k;

        if (v_gc < -50.0) { // Tube is in cutoff, effectively an open circuit
             return; // No current, no derivatives to stamp
        }

        double v_d_sq = Tube::KVB + v_gp * v_gp;
        if (v_d_sq < 1e-6) v_d_sq = 1e-6; // Avoid division by zero
        double v_d = sqrt(v_d_sq);

        double e1_arg = Tube::KP * (1.0 / Tube::MU + v_gc / v_d);
        
        // Use log1p for numerical stability when argument is near zero
        double e1 = (e1_arg > 30.0) ? e1_arg / Tube::KP : log1p(exp(e1_arg)) / Tube::KP;

        // Plate current (using a safe power function)
        auto safe_pow = [](double base, double exp) { return (base > 0) ? pow(base, exp) : 0.0; };
        double g1 = safe_pow(e1, Tube::EX) / Tube::KG1;

        // --- Calculate Partial Derivatives via Chain Rule ---
        double exp_e1_arg = exp(e1_arg);
        double dE1_dArg = (e1_arg > 30.0) ? 1.0/Tube::KP : exp_e1_arg / (1.0 + exp_e1_arg) / Tube::KP;

        double dArg_dVgc = Tube::KP / v_d;
        double dArg_dVgp = -Tube::KP * v_gc * v_gp / (v_d * v_d_sq);

        double dE1_dVgc = dE1_dArg * dArg_dVgc;
        double dE1_dVgp = dE1_dArg * dArg_dVgp;

        double dG1_dE1 = Tube::EX * safe_pow(e1, Tube::EX - 1.0) / Tube::KG1;

        double dG1_dVgc = dG1_dE1 * dE1_dVgc;
        double dG1_dVgp = dG1_dE1 * dE1_dVgp;

        // --- Stamp values onto J and b ---
        // Current Source: I(p,c) = g1
        b[plate_idx] -= g1;
        b[cathode_idx] += g1;

        // Derivatives of the current source
        double dG1_dVp = -dG1_dVgp;
        double dG1_dVc = -dG1_dVgc;
        double dG1_dVg = dG1_dVgp + dG1_dVgc;

        J[plate_idx][plate_idx] += -dG1_dVp;
        J[plate_idx][cathode_idx] += -dG1_dVc;
        J[cathode_idx][plate_idx] += dG1_dVp;
        J[cathode_idx][cathode_idx] += dG1_dVc;

        // The grid voltage is an independent variable (input), so its
        // contribution goes to the RHS vector b.
        // Equivalent to moving J_ig * u_n to the other side: -J_ig * u_n
        b[plate_idx] -= dG1_dVg * v_g_k;
        b[cathode_idx] += dG1_dVg * v_g_k;
    }
}
