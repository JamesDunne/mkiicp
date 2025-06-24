// PreampSimulator.hpp
// Header-only vacuum tube preamp simulator for real-time audio
// Derived from LTspice netlist and TriodeK model of 12AX7
#pragma once

#include <array>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <vector>

namespace PreampSim {

constexpr int kNumNodes = 40; // Number of voltage nodes (adjust as needed)
constexpr int kNumVSources = 3; // Number of independent voltage sources
constexpr int kSystemSize = kNumNodes + kNumVSources;
constexpr int kMaxNRIterations = 20;
constexpr float kTolerance = 1e-6f;

struct Params {
    float vol1;   // [0.0, 1.0]
    float treble; // [0.0, 1.0]
    float bass;   // [0.0, 1.0]
    float mid;    // [0.0, 1.0]
    float gain;   // [0.0, 1.0]
};

struct Triode {
    int anode;
    int grid;
    int cathode;
};

class PreampSimulator {
public:
    PreampSimulator(float sampleRate, int numNodes, int numVSources)
        : sampleRate(sampleRate),
          numNodes(numNodes),
          numVSources(numVSources),
          kSystemSize(numNodes + numVSources),
          G(kSystemSize * kSystemSize, 0.0f),
          C(kSystemSize * kSystemSize, 0.0f),
          x(kSystemSize, 0.0f),
          b(kSystemSize, 0.0f) {
        dt = 1.0f / sampleRate;
        reset();
    }

    void setParams(const Params& p) {
        this->params = p;
        updateParamDependentElements();
    }

    float process(float inSample) {
        newtonRaphsonSolve();
        return x[0]; // Replace with actual output node
    }

    void addTriode(int anode, int grid, int cathode) {
        triodes.push_back({anode, grid, cathode});
    }

    void stampDCVoltageSource(int n_plus, int n_minus, int vIndex, float voltage) {
        int k = numNodes + vIndex;
        if (n_plus != 0) {
            G[n_plus * kSystemSize + k] += 1.0f;
            G[k * kSystemSize + n_plus] += 1.0f;
        }
        if (n_minus != 0) {
            G[n_minus * kSystemSize + k] -= 1.0f;
            G[k * kSystemSize + n_minus] -= 1.0f;
        }
        b[k] = voltage;
    }

    void stampResistor(int n1, int n2, float resistance) {
        float g = 1.0f / resistance;
        if (n1 != 0) G[n1 * kSystemSize + n1] += g;
        if (n2 != 0) G[n2 * kSystemSize + n2] += g;
        if (n1 != 0 && n2 != 0) {
            G[n1 * kSystemSize + n2] -= g;
            G[n2 * kSystemSize + n1] -= g;
        }
    }

    void stampCapacitor(int n1, int n2, float capacitance) {
        float c = capacitance / dt;
        if (n1 != 0) C[n1 * kSystemSize + n1] += c;
        if (n2 != 0) C[n2 * kSystemSize + n2] += c;
        if (n1 != 0 && n2 != 0) {
            C[n1 * kSystemSize + n2] -= c;
            C[n2 * kSystemSize + n1] -= c;
        }
    }

private:
    float sampleRate;
    float dt;
    int numNodes;
    int numVSources;
    int kSystemSize;
    Params params;

    std::vector<float> G;
    std::vector<float> C;
    std::vector<float> x;
    std::vector<float> b;
    std::vector<Triode> triodes;

    void reset() {
        std::fill(x.begin(), x.end(), 0.0f);
    }

    void updateParamDependentElements() {
        // Update elements that depend on vol1, treble, bass, mid, gain
    }

    void newtonRaphsonSolve() {
        std::vector<float> x_guess = x;

        for (int iter = 0; iter < kMaxNRIterations; ++iter) {
            std::vector<float> G_backup = G;
            std::vector<float> b_backup = b;

            for (const auto& triode : triodes) {
                float Vak = x_guess[triode.anode] - x_guess[triode.cathode];
                float Vgk = x_guess[triode.grid] - x_guess[triode.cathode];
                float I = triodeCurrent(Vak, Vgk);

                float dVak = 1e-3f;
                float dVgk = 1e-3f;
                float dIa_dVak = (triodeCurrent(Vak + dVak, Vgk) - I) / dVak;
                float dIa_dVgk = (triodeCurrent(Vak, Vgk + dVgk) - I) / dVgk;

                if (triode.anode != 0) G[triode.anode * kSystemSize + triode.anode] += dIa_dVak;
                if (triode.anode != 0 && triode.cathode != 0) {
                    G[triode.anode * kSystemSize + triode.cathode] -= dIa_dVak + dIa_dVgk;
                    G[triode.cathode * kSystemSize + triode.anode] -= dIa_dVak;
                    G[triode.cathode * kSystemSize + triode.cathode] += dIa_dVak + dIa_dVgk;
                }
                if (triode.grid != 0 && triode.cathode != 0) {
                    G[triode.grid * kSystemSize + triode.cathode] -= dIa_dVgk;
                    G[triode.cathode * kSystemSize + triode.grid] -= dIa_dVgk;
                    G[triode.grid * kSystemSize + triode.grid] += dIa_dVgk;
                }

                if (triode.anode != 0) b[triode.anode] -= I;
                if (triode.cathode != 0) b[triode.cathode] += I;
            }

            solve();

            float maxDelta = 0.0f;
            for (int i = 0; i < kSystemSize; ++i) {
                float delta = std::fabs(x[i] - x_guess[i]);
                maxDelta = std::max(maxDelta, delta);
                x_guess[i] = x[i];
            }

            if (maxDelta < kTolerance) break;

            G = G_backup;
            b = b_backup;
        }
    }

    void solve() {
        std::vector<float> A(kSystemSize * kSystemSize);
        for (int i = 0; i < kSystemSize * kSystemSize; ++i) A[i] = G[i] + C[i];

        std::vector<float> y = b;
        std::vector<float> x_out(kSystemSize);

        for (int k = 0; k < kSystemSize; ++k) {
            float pivot = A[k * kSystemSize + k];
            assert(std::fabs(pivot) > 1e-12f);
            for (int i = k + 1; i < kSystemSize; ++i) {
                float factor = A[i * kSystemSize + k] / pivot;
                for (int j = k; j < kSystemSize; ++j)
                    A[i * kSystemSize + j] -= factor * A[k * kSystemSize + j];
                y[i] -= factor * y[k];
            }
        }

        for (int i = kSystemSize - 1; i >= 0; --i) {
            float sum = y[i];
            for (int j = i + 1; j < kSystemSize; ++j)
                sum -= A[i * kSystemSize + j] * x_out[j];
            x_out[i] = sum / A[i * kSystemSize + i];
        }

        x = x_out;
    }

    [[nodiscard]] static float triodeCurrent(float Vak, float Vgk) {
        constexpr float MU = 96.2f;
        constexpr float EX = 1.437f;
        constexpr float KG1 = 613.4f;
        constexpr float KP = 740.3f;
        constexpr float KVB = 1672.0f;

        float E1 = Vak / KP * std::log1p(std::exp(KP * (1.0f / MU + Vgk / std::sqrt(KVB + Vak * Vak))));
        float I = 0.5f * (std::pow(E1, EX) + std::pow(std::fabs(E1), EX)) / KG1;
        return I;
    }
};

} // namespace PreampSim
