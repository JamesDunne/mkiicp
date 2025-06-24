#pragma once

#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <memory>

class MKIICSimulator {
private:
    // Circuit parameters
    struct Parameters {
        double vol1 = 0.75;
        double treble = 0.8;
        double mid = 0.33;
        double bass = 0.05;
        double gain = 0.5;
        double master = 0.5;
    } params;

    // Node mapping from SPICE netlist
    std::map<std::string, int> nodeMap;
    int numNodes, numVSources;
    int kSystemSize;

    // MNA matrices
    std::vector<std::vector<double>> G; // Conductance matrix
    std::vector<std::vector<double>> C; // Capacitance matrix
    std::vector<double> b;              // Current source vector
    std::vector<double> x;              // Solution vector (node voltages)
    std::vector<double> x_prev;         // Previous solution for integration

    // Track variable resistor stamps for removal
    struct ResistorStamp {
        int node1, node2;
        double conductance;
    };
    std::map<std::string, ResistorStamp> variableResistors;

    // Tube model parameters (12AX7)
    struct TubeParams {
        double MU = 96.20;
        double EX = 1.437;
        double KG1 = 613.4;
        double KP = 740.3;
        double KVB = 1672.0;
        double RGI = 2000.0;
    } tubeParams;

    // Tube instances
    struct Triode {
        int anode, grid, cathode;
    };

    std::vector<Triode> tubes;

    // Sample rate and time step
    double sampleRate = 48000.0;
    double dt;

    // LU decomposition storage
    std::vector<std::vector<double>> LU;
    std::vector<int> pivot;

public:
    MKIICSimulator(double fs = 48000.0) : sampleRate(fs), dt(1.0/fs) {
        initializeCircuit();
    }

    // Parameter control methods
    void setVolume1(double vol) { params.vol1 = std::clamp(vol, 0.0, 1.0); updateToneStack(); }
    void setTreble(double treble) { params.treble = std::clamp(treble, 0.0, 1.0); updateToneStack(); }
    void setMid(double mid) { params.mid = std::clamp(mid, 0.0, 1.0); updateToneStack(); }
    void setBass(double bass) { params.bass = std::clamp(bass, 0.0, 1.0); updateToneStack(); }
    void setGain(double gain) { params.gain = std::clamp(gain, 0.0, 1.0); updateGainControl(); }
    void setMaster(double master) { params.master = std::clamp(master, 0.0, 1.0); updateMasterControl(); }

    // Main processing function
    double processSample(double input) {
        // Set input voltage source
        int inputNode = nodeMap["Vin"];
        b[inputNode] = input;

        newtonRaphsonSolve();

        // Store previous solution for next time step
        x_prev = x;

        // Extract output (scaled by 1/1400 as in SPICE)
        int outputNode = nodeMap["N014"];
        return x[outputNode] / 1400.0;
    }

private:
    void initializeCircuit() {
        // Initialize node mapping (0 = ground)
        nodeMap["0"] = 0;

        // Add all nodes from the netlist
        std::vector<std::string> nodes = {
            "N001", "N002", "N003", "N004", "N005", "N006", "N007", "N008", "N009", "N010",
            "N011", "N012", "N013", "N014", "N015", "N016", "N017", "N018", "N019", "N020",
            "N021", "N022", "N023", "N024", "N025", "N026", "N027", "N028", "N029", "N030",
            "N031", "N032", "N033", "N034", "N035", "N036", "P001"
        };
        std::vector<std::string> vsources = {
            "Vin", "VE", "VC", "VC2"
        };

        int i;
        for (i = 0; i < nodes.size(); i++) {
            nodeMap[nodes[i]] = i + 1;
        }
        numNodes = nodes.size() + 1; // +1 for ground

        for (size_t j = 0; j < vsources.size(); j++, i++) {
            nodeMap[vsources[j]] = i + 1;
        }
        numVSources = vsources.size();

        kSystemSize = numNodes + numVSources;

        // Initialize matrices
        G.assign(kSystemSize, std::vector<double>(kSystemSize, 0.0));
        C.assign(kSystemSize, std::vector<double>(kSystemSize, 0.0));
        b.assign(kSystemSize, 0.0);
        x.assign(kSystemSize, 0.0);
        x_prev.assign(kSystemSize, 0.0);

        // Initialize LU decomposition storage
        LU.assign(kSystemSize, std::vector<double>(kSystemSize, 0.0));
        pivot.assign(kSystemSize, 0);

        // Stamp passive components
        stampPassiveComponents();

        // Initialize tone stack and gain controls
        updateToneStack();
        updateGainControl();
        updateMasterControl();

        // Initialize tube stages
        initializeTubeStages();

        // Perform initial LU decomposition
        updateSystemMatrix();
    }

    void stampPassiveComponents() {
        // Resistors - stamp conductances
        stampResistor("N004", "N015", 100e3);     // R5
        stampResistor("N008", "N007", 100e3);     // R5A
        stampResistor("N003", "N004", 150e3);     // R4
        stampResistor("N019", "0", 1e6);          // R1
        stampResistor("N033", "0", 1.5e3);        // R2
        stampResistor("N028", "0", 1.5e3);        // R7
        stampResistor("N003", "N009", 100e3);     // R8
        stampResistor("N001", "0", 100e3);        // R9
        stampResistor("N011", "N025", 680e3);     // R21
        stampResistor("N029", "0", 475e3);        // R22
        stampResistor("N035", "0", 1.5e3);        // R23
        stampResistor("N002", "N001", 3.3e6);     // R10
        stampResistor("N002", "0", 680e3);        // R11
        stampResistor("N017", "N010", 82e3);      // R26
        stampResistor("N030", "0", 68e3);         // R24
        stampResistor("N030", "N018", 270e3);     // R25
        stampResistor("N034", "0", 3.3e3);        // R30
        stampResistor("N026", "N010", 274e3);     // R27
        stampResistor("N002", "N027", 220e3);     // R31
        stampResistor("N036", "0", 1.5e3);        // R16
        stampResistor("N006", "N021", 100e3);     // R13
        stampResistor("N022", "P001", 47e3);      // R105
        stampResistor("N022", "0", 47e3);         // R46
        stampResistor("N022", "N032", 150e3);     // R102
        stampResistor("N032", "0", 4.7e3);        // R101
        stampResistor("N006", "N012", 120e3);     // R19
        stampResistor("N023", "0", 47e3);         // R103
        stampResistor("N031", "0", 1e3);          // R104
        stampResistor("N023", "N032", 2.2e3);     // R12
        stampResistor("N014", "N013", 15e3);      // R106
        stampResistor("N005", "N005", 10e6);      // R6
        stampResistor("N002", "0", 100e3);        // R32

        // Capacitors
        stampCapacitor("N005", "N004", 750e-12);  // C6
        stampCapacitor("N005", "N004", 250e-12);  // C5 (parallel with C6)
        stampCapacitor("N016", "N015", 0.1e-6);   // C4
        stampCapacitor("N024", "N015", 0.047e-6); // C3
        stampCapacitor("N020", "N008", 180e-12);  // C13B
        stampCapacitor("N033", "0", 0.47e-6);     // C1
        stampCapacitor("N033", "0", 22e-6);       // C2
        stampCapacitor("N028", "0", 22e-6);       // C13
        stampCapacitor("N001", "N009", 0.1e-6);   // C7
        stampCapacitor("N035", "N029", 120e-12);  // C22
        stampCapacitor("N035", "0", 2.2e-6);      // C23
        stampCapacitor("N002", "N001", 20e-12);   // C10
        stampCapacitor("N030", "0", 1000e-12);    // C24
        stampCapacitor("N034", "0", 0.22e-6);     // C29
        stampCapacitor("N027", "N026", 0.047e-6); // C30
        stampCapacitor("N002", "N027", 250e-12);  // C31
        stampCapacitor("P001", "N021", 0.047e-6); // C9
        stampCapacitor("N013", "N012", 0.047e-6); // C12
        stampCapacitor("N002", "0", 47e-12);      // C11
        stampCapacitor("N031", "0", 0.47e-6);     // C16
        stampCapacitor("N031", "0", 15e-6);       // C15
        stampCapacitor("N001", "N011", 0.02e-6);  // C21
        stampCapacitor("N002", "0", 500e-12);     // C32
        stampCapacitor("N018", "N017", 0.022e-6); // C25

        // input:
        stampVoltageSource("Vin", "N019", "0", 1.0); // Vin

        // Voltage sources (DC bias)
        stampVoltageSource("VE", "N003", "0", 405);      // VE
        stampVoltageSource("VC", "N010", "0", 410);      // VC
        stampVoltageSource("VC2", "N006", "0", 410);     // VC2
    }

    void stampResistor(const std::string& n1, const std::string& n2, double R) {
        int node1 = nodeMap[n1];
        int node2 = nodeMap[n2];
        double g = 1.0 / R;

        if (node1 != 0) {
            G[node1][node1] += g;
            if (node2 != 0) G[node1][node2] -= g;
        }
        if (node2 != 0) {
            G[node2][node2] += g;
            if (node1 != 0) G[node2][node1] -= g;
        }
    }

    void stampCapacitor(const std::string& n1, const std::string& n2, double Cap) {
        int node1 = nodeMap[n1];
        int node2 = nodeMap[n2];

        if (node1 != 0) {
            C[node1][node1] += Cap;
            if (node2 != 0) C[node1][node2] -= Cap;
        }
        if (node2 != 0) {
            C[node2][node2] += Cap;
            if (node1 != 0) C[node2][node1] -= Cap;
        }
    }

    void stampVariableResistor(const std::string& name, const std::string& n1, const std::string& n2, double R) {
        // Remove previous stamp if it exists
        removeVariableResistor(name);

        // Stamp new resistor
        int node1 = nodeMap[n1];
        int node2 = nodeMap[n2];
        double g = (R > 1e-12) ? 1.0 / R : 0.0; // Avoid division by zero

        if (node1 != 0) {
            G[node1][node1] += g;
            if (node2 != 0) G[node1][node2] -= g;
        }
        if (node2 != 0) {
            G[node2][node2] += g;
            if (node1 != 0) G[node2][node1] -= g;
        }

        // Store the stamp for later removal
        variableResistors[name] = {node1, node2, g};
    }

    void removeVariableResistor(const std::string& name) {
        auto it = variableResistors.find(name);
        if (it != variableResistors.end()) {
            const ResistorStamp& stamp = it->second;
            int node1 = stamp.node1;
            int node2 = stamp.node2;
            double g = stamp.conductance;

            // Remove the stamp by subtracting it
            if (node1 != 0) {
                G[node1][node1] -= g;
                if (node2 != 0) G[node1][node2] += g;
            }
            if (node2 != 0) {
                G[node2][node2] -= g;
                if (node1 != 0) G[node2][node1] += g;
            }

            variableResistors.erase(it);
        }
    }

    void stampVoltageSource(const std::string& name, const std::string& nPos, const std::string& nNeg, double voltage) {
        int k = nodeMap[name];
        int n_plus = nodeMap[nPos];
        if (n_plus != 0) {
            G[n_plus][k] += 1.0;
            G[k][n_plus] += 1.0;
        }
        int n_minus = nodeMap[nNeg];
        if (n_minus != 0) {
            G[n_minus][k] -= 1.0;
            G[k][n_minus] -= 1.0;
        }
        b[k] = voltage;
    }

    void updateToneStack() {
        // Volume control - potentiometer between N020 and N008
        double vol1 = params.vol1;
        if (vol1 < 1e-6) vol1 = 1e-6; // Prevent zero resistance
        if (vol1 > 1.0 - 1e-6) vol1 = 1.0 - 1e-6;

        stampVariableResistor("RA_VOLUME1", "N020", "N008", 1e6 * (1 - vol1));
        stampVariableResistor("RC_VOLUME1", "0", "N020", 1e6 * vol1);

        // Treble control - potentiometer between N005 and N016 via N007
        double treble = params.treble;
        if (treble < 1e-6) treble = 1e-6;
        if (treble > 1.0 - 1e-6) treble = 1.0 - 1e-6;

        stampVariableResistor("RA_TREBLE", "N005", "N007", 250e3 * (1 - treble));
        stampVariableResistor("RC_TREBLE", "N007", "N016", 250e3 * treble);

        // Bass control - potentiometer between N016 and N024
        double bass = params.bass;
        if (bass < 1e-6) bass = 1e-6;
        if (bass > 1.0 - 1e-6) bass = 1.0 - 1e-6;

        stampVariableResistor("RA_BASS", "N016", "N024", 250e3 * (1 - bass));
        stampVariableResistor("RC_BASS", "N024", "0", 250e3 * bass); // Fixed: was N024 to N024

        // Mid control - potentiometer from N024 to ground
        double mid = params.mid;
        if (mid < 1e-6) mid = 1e-6;
        if (mid > 1.0 - 1e-6) mid = 1.0 - 1e-6;

        stampVariableResistor("RA_MID", "N024", "0", 10e3 * (1 - mid));
        // Note: RC_MID connects 0 to 0 which is meaningless, so we skip it
    }

    void updateGainControl() {
        double gain = params.gain;
        if (gain < 1e-6) gain = 1e-6;
        if (gain > 1.0 - 1e-6) gain = 1.0 - 1e-6;

        stampVariableResistor("RA_GAIN", "N025", "N029", 1e6 * (1 - gain));
        stampVariableResistor("RC_GAIN", "N029", "0", 1e6 * gain);
    }

    void updateMasterControl() {
        double master = params.master;
        if (master < 1e-6) master = 1e-6;
        if (master > 1.0 - 1e-6) master = 1.0 - 1e-6;

        stampVariableResistor("RA_MASTER", "N014", "0", 1e6 * (1 - master));
        // Note: RC_MASTER connects 0 to 0 which is meaningless, so we skip it
    }

    void initializeTubeStages() {
        // Add tube stages based on SPICE netlist
        tubes.resize(6);
        
        // XV1A: 12AX7 N004 N019 N033
        tubes[0] = {nodeMap["N004"], nodeMap["N019"], nodeMap["N033"]};
        
        // XV1B: 12AX7 N009 N020 N028  
        tubes[1] = {nodeMap["N009"], nodeMap["N020"], nodeMap["N028"]};
        
        // XV3B: 12AX7 N017 N029 N035
        tubes[2] = {nodeMap["N017"], nodeMap["N029"], nodeMap["N035"]};
        
        // XV4A: 12AX7 N026 N030 N034
        tubes[3] = {nodeMap["N026"], nodeMap["N030"], nodeMap["N034"]};
        
        // XV2B: 12AX7 N021 N002 N036
        tubes[4] = {nodeMap["N021"], nodeMap["N002"], nodeMap["N036"]};
        
        // XV2A: 12AX7 N012 N023 N031
        tubes[5] = {nodeMap["N012"], nodeMap["N023"], nodeMap["N031"]};
    }
    
    void updateCapacitorStamps() {
        // Update the system matrix to include capacitor contributions for current timestep
        updateSystemMatrix();
    }
    
    void updateSystemMatrix() {
        // Combine G and C matrices with backward Euler integration: (G + C/dt) * x = b + C/dt * x_prev
        for (int i = 0; i < kSystemSize; i++) {
            for (int j = 0; j < kSystemSize; j++) {
                LU[i][j] = G[i][j] + C[i][j] / dt;
            }
        }
        
        // Perform LU decomposition
        luDecompose();
    }
    
    void luDecompose() {
        for (int i = 0; i < kSystemSize; i++) {
            pivot[i] = i;
        }
        
        for (int k = 0; k < kSystemSize - 1; k++) {
            // Find pivot
            int maxRow = k;
            for (int i = k + 1; i < kSystemSize; i++) {
                if (std::abs(LU[i][k]) > std::abs(LU[maxRow][k])) {
                    maxRow = i;
                }
            }
            
            // Swap rows
            if (maxRow != k) {
                std::swap(LU[k], LU[maxRow]);
                std::swap(pivot[k], pivot[maxRow]);
            }
            
            // Elimination
            for (int i = k + 1; i < kSystemSize; i++) {
                if (std::abs(LU[k][k]) > 1e-12) {
                    LU[i][k] /= LU[k][k];
                    for (int j = k + 1; j < kSystemSize; j++) {
                        LU[i][j] -= LU[i][k] * LU[k][j];
                    }
                }
            }
        }
    }
    
    void solveLU() {
        // Add capacitor current contribution to RHS
        std::vector<double> rhs = b;
        for (int i = 0; i < kSystemSize; i++) {
            for (int j = 0; j < kSystemSize; j++) {
                rhs[i] += C[i][j] * x_prev[j] / dt;
            }
        }
        
        // Forward substitution
        for (int i = 0; i < kSystemSize; i++) {
            x[i] = rhs[pivot[i]];
            for (int j = 0; j < i; j++) {
                x[i] -= LU[i][j] * x[j];
            }
        }
        
        // Back substitution
        for (int i = kSystemSize - 1; i >= 0; i--) {
            for (int j = i + 1; j < kSystemSize; j++) {
                x[i] -= LU[i][j] * x[j];
            }
            if (std::abs(LU[i][i]) > 1e-12) {
                x[i] /= LU[i][i];
            }
        }
    }

    void solve() {
        std::vector<std::vector<double>> A;
        A.assign(kSystemSize, std::vector<double>(kSystemSize, 0.0));

        for (int i = 0; i < kSystemSize; ++i)
            for (int j = 0; j < kSystemSize; ++j)
                A[i][j] = G[i][j] + C[i][j] / dt;

        std::vector<double> y = b;
        std::vector<double> x_out(kSystemSize);

        for (int k = 0; k < kSystemSize; ++k) {
            double pivot = A[k][k];
            if (std::fabs(pivot) <= 1e-12) continue;
            for (int i = k + 1; i < kSystemSize; ++i) {
                double factor = A[i][k] / pivot;
                for (int j = k; j < kSystemSize; ++j)
                    A[i][j] -= factor * A[k][j];
                y[i] -= factor * y[k];
            }
        }

        for (int i = kSystemSize - 1; i >= 0; --i) {
            double sum = y[i];
            for (int j = i + 1; j < kSystemSize; ++j)
                sum -= A[i][j] * x_out[j];
            x_out[i] = sum / A[i][i];
        }

        x = x_out;
    }

    void newtonRaphsonSolve() {
        constexpr int kMaxNRIterations = 20;
        constexpr double kTolerance = 1e-6;

        std::vector<double> x_guess = x;

        for (int iter = 0; iter < kMaxNRIterations; ++iter) {
            std::vector<std::vector<double>> G_backup = G;
            std::vector<double> b_backup = b;

            for (const auto& triode : tubes) {
                double Vak = x_guess[triode.anode] - x_guess[triode.cathode];
                double Vgk = x_guess[triode.grid] - x_guess[triode.cathode];
                double I = triodeCurrent(Vak, Vgk);

                double dVak = 1e-3;
                double dVgk = 1e-3;
                double dIa_dVak = (triodeCurrent(Vak + dVak, Vgk) - I) / dVak;
                double dIa_dVgk = (triodeCurrent(Vak, Vgk + dVgk) - I) / dVgk;

                if (triode.anode != 0) G[triode.anode][triode.anode] += dIa_dVak;
                if (triode.anode != 0 && triode.cathode != 0) {
                    G[triode.anode][triode.cathode] -= dIa_dVak + dIa_dVgk;
                    G[triode.cathode][triode.anode] -= dIa_dVak;
                    G[triode.cathode][triode.cathode] += dIa_dVak + dIa_dVgk;
                }
                if (triode.grid != 0 && triode.cathode != 0) {
                    G[triode.grid][triode.cathode] -= dIa_dVgk;
                    G[triode.cathode][triode.grid] -= dIa_dVgk;
                    G[triode.grid][triode.grid] += dIa_dVgk;
                }

                if (triode.anode != 0) b[triode.anode] -= I;
                if (triode.cathode != 0) b[triode.cathode] += I;
            }

            updateSystemMatrix();
            solveLU();
            // solve();

            double maxDelta = 0.0;
            for (int i = 0; i < kSystemSize; ++i) {
                double delta = std::fabs(x[i] - x_guess[i]);
                maxDelta = std::max(maxDelta, delta);
                x_guess[i] = x[i];
            }

            if (maxDelta < kTolerance) break;

            G = G_backup;
            b = b_backup;
        }
    }

    [[nodiscard]] static double triodeCurrent(double Vak, double Vgk) {
        constexpr double MU = 96.2;
        constexpr double EX = 1.437;
        constexpr double KG1 = 613.4;
        constexpr double KP = 740.3;
        constexpr double KVB = 1672.0;

        double E1 = Vak / KP * std::log1p(std::exp(KP * (1.0 / MU + Vgk / std::sqrt(KVB + Vak * Vak))));

        // Clamp E1 to prevent overflow
        E1 = std::max(0.0, std::min(E1, 10.0));

        double I = 0.5 * (std::pow(E1, EX) + std::pow(std::abs(E1), EX)) / KG1;

        // Realistic current limiting
        return std::max(0.0, std::min(I, 0.01)); // Max 10mA
    }
};
