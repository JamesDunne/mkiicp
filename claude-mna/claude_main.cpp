#pragma once

#include <vector>
#include <map>
#include <string>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <memory>
#include <ostream>

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
    int numNodes;
    
    // MNA matrices
    std::vector<std::vector<double>> G; // Conductance matrix
    std::vector<std::vector<double>> C; // Capacitance matrix
    std::vector<double> b;              // Current source vector
    std::vector<double> x;              // Solution vector (node voltages)
    std::vector<double> x_prev;         // Previous solution for integration
    
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
    struct TubeStage {
        int anode, grid, cathode;
        double va_prev = 0, vg_prev = 0, vc_prev = 0;
        double ia_prev = 0;
    };
    
    std::vector<TubeStage> tubes;
    
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
    
    void setSampleRate(double fs) {
        sampleRate = fs;
        dt = 1.0 / fs;
        updateCapacitorStamps();
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
        int inputNode = nodeMap["N019"];
        b[inputNode] = input;
        
        // Newton-Raphson iteration for nonlinear elements (tubes)
        const int maxIterations = 5;
        const double tolerance = 1e-6;
        
        for (int iter = 0; iter < maxIterations; iter++) {
            // Solve linear system Ax = b
            solveLU();
            
            // Update tube nonlinearities
            bool converged = updateTubeStages(tolerance);
            if (converged) break;
        }
        
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
        
        for (size_t i = 0; i < nodes.size(); i++) {
            nodeMap[nodes[i]] = i + 1;
        }
        
        numNodes = nodes.size() + 1; // +1 for ground
        
        // Initialize matrices
        G.assign(numNodes, std::vector<double>(numNodes, 0.0));
        C.assign(numNodes, std::vector<double>(numNodes, 0.0));
        b.assign(numNodes, 0.0);
        x.assign(numNodes, 0.0);
        x_prev.assign(numNodes, 0.0);
        
        // Initialize LU decomposition storage
        LU.assign(numNodes, std::vector<double>(numNodes, 0.0));
        pivot.assign(numNodes, 0);
        
        // Stamp passive components
        stampPassiveComponents();
        
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
        stampResistor("N006", "0", 410e3);        // VC2 equivalent
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
        
        // Voltage sources (DC bias)
        stampVoltageSource("N003", "0", 405);     // VE
        stampVoltageSource("N010", "0", 410);     // VC
        
        // Tone stack and gain controls (initial values)
        updateToneStack();
        updateGainControl();
        updateMasterControl();
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
    
    void stampVoltageSource(const std::string& nPos, const std::string& nNeg, double voltage) {
        int pos = nodeMap[nPos];
        int neg = nodeMap[nNeg];
        
        // Simple voltage source implementation - add large conductance
        double largeG = 1e6;
        if (pos != 0) {
            G[pos][pos] += largeG;
            b[pos] += largeG * voltage;
        }
        if (neg != 0) {
            G[neg][neg] += largeG;
            b[neg] -= largeG * voltage;
        }
    }
    
    void updateToneStack() {
        // Clear previous tone stack stamps
        clearResistorRange("RA_VOLUME1", "RC_VOLUME1");
        clearResistorRange("RA_TREBLE", "RC_TREBLE");
        clearResistorRange("RA_BASS", "RC_BASS");
        clearResistorRange("RA_MID", "RC_MID");
        
        // Volume control
        double vol1 = params.vol1;
        stampResistor("N020", "N008", 1e6 * (1 - vol1));  // RA_VOLUME1
        stampResistor("0", "N020", 1e6 * vol1);           // RC_VOLUME1
        
        // Treble control
        double treble = params.treble;
        stampResistor("N005", "N007", 250e3 * (1 - treble)); // RA_TREBLE
        stampResistor("N007", "N016", 250e3 * treble);        // RC_TREBLE
        
        // Bass control
        double bass = params.bass;
        stampResistor("N016", "N024", 250e3 * (1 - bass)); // RA_BASS
        stampResistor("N024", "N024", 250e3 * bass);        // RC_BASS
        
        // Mid control
        double mid = params.mid;
        stampResistor("N024", "0", 10e3 * (1 - mid)); // RA_MID
        stampResistor("0", "0", 10e3 * mid);           // RC_MID
    }
    
    void updateGainControl() {
        clearResistorRange("RA_GAIN", "RC_GAIN");
        
        double gain = params.gain;
        stampResistor("N025", "N029", 1e6 * (1 - gain)); // RA_GAIN
        stampResistor("N029", "0", 1e6 * gain);          // RC_GAIN
    }
    
    void updateMasterControl() {
        clearResistorRange("RA_MASTER", "RC_MASTER");
        
        double master = params.master;
        stampResistor("N014", "0", 1e6 * (1 - master)); // RA_MASTER
        stampResistor("0", "0", 1e6 * master);          // RC_MASTER
    }
    
    void clearResistorRange(const std::string& start, const std::string& end) {
        // Simple implementation - just rebuild the entire G matrix
        // In practice, you'd want to track and remove specific stamps
        for (auto& row : G) {
            std::fill(row.begin(), row.end(), 0.0);
        }
        stampPassiveComponents();
    }
    
    void initializeTubeStages() {
        // Add tube stages based on SPICE netlist
        tubes.resize(6);
        
        // XV1A: 12AX7 N004 N019 N033
        tubes[0] = {nodeMap["N004"], nodeMap["N019"], nodeMap["N033"], 0, 0, 0, 0};
        
        // XV1B: 12AX7 N009 N020 N028  
        tubes[1] = {nodeMap["N009"], nodeMap["N020"], nodeMap["N028"], 0, 0, 0, 0};
        
        // XV3B: 12AX7 N017 N029 N035
        tubes[2] = {nodeMap["N017"], nodeMap["N029"], nodeMap["N035"], 0, 0, 0, 0};
        
        // XV4A: 12AX7 N026 N030 N034
        tubes[3] = {nodeMap["N026"], nodeMap["N030"], nodeMap["N034"], 0, 0, 0, 0};
        
        // XV2B: 12AX7 N021 N002 N036
        tubes[4] = {nodeMap["N021"], nodeMap["N002"], nodeMap["N036"], 0, 0, 0, 0};
        
        // XV2A: 12AX7 N012 N023 N031
        tubes[5] = {nodeMap["N012"], nodeMap["N023"], nodeMap["N031"], 0, 0, 0, 0};
    }
    
    void updateCapacitorStamps() {
        // Update the system matrix to include capacitor contributions for current timestep
        updateSystemMatrix();
    }
    
    void updateSystemMatrix() {
        // Combine G and C matrices with backward Euler integration: (G + C/dt) * x = b + C/dt * x_prev
        for (int i = 0; i < numNodes; i++) {
            for (int j = 0; j < numNodes; j++) {
                LU[i][j] = G[i][j] + C[i][j] / dt;
            }
        }
        
        // Perform LU decomposition
        luDecompose();
    }
    
    void luDecompose() {
        for (int i = 0; i < numNodes; i++) {
            pivot[i] = i;
        }
        
        for (int k = 0; k < numNodes - 1; k++) {
            // Find pivot
            int maxRow = k;
            for (int i = k + 1; i < numNodes; i++) {
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
            for (int i = k + 1; i < numNodes; i++) {
                if (std::abs(LU[k][k]) > 1e-12) {
                    LU[i][k] /= LU[k][k];
                    for (int j = k + 1; j < numNodes; j++) {
                        LU[i][j] -= LU[i][k] * LU[k][j];
                    }
                }
            }
        }
    }
    
    void solveLU() {
        // Add capacitor current contribution to RHS
        std::vector<double> rhs = b;
        for (int i = 0; i < numNodes; i++) {
            for (int j = 0; j < numNodes; j++) {
                rhs[i] += C[i][j] * x_prev[j] / dt;
            }
        }
        
        // Forward substitution
        for (int i = 0; i < numNodes; i++) {
            x[i] = rhs[pivot[i]];
            for (int j = 0; j < i; j++) {
                x[i] -= LU[i][j] * x[j];
            }
        }
        
        // Back substitution
        for (int i = numNodes - 1; i >= 0; i--) {
            for (int j = i + 1; j < numNodes; j++) {
                x[i] -= LU[i][j] * x[j];
            }
            if (std::abs(LU[i][i]) > 1e-12) {
                x[i] /= LU[i][i];
            }
        }
    }
    
    bool updateTubeStages(double tolerance) {
        bool converged = true;
        
        for (auto& tube : tubes) {
            double va = x[tube.anode];
            double vg = x[tube.grid];
            double vc = x[tube.cathode];
            
            // 12AX7 triode model equations
            double vgk = vg - vc;
            double vak = va - vc;
            
            // Grid current (diode model for grid conduction)
            double ig = 0.0;
            if (vgk > 0.7) {
                ig = (vgk - 0.7) / tubeParams.RGI;
            }
            
            // Plate current using Koren model
            double E1 = vak / tubeParams.KP * std::log(1.0 + std::exp(tubeParams.KP * 
                (1.0 / tubeParams.MU + vgk / std::sqrt(tubeParams.KVB + vak * vak))));
            
            double ia = 0.0;
            if (E1 > 0) {
                ia = std::pow(E1, tubeParams.EX) / tubeParams.KG1;
            }
            
            // Check convergence
            if (std::abs(ia - tube.ia_prev) > tolerance) {
                converged = false;
            }
            
            // Update nonlinear stamps
            updateTubeStamp(tube, ia, ig);
            
            // Store for next iteration
            tube.va_prev = va;
            tube.vg_prev = vg;
            tube.vc_prev = vc;
            tube.ia_prev = ia;
        }
        
        return converged;
    }
    
    void updateTubeStamp(const TubeStage& tube, double ia, double ig) {
        // Norton equivalent: current source + conductance
        // Simplified linearization around operating point
        double gm = 1e-3; // Simplified transconductance
        double rp = 100e3; // Simplified plate resistance
        
        // Stamp current sources
        if (tube.anode != 0) b[tube.anode] -= ia;
        if (tube.cathode != 0) b[tube.cathode] += ia;
        if (tube.grid != 0) b[tube.grid] -= ig;
        
        // Stamp linearized conductances (simplified)
        double gp = 1.0 / rp;
        if (tube.anode != 0) {
            G[tube.anode][tube.anode] += gp;
            if (tube.cathode != 0) G[tube.anode][tube.cathode] -= gp;
        }
        if (tube.cathode != 0) {
            G[tube.cathode][tube.cathode] += gp;
            if (tube.anode != 0) G[tube.cathode][tube.anode] -= gp;
        }
    }
};

int main(int argc, char *argv[]) {
    const double sampleRate = 48000.0;
    MKIICSimulator sim { sampleRate };

    double startFreq = 20.0f;        // Start frequency in Hz
    double endFreq = 20000.0f;       // End frequency in Hz
    double duration = 5.0f;          // Sweep duration in seconds

    double k = (endFreq / startFreq);
    // float L = duration / std::log(k);

    double phase = 0.0;
    for (long i = 0; i < 48000 * 5; ++i) {
        double t = static_cast<double>(i) / sampleRate;
        double instantaneousFreq = startFreq * std::pow(k, t / duration);

        // Calculate phase increment for this sample
        double phaseIncrement = 2.0f * static_cast<double>(M_PI) * instantaneousFreq / sampleRate;

        // Add to accumulated phase and generate sample
        phase += phaseIncrement;
        double x = std::sin(phase);

        // process sample with first gain stage and tone stack:
        x = sim.processSample(x);

        std::cout
            << t << ","
            << x
            << std::endl;
    }

    return 0;
}