#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <memory>

class TubeTriode {
private:
    // 12AX7 parameters from the model
    static constexpr double MU = 96.20;
    static constexpr double EX = 1.437;
    static constexpr double KG1 = 613.4;
    static constexpr double KP = 740.3;
    static constexpr double KVB = 1672.0;
    static constexpr double RGI = 2000.0;
    
    // Internal state for numerical integration
    double lastPlateVoltage;
    double lastGridVoltage;
    double lastCathodeVoltage;
    
public:
    TubeTriode() : lastPlateVoltage(0), lastGridVoltage(0), lastCathodeVoltage(0) {}
    
    double computePlateCurrent(double plateVoltage, double gridVoltage, double cathodeVoltage) {
        double Vpk = plateVoltage - cathodeVoltage;
        double Vgk = gridVoltage - cathodeVoltage;
        
        // Prevent numerical issues
        if (Vpk < 0.1) Vpk = 0.1;
        
        // Koren model implementation
        double E1 = Vpk / KP * log(1.0 + exp(KP * (1.0/MU + Vgk / sqrt(KVB + Vpk * Vpk))));
        
        // Clamp E1 to prevent overflow
        E1 = std::max(0.0, std::min(E1, 10.0));
        
        double plateCurrent = 0.5 * (pow(E1, EX) + pow(std::abs(E1), EX)) / KG1;
        
        // Realistic current limiting
        return std::max(0.0, std::min(plateCurrent, 0.01)); // Max 10mA
    }
    
    double computeGridCurrent(double gridVoltage, double cathodeVoltage) {
        double Vgk = gridVoltage - cathodeVoltage;
        if (Vgk > 0.7) {
            // Simple diode model for grid conduction
            return (exp(Vgk - 0.7) - 1.0) * 1e-12;
        }
        return 0.0;
    }
};

class CapacitorFilter {
private:
    double capacitance;
    double lastVoltage;
    double sampleRate;
    
public:
    CapacitorFilter(double cap, double fs) : capacitance(cap), lastVoltage(0), sampleRate(fs) {}
    
    double process(double input, double resistance) {
        double dt = 1.0 / sampleRate;
        double tau = resistance * capacitance;
        double alpha = dt / (tau + dt);
        
        lastVoltage = lastVoltage * (1.0 - alpha) + input * alpha;
        return lastVoltage;
    }
    
    void reset() { lastVoltage = 0; }
};

class ToneStack {
private:
    CapacitorFilter C3, C4, C5, C6;
    double treble, mid, bass;
    
public:
    ToneStack(double fs) : 
        C3(0.047e-6, fs), C4(0.1e-6, fs), C5(250e-12, fs), C6(750e-12, fs),
        treble(0.8), mid(0.33), bass(0.05) {}
    
    void setControls(double t, double m, double b) {
        treble = std::max(0.0, std::min(1.0, t));
        mid = std::max(0.0, std::min(1.0, m));
        bass = std::max(0.0, std::min(1.0, b));
    }
    
    double process(double input) {
        // Simplified tone stack simulation
        double trebleR = 250000.0 * (1.0 - treble);
        double midR = 10000.0 * (1.0 - mid);
        double bassR = 250000.0 * (1.0 - bass);
        
        // High frequency response (treble)
        double highFreq = C5.process(input, trebleR + 1000.0);
        
        // Mid frequency response
        double midFreq = C4.process(input, midR + 10000.0);
        
        // Low frequency response (bass)
        double lowFreq = C3.process(input, bassR + 10000.0);
        
        // Combine responses with frequency weighting
        return (highFreq * 0.3 + midFreq * 0.4 + lowFreq * 0.3);
    }
};

class TubeAmplifier {
private:
    // Tube stages
    std::vector<std::unique_ptr<TubeTriode>> tubes;
    
    // Filters for coupling capacitors and tone stack
    std::vector<std::unique_ptr<CapacitorFilter>> couplingCaps;
    std::unique_ptr<ToneStack> toneStack;
    
    // Control parameters
    double volume1, gain, master;
    double sampleRate;
    
    // Node voltages (simplified - key nodes only)
    std::vector<double> nodeVoltages;
    
    // Supply voltages
    static constexpr double VE = 405.0;  // Plate supply
    static constexpr double VC = 410.0;  // Screen supply
    static constexpr double VC2 = 410.0; // Additional supply
    
public:
    TubeAmplifier(double fs = 48000.0) : sampleRate(fs) {
        // Initialize 5 tube stages (V1A, V1B, V2A, V2B, V3B, V4A)
        for (int i = 0; i < 6; i++) {
            tubes.push_back(std::make_unique<TubeTriode>());
        }
        
        // Initialize coupling capacitors
        couplingCaps.push_back(std::make_unique<CapacitorFilter>(0.47e-6, fs));  // C1
        couplingCaps.push_back(std::make_unique<CapacitorFilter>(22e-6, fs));    // C2
        couplingCaps.push_back(std::make_unique<CapacitorFilter>(0.1e-6, fs));   // C7
        couplingCaps.push_back(std::make_unique<CapacitorFilter>(0.047e-6, fs)); // C9
        couplingCaps.push_back(std::make_unique<CapacitorFilter>(0.047e-6, fs)); // C12
        couplingCaps.push_back(std::make_unique<CapacitorFilter>(0.022e-6, fs)); // C25
        
        // Initialize tone stack
        toneStack = std::make_unique<ToneStack>(fs);
        
        // Initialize node voltages
        nodeVoltages.resize(50, 0.0);
        
        // Set default control values
        volume1 = 0.75;
        gain = 0.5;
        master = 0.5;
    }
    
    void setControls(double vol1, double g, double m, double treble, double mid, double bass) {
        volume1 = std::max(0.0, std::min(1.0, vol1));
        gain = std::max(0.0, std::min(1.0, g));
        master = std::max(0.0, std::min(1.0, m));
        toneStack->setControls(treble, mid, bass);
    }
    
    double processample(double input) {
        // Input stage with volume control
        double inputGain = volume1 * 2.0; // Scaling factor
        double scaledInput = input * inputGain;
        
        // V1A - First preamp stage
        double v1aGrid = scaledInput;
        double v1aPlate = VE - tubes[0]->computePlateCurrent(VE, v1aGrid, 0) * 150000.0; // R4
        double v1aCoupled = couplingCaps[0]->process(v1aPlate, 1000000.0); // R1
        
        // Apply tone stack
        double tonedSignal = toneStack->process(v1aCoupled);
        
        // V1B - Second preamp stage  
        double v1bGrid = tonedSignal;
        double v1bPlate = VE - tubes[1]->computePlateCurrent(VE, v1bGrid, 0) * 100000.0; // R8
        double v1bCoupled = couplingCaps[1]->process(v1bPlate, 100000.0); // R9
        
        // V2A - Third stage with gain control
        double gainFactor = gain * 3.0 + 0.5; // Gain control scaling
        double v2aGrid = v1bCoupled * gainFactor;
        double v2aPlate = VC2 - tubes[2]->computePlateCurrent(VC2, v2aGrid, 0) * 120000.0; // R19
        double v2aCoupled = couplingCaps[2]->process(v2aPlate, 47000.0); // R103
        
        // V2B - Fourth stage
        double v2bGrid = v2aCoupled;
        double v2bPlate = VC2 - tubes[3]->computePlateCurrent(VC2, v2bGrid, 0) * 100000.0; // R13
        double v2bCoupled = couplingCaps[3]->process(v2bPlate, 47000.0); // R105
        
        // V3B - Fifth stage (lead drive)
        double v3bGrid = v2bCoupled;
        double v3bPlate = VC - tubes[4]->computePlateCurrent(VC, v3bGrid, 0) * 82000.0; // R26
        double v3bCoupled = couplingCaps[4]->process(v3bPlate, 68000.0); // R24
        
        // V4A - Final preamp stage
        double v4aGrid = v3bCoupled;
        double v4aPlate = VC - tubes[5]->computePlateCurrent(VC, v4aGrid, 0) * 274000.0; // R27
        double v4aCoupled = couplingCaps[5]->process(v4aPlate, 220000.0); // R31
        
        // Master volume control and output scaling
        double output = v4aCoupled * master;
        
        // Apply soft clipping for tube saturation
        output = softClip(output);
        
        // Final output scaling to match expected levels
        return output / 1400.0; // Similar to E1 scaling in original netlist
    }
    
private:
    double softClip(double input) {
        // Tube-like soft clipping
        double threshold = 2.0;
        if (std::abs(input) < threshold) {
            return input;
        } else {
            double sign = (input > 0) ? 1.0 : -1.0;
            double excess = std::abs(input) - threshold;
            return sign * (threshold + excess / (1.0 + excess));
        }
    }
};

// Example usage and real-time processing framework
class AudioProcessor {
private:
    TubeAmplifier amplifier;
    std::vector<double> inputBuffer;
    std::vector<double> outputBuffer;
    size_t bufferSize;
    
public:
    AudioProcessor(size_t bufSize = 256) : bufferSize(bufSize) {
        inputBuffer.resize(bufferSize);
        outputBuffer.resize(bufferSize);
        
        // Set up amplifier with typical settings
        amplifier.setControls(
            0.75,  // volume1
            0.5,   // gain
            0.5,   // master
            0.8,   // treble
            0.33,  // mid
            0.05   // bass
        );
    }
    
    void processBlock(const double* input, double* output, size_t numSamples) {
        for (size_t i = 0; i < numSamples; ++i) {
            output[i] = amplifier.processample(input[i]);
        }
    }
    
    void setAmpControls(double volume1, double gain, double master, 
                       double treble, double mid, double bass) {
        amplifier.setControls(volume1, gain, master, treble, mid, bass);
    }
    
    // Demo function to show processing
    void demonstrateProcessing() {
        std::cout << "Mesa Boogie Mark IIC+ Tube Amplifier Simulator\n";
        std::cout << "==============================================\n\n";
        
        // Generate test signal (sine wave)
        const size_t testLength = 1000;
        std::vector<double> testInput(testLength);
        std::vector<double> testOutput(testLength);
        
        for (size_t i = 0; i < testLength; ++i) {
            // 440 Hz sine wave
            testInput[i] = 0.75 * sin(2.0 * M_PI * 440.0 * static_cast<double>(i) / 48000.0);
        }

        // Process the test signal
        processBlock(testInput.data(), testOutput.data(), testLength);
        
        // Show some statistics
        double maxInput = *std::max_element(testInput.begin(), testInput.end());
        double maxOutput = *std::max_element(testOutput.begin(), testOutput.end());
        double minOutput = *std::min_element(testOutput.begin(), testOutput.end());
        
        std::cout << "Test Signal Processing Results:\n";
        std::cout << "Input peak: " << maxInput << "\n";
        std::cout << "Output peak: " << maxOutput << "\n";
        std::cout << "Output range: " << minOutput << " to " << maxOutput << "\n";
        std::cout << "Gain: " << (maxOutput / maxInput) << "x\n\n";
        
        // Test different settings
        std::cout << "Testing different amp settings:\n";
        std::cout << "--------------------------------\n";
        
        struct Settings {
            const char* name;
            double vol, gain, master, treble, mid, bass;
        };
        
        Settings presets[] = {
            {"Clean", 0.3, 0.2, 0.4, 0.6, 0.5, 0.4},
            {"Crunch", 0.6, 0.5, 0.5, 0.7, 0.6, 0.3},
            {"Lead", 0.8, 0.8, 0.6, 0.8, 0.4, 0.2},
            {"High Gain", 0.9, 0.9, 0.7, 0.9, 0.3, 0.1}
        };
        
        for (const auto& preset : presets) {
            setAmpControls(preset.vol, preset.gain, preset.master, 
                          preset.treble, preset.mid, preset.bass);
            
            // Process a single sample to show response
            double testSample = 0.75;
            double result = amplifier.processample(testSample);
            
            std::cout << preset.name << " setting: " 
                      << "Input=" << testSample << " -> Output=" << result 
                      << " (Gain=" << (result/testSample) << "x)\n";
        }
    }
};

int main() {
    //AudioProcessor processor;
    //processor.demonstrateProcessing();
    
    std::cout << "\nSimulator ready for real-time audio processing!\n";
    std::cout << "Use processBlock() method for real-time audio.\n";

    float sampleRate = 48000.0f;    // Audio sample rate in Hz

    TubeAmplifier amplifier(sampleRate);

    float startFreq = 20.0f;        // Start frequency in Hz
    float endFreq = 20000.0f;       // End frequency in Hz
    float duration = 5.0f;          // Sweep duration in seconds

    float k = (endFreq / startFreq);
    // float L = duration / std::log(k);

    double phase = 0.0;
    for (long i = 0; i < 48000 * 5; ++i) {
        float t = static_cast<float>(i) / sampleRate;
        float instantaneousFreq = startFreq * std::pow(k, t / duration);

        // Calculate phase increment for this sample
        float phaseIncrement = 2.0f * static_cast<float>(M_PI) * instantaneousFreq / sampleRate;

        // Add to accumulated phase and generate sample
        phase += phaseIncrement;
        float x = std::sin(phase);

        // process sample with first gain stage and tone stack:
        x = amplifier.processample(x);

        std::cout
            << t << ","
            << x
            << std::endl;
    }

    return 0;
}
