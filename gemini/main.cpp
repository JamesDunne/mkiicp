
#include <cassert>
#include <iostream>
#include <ranges>

#include "sim.hpp"
#include "tonestack.h"
#include "tubemodel.h"

using StageV1A = TubeModel<TubeStageV1A>;
using StageV1B = TubeModel<TubeStageV1B>;
using StageV3B = TubeModel<TubeStageV3B>;

/**
 * @brief Parses a RIFF WAVE file, processes its samples, and writes a new WAVE file.
 *
 * @param inputFilename The path to the input .wav file.
 * @param outputFilename The path to write the output .wav file.
 * @param processor A function that takes a sample as a double [-1.0, 1.0] and returns a processed sample.
 * @throws std::runtime_error on file I/O errors or if the WAVE format is unsupported.
 */
void processWavFile(
    const std::string& inputFilename,
    const std::string& outputFilename,
    const std::function<double(double)>& processor
);

void setup_markiicp(RealtimeTubeSim& sim) {
    std::vector<std::string> nodes = {
        "N001", "N002", "N003", "N004", "N005", "N006", "N007", "N008", "N009", "N010",
        "N011", "N012", "N013", "N014", "N015", "N016", "N017", "N018", "N019", "N020",
        "N021", "N022", "N023", "N024", "N025", "N026", "N027", "N028", "N029", "N030",
        "N031", "N032", "N033", "N034", "N035", "N036", "P001"
    };

    for (const auto & node : nodes) {
        sim.add_node(node);
    }

    // Resistors
    sim.add_resistor("R5", "N004", "N015", 100e3);     // R5
    sim.add_resistor("R5A", "N008", "N007", 100e3);     // R5A
    sim.add_resistor("R4", "N003", "N004", 150e3);     // R4
    sim.add_resistor("R1", "N019", "0", 1e6);          // R1
    sim.add_resistor("R2", "N033", "0", 1.5e3);        // R2
    sim.add_resistor("R7", "N028", "0", 1.5e3);        // R7
    sim.add_resistor("R8", "N003", "N009", 100e3);     // R8
    sim.add_resistor("R9", "N001", "0", 100e3);        // R9
    sim.add_resistor("R21", "N011", "N025", 680e3);     // R21
    sim.add_resistor("R22", "N029", "0", 475e3);        // R22
    sim.add_resistor("R23", "N035", "0", 1.5e3);        // R23
    sim.add_resistor("R10", "N002", "N001", 3.3e6);     // R10
    sim.add_resistor("R11", "N002", "0", 680e3);        // R11
    sim.add_resistor("R26", "N017", "N010", 82e3);      // R26
    sim.add_resistor("R24", "N030", "0", 68e3);         // R24
    sim.add_resistor("R25", "N030", "N018", 270e3);     // R25
    sim.add_resistor("R30", "N034", "0", 3.3e3);        // R30
    sim.add_resistor("R27", "N026", "N010", 274e3);     // R27
    sim.add_resistor("R31", "N002", "N027", 220e3);     // R31
    sim.add_resistor("R16", "N036", "0", 1.5e3);        // R16
    sim.add_resistor("R13", "N006", "N021", 100e3);     // R13
    sim.add_resistor("R105", "N022", "P001", 47e3);      // R105
    // sim.add_resistor("R46", "N022", "0", 47e3);         // R46 = lead_master
    sim.add_resistor("R102", "N022", "N032", 150e3);     // R102
    sim.add_resistor("R101", "N032", "0", 4.7e3);        // R101
    sim.add_resistor("R19", "N006", "N012", 120e3);     // R19
    sim.add_resistor("R103", "N023", "0", 47e3);         // R103
    sim.add_resistor("R104", "N031", "0", 1e3);          // R104
    sim.add_resistor("R12", "N023", "N032", 2.2e3);     // R12
    sim.add_resistor("R106", "N014", "N013", 15e3);      // R106
    //sim.add_resistor("R6", "N005", "N005", 10e6);      // R6 (bypassed)
    sim.add_resistor("R32", "N002", "0", 100e3);        // R32

    // Capacitors
    sim.add_capacitor("C6", "N005", "N004", 750e-12);  // C6
    sim.add_capacitor("C5", "N005", "N004", 250e-12);  // C5 (parallel with C6)
    sim.add_capacitor("C4", "N016", "N015", 0.1e-6);   // C4
    sim.add_capacitor("C3", "N024", "N015", 0.047e-6); // C3
    sim.add_capacitor("C13B", "N020", "N008", 180e-12);  // C13B
    sim.add_capacitor("C1", "N033", "0", 0.47e-6);     // C1
    sim.add_capacitor("C2", "N033", "0", 22e-6);       // C2
    sim.add_capacitor("C13", "N028", "0", 22e-6);       // C13
    sim.add_capacitor("C7", "N001", "N009", 0.1e-6);   // C7
    sim.add_capacitor("C22", "N035", "N029", 120e-12);  // C22
    sim.add_capacitor("C23", "N035", "0", 2.2e-6);      // C23
    sim.add_capacitor("C10", "N002", "N001", 20e-12);   // C10
    sim.add_capacitor("C24", "N030", "0", 1000e-12);    // C24
    sim.add_capacitor("C29", "N034", "0", 0.22e-6);     // C29
    sim.add_capacitor("C30", "N027", "N026", 0.047e-6); // C30
    sim.add_capacitor("C31", "N002", "N027", 250e-12);  // C31
    sim.add_capacitor("C9", "P001", "N021", 0.047e-6); // C9
    sim.add_capacitor("C12", "N013", "N012", 0.047e-6); // C12
    sim.add_capacitor("C11", "N002", "0", 47e-12);      // C11
    sim.add_capacitor("C16", "N031", "0", 0.47e-6);     // C16
    sim.add_capacitor("C15", "N031", "0", 15e-6);       // C15
    sim.add_capacitor("C21", "N001", "N011", 0.02e-6);  // C21
    sim.add_capacitor("C32", "N002", "0", 500e-12);     // C32
    sim.add_capacitor("C25", "N018", "N017", 0.022e-6); // C25

    // input:
    sim.add_voltage_source("Vin", "N019", "0", 0.0, true);  // Vin

    // Voltage sources (DC bias)
    sim.add_voltage_source("VE", "N003", "0", 405, false); // VE
    sim.add_voltage_source("VC", "N010", "0", 410, false); // VC
    sim.add_voltage_source("VC2", "N006", "0", 410, false); // VC2

    // tone stack:
    sim.add_potentiometer("treble", "N005", "N016", "N007", 250e3, 'A');
    sim.add_variable_resistor("bass", "N016", "N024", 250e3, 'A');
    sim.add_variable_resistor("mid", "N024", "0", 10e3, 'A');

    // volume1:
    sim.add_potentiometer("volume1", "N008", "0", "N020", 1e6, 'A');
    // lead_drive:
    sim.add_potentiometer("lead_drive", "N025", "0", "N029", 1e6, 'A');
    // lead_master:
    sim.add_variable_resistor("lead_master", "N022", "0", 250e3, 'A');

    // master volume:
    sim.add_variable_resistor("master", "N014", "0", 1e6, 'A');

    // tube stages:
    // mu = 96.2
    // ex = 1.437
    // kg1 = 613.4
    // kp = 740.3
    // kvb = 1672.0
    // rgi = 2000.0

    // XV1A: 12AX7 N004 N019 N033
    sim.add_triode("XV1A", "N004", "N019", "N033");

    // XV1B: 12AX7 N009 N020 N028
    sim.add_triode("XV1B", "N009", "N020", "N028");

    // XV3B: 12AX7 N017 N029 N035
    sim.add_triode("XV3B", "N017", "N029", "N035");

    // XV4A: 12AX7 N026 N030 N034
    sim.add_triode("XV4A", "N026", "N030", "N034");

    // XV2B: 12AX7 N021 N002 N036
    sim.add_triode("XV2B", "N021", "N002", "N036");

    // XV2A: 12AX7 N012 N023 N031
    sim.add_triode("XV2A", "N012", "N023", "N031");

    // capture output from N014 after master volume:
    sim.set_output_node("N014");
}

void setup_basic_tube_preamp(RealtimeTubeSim& sim) {
    // single tube stage preamp with tone stack:

    // Define all nodes first
    sim.add_node("in"); sim.add_node("g"); sim.add_node("p"); sim.add_node("k");
    sim.add_node("V_P"); sim.add_node("ts1"); sim.add_node("ts2"); sim.add_node("ts3");
    sim.add_node("out"); sim.add_node("out_final");

    // --- Static Components ---
    sim.add_resistor("R1", "g", "gnd", 1e6);        // Grid leak
    sim.add_resistor("R2", "V_P", "p", 100e3);      // Plate load
    sim.add_resistor("R3", "k", "gnd", 1.5e3);       // Cathode resistor
    sim.add_resistor("R4", "out_final", "gnd", 1e6);// Final load
    sim.add_capacitor("C1", "in", "g", 22e-9);      // Input cap
    sim.add_capacitor("C2", "k", "gnd", 22e-6);     // Cathode bypass
    sim.add_capacitor("C3", "p", "ts1", 22e-9);     // Coupling cap to tone stack
    sim.add_capacitor("C4", "ts1", "ts3", 250e-12);  // Treble cap
    sim.add_capacitor("C5", "ts2", "ts3", 22e-9);    // Mid cap
    sim.add_capacitor("C6", "out", "out_final", 100e-9); // Output cap

    // --- Dynamic Components (Pots/Variable Resistors) ---
    sim.add_variable_resistor("treble", "ts1", "ts2", 250e3, 'L');
    sim.add_variable_resistor("bass", "ts2", "gnd", 1e6, 'A'); // Bass pot wired as rheostat
    sim.add_variable_resistor("mid", "ts3", "gnd", 25e3, 'L');
    sim.add_potentiometer("master", "ts3", "gnd", "out", 1e6, 'A');

    // --- Active Components ---
    sim.add_triode("V1A", "p", "g", "k", 96.2, 1.437, 613.4, 740.3, 1672.0, 2000.0);

    // --- Sources ---
    sim.add_voltage_source("V_P", "V_P", "gnd", 300.0, false);
    sim.add_voltage_source("Vin", "in", "gnd", 0.0, true);
}

void simulate_sine_sweep(RealtimeTubeSim& sim, const double sampleRate) {
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
        double phaseIncrement = 2.0 * M_PI * instantaneousFreq / sampleRate;

        // Add to accumulated phase and generate sample
        phase += phaseIncrement;
        double x = std::sin(phase);

        x = sim.process_sample(x);

        std::cout
                << t << ","
                << x
                << std::endl;
    }
}

void dump_netlist(const RealtimeTubeSim & sim) {
    for (const auto& r : sim.m_resistors) {
        assert(r.name.at(0) == 'R');
        std::cout << r.name << " " << sim.name_of(r.n1) << " " << sim.name_of(r.n2) << " " << r.R << std::endl;
    }
    for (const auto& c : sim.m_capacitors) {
        assert(c.name.at(0) == 'C');
        std::cout << c.name << " " << sim.name_of(c.n1) << " " << sim.name_of(c.n2) << " " << c.C << std::endl;
    }
    for (const auto &r: sim.m_variable_resistors | std::views::values) {
        std::cout << "Rv_" << r.name << " " << sim.name_of(r.n1) << " " << sim.name_of(r.n2) << " {" << r.R_max << "*(1-" << r.name << "}" << std::endl;
        double w = r.value;
        if (r.taper == 'A' || r.taper == 'a') {
            w = w * w;
        }
        std::cout << ".param " << r.name << " " << w << std::endl;
    }
    for (const auto &p: sim.m_potentiometers | std::views::values) {
        std::cout << "Rpa_" << p.name << " " << sim.name_of(p.n1) << " " << sim.name_of(p.wiper) << " {" << p.R_total << "*(1-" << p.name << "}" << std::endl;
        std::cout << "Rpc_" << p.name << " " << sim.name_of(p.wiper) << " " << sim.name_of(p.n2) << " {" << p.R_total << "*(" << p.name << "}" << std::endl;
        double w = p.value;
        if (p.taper == 'A' || p.taper == 'a') {
            w = w * w;
        }
        std::cout << ".param " << p.name << " " << w << std::endl;
    }
    for (const auto &t: sim.m_triodes) {
        assert(t.name.at(0) == 'X');
        std::cout << t.name << " " << sim.name_of(t.p_node) << " " << sim.name_of(t.g_node) << " " << sim.name_of(t.k_node) << " 12AX7" << std::endl;
    }
    for (const auto &v : sim.m_v_sources) {
        assert(v.name.at(0) == 'V');
        std::cout << v.name << " " << sim.name_of(v.node) << " 0";
        if (v.is_time_varying) {
            // input wave file:f
            std::cout << " wavefile=di-cut.wav";
        } else {
            std::cout << " " << v.dc_voltage;
        }
        std::cout << std::endl;
    }

    std::cout << "Ewout wout 0 " << sim.name_of(sim.m_output_node) << " 0 {1/1400}" << std::endl;
    // output wave file:
    std::cout << ".wave mkiicp-out.wav 16 48000 V(wout)" << std::endl;

    std::cout << ".inc 12AX7.cir" << std::endl;

    std::cout << ".tran 0.0000026 10" << std::endl;
    std::cout << ".backanno" << std::endl;
    std::cout << ".end" << std::endl;
}

int mnasim_main(int argc, char *argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input.wav> <output.wav>" << std::endl;
        return 1;
    }

    const double sampleRate = 48000.0;
    RealtimeTubeSim sim { sampleRate };

#if 0
    setup_basic_tube_preamp(sim);
#else
    setup_markiicp(sim);
#endif

    // Set controls for IIC+:
    // sim.set_parameter("volume1", 0.68556546); // hard-coded from Mark V IIC+ mode (470kOhm)
    sim.set_parameter("volume1", 0.8660254);
    sim.set_parameter("lead_drive", 0.70710678);
    sim.set_parameter("treble", 0.89442719);
    sim.set_parameter("mid", 0.57445626);
    sim.set_parameter("bass", 0.2236068);
    sim.set_parameter("lead_master", 0.90111043);
    sim.set_parameter("master", 0.70710678);

    sim.prepare_to_play();

    dump_netlist(sim);

    //simulate_sine_sweep(sim, sampleRate);

    try {
        const std::string input_filename = argv[1];
        const std::string output_filename = argv[2];

        double minV = 1.0;
        double maxV = -0.0;
        double t = 0.0;
        auto ampsim_process = [&](double sample) -> double {
            double out = sim.process_sample(sample) / 200.0;
#if MINMAX
            if (out > maxV) {
                maxV = out;
                std::cout << "min=" << minV << " max=" << maxV << std::endl;
            }
            if (out < minV) {
                minV = out;
                std::cout << "min=" << minV << " max=" << maxV << std::endl;
            }
#endif
#if CSV
            std::cout << t
                // << "," << sim.m_x.at(sim.m_node("N019"))
                // << "," << sim.m_x.at(sim.m_node("N007"))
                << "," << sim.m_x.at(sim.m_node("N004"))
                // << "," << sim.m_x.at(sim.m_node("N020"))
                // << "," << sim.m_x.at(sim.m_node("N029"))
                // << "," << sim.m_x.at(sim.m_node("N030"))
                // << "," << sim.m_x.at(sim.m_node("N002"))
                // << "," << sim.m_x.at(sim.m_node("N023"))
                // << "," << out
                << std::endl;
            t += 1.0 / sampleRate;
#endif
            return out;
        };

#if CSV
        std::cout << "t,N019,N007,N020,N029,N030,N002,N023,N014" << std::endl;
#endif
        processWavFile(input_filename, output_filename, ampsim_process);
#if MINMAX
        std::cout << "min=" << minV << " max=" << maxV << std::endl;
#endif
        std::cout << "Successfully processed " << input_filename << " to " << output_filename << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

int tonestack_main(int argc, char *argv[]) {
    const std::string input_filename = argv[1];
    const std::string output_filename = argv[2];

    PassiveToneStack tone;
    tone.setParameters(48000.0, 0.8, 0.1, 0.25, 0.9);

    // process input file and produce output file:
    double min = 1.0, max = -1.0;
    processWavFile(
        input_filename,
        output_filename,
        [&](double sample) -> double {
            double vout = tone.processSample(sample);
            double ac_out = vout;
            if (ac_out > max) { max = ac_out; }
            if (ac_out < min) { min = ac_out; }
            return ac_out / 1.0;
        }
    );
    std::cout << min << " " << max << std::endl;

    return 0;
}

int tubemodel_main(int argc, char *argv[]) {
    // Set precision for printing doubles
    std::cout << std::fixed << std::setprecision(6);

    std::cout << "--- TubeModel Verification ---" << std::endl;

    // 1. Test clamping regions
    double vg_below_clamp = -10.0;
    double vg_above_clamp = 20.0;
    std::cout << "Vg = " << vg_below_clamp << " -> Vp = " << StageV1A::processSample(vg_below_clamp) << " (Expected High Clamp)" << std::endl;
    std::cout << "Vg = " << vg_above_clamp << " -> Vp = " << StageV1A::processSample(vg_above_clamp) << " (Expected Low Clamp)" << std::endl;
    std::cout << std::endl;

    // 2. Test values within each segment
    double vg_seg1 = -6.0;
    double vg_seg2 = -4.0;
    double vg_seg3 = 0.0;
    double vg_seg4 = 4.0;
    double vg_seg5 = 10.0;
    std::cout << "--- Testing Mid-Segment Values ---" << std::endl;
    std::cout << "Vg = " << vg_seg1 << " (Seg 1) -> Vp = " << StageV1A::processSample(vg_seg1) << std::endl;
    std::cout << "Vg = " << vg_seg2 << " (Seg 2) -> Vp = " << StageV1A::processSample(vg_seg2) << std::endl;
    std::cout << "Vg = " << vg_seg3 << " (Seg 3) -> Vp = " << StageV1A::processSample(vg_seg3) << std::endl;
    std::cout << "Vg = " << vg_seg4 << " (Seg 4) -> Vp = " << StageV1A::processSample(vg_seg4) << std::endl;
    std::cout << "Vg = " << vg_seg5 << " (Seg 5) -> Vp = " << StageV1A::processSample(vg_seg5) << std::endl;
    std::cout << std::endl;

    // 3. Test values at the boundaries to check for continuity
    // The values should be very close, proving the C0 continuity constraint worked.
    double b1 = -5.677862;
    double b2 = -2.481141;
    double b3 = 1.722849;
    double b4 = 6.116238;
    std::cout << "--- Testing Boundary Values ---" << std::endl;
    std::cout << "Vg approaching " << b1 << " -> Vp = " << StageV1A::processSample(b1 - 1e-9) << std::endl;
    std::cout << "Vg at " << b1 << "         -> Vp = " << StageV1A::processSample(b1) << std::endl;
    std::cout << "Vg approaching " << b2 << " -> Vp = " << StageV1A::processSample(b2 - 1e-9) << std::endl;
    std::cout << "Vg at " << b2 << "         -> Vp = " << StageV1A::processSample(b2) << std::endl;
    std::cout << "Vg approaching " << b3 << " -> Vp = " << StageV1A::processSample(b3 - 1e-9) << std::endl;
    std::cout << "Vg at " << b3 << "         -> Vp = " << StageV1A::processSample(b3) << std::endl;
    std::cout << "Vg approaching " << b4 << " -> Vp = " << StageV1A::processSample(b4 - 1e-9) << std::endl;
    std::cout << "Vg at " << b4 << "         -> Vp = " << StageV1A::processSample(b4) << std::endl;

    return 0;
}

int tonestack_test_main() {
    const double sampleRate = 48000.0;
    const int numSamples = 10;

    // 1. Create the filter
    PassiveToneStack toneStack;

    // 2. Set initial parameters (e.g., all knobs at noon)
    std::cout << "Setting parameters: Treble=0.5, Bass=0.5, Mid=0.5, Volume=1.0\n";
    toneStack.setParameters(sampleRate, 0.5, 0.5, 0.5, 0.9);

    // Create a simple impulse signal (1.0 followed by zeros)
    std::vector<float> inputSignal(numSamples, 0.0f);
    inputSignal[0] = 1.0f;

    // 3. Process the signal
    std::cout << "Impulse Response:\n";
    for (int i = 0; i < numSamples; ++i) {
        float inputSample = inputSignal[i];
        float outputSample = toneStack.processSample(inputSample);
        std::cout << "Sample " << i << ": " << outputSample << std::endl;
    }

    // --- Example of changing a parameter ---
    std::cout << "\nScooping the mids: Treble=0.8, Bass=0.8, Mid=0.1\n";
    toneStack.setParameters(sampleRate, 0.8, 0.8, 0.1, 0.9);

    // Resetting the filter state is good practice after a large parameter change
    toneStack.reset();

    std::cout << "New Impulse Response:\n";
    for (int i = 0; i < numSamples; ++i) {
        float inputSample = inputSignal[i];
        float outputSample = toneStack.processSample(inputSample);
        std::cout << "Sample " << i << ": " << outputSample << std::endl;
    }

    return 0;
}

int stage1_main(int argc, char *argv[]) {
    PassiveToneStack tone;
    tone.setParameters(48000.0, 0.67, 0.15, 0.25, 0.86);

    const std::string input_filename = argv[1];
    const std::string output_filename = argv[2];

    // process input file and produce output file:
    double min = 405.0, max = -405.0;
    processWavFile(
        input_filename,
        output_filename,
        [&](double sample) -> double {
            // double vout = (StageV1A::processSample(sample) - 2.08157285e+02) / 45.0;
            double vout = StageV1A::processSample(sample);
            vout -= StageV1A::Params::SEG3_COEFFS[2]; // remove DC offset
            vout = tone.processSample(vout); // apply tone stack
            vout = StageV1B::processSample(vout);
            vout -= StageV1B::Params::SEG3_COEFFS[2]; // remove DC offset
            vout /= 47.0;
            vout = StageV3B::processSample(vout);
            vout = (vout - StageV3B::Params::SEG3_COEFFS[2]) / 200.0; // remove DC offset; scale to unity range
            double ac_out = vout;
            if (ac_out > max) { max = ac_out; }
            if (ac_out < min) { min = ac_out; }
            return ac_out;
        }
    );
    std::cout << min << " " << max << std::endl;

    return 0;
}

int main(int argc, char *argv[]) {
    // return tubemodel_main(argc, argv);
    // return tonestack_test_main();
    // return tonestack_main(argc, argv);
    return stage1_main(argc, argv);
}
