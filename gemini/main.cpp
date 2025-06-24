
#include <iostream>

#include "sim.hpp"

int main(int argc, char *argv[]) {
    const double sampleRate = 48000.0;
    RealtimeTubeSim sim { sampleRate };

#if 0
    // single tube stage preamp with tone stack:

    // Define all nodes first
    sim.add_node("in"); sim.add_node("g"); sim.add_node("p"); sim.add_node("k");
    sim.add_node("V_P"); sim.add_node("ts1"); sim.add_node("ts2"); sim.add_node("ts3");
    sim.add_node("out"); sim.add_node("out_final");

    // --- Static Components ---
    sim.add_resistor("g", "gnd", 1e6);        // Grid leak
    sim.add_resistor("V_P", "p", 100e3);      // Plate load
    sim.add_resistor("k", "gnd", 1.5e3);       // Cathode resistor
    sim.add_resistor("out_final", "gnd", 1e6);// Final load
    sim.add_capacitor("in", "g", 22e-9);      // Input cap
    sim.add_capacitor("k", "gnd", 22e-6);     // Cathode bypass
    sim.add_capacitor("p", "ts1", 22e-9);     // Coupling cap to tone stack
    sim.add_capacitor("ts1", "ts3", 250e-12);  // Treble cap
    sim.add_capacitor("ts2", "ts3", 22e-9);    // Mid cap
    sim.add_capacitor("out", "out_final", 100e-9); // Output cap

    // --- Dynamic Components (Pots/Variable Resistors) ---
    sim.add_variable_resistor("treble", "ts1", "ts2", 250e3, 'L');
    sim.add_variable_resistor("bass", "ts2", "gnd", 1e6, 'A'); // Bass pot wired as rheostat
    sim.add_variable_resistor("mid", "ts3", "gnd", 25e3, 'L');
    sim.add_potentiometer("master", "ts3", "gnd", "out", 1e6, 'A');

    // --- Active Components ---
    sim.add_triode("p", "g", "k", 96.2, 1.437, 613.4, 740.3, 1672.0, 2000.0);

    // --- Sources ---
    sim.add_voltage_source("V_P", "gnd", 300.0, false);
    sim.add_voltage_source("in", "gnd", 0.0, true);
#else
    std::vector<std::string> nodes = {
        "N001", "N002", "N003", "N004", "N005", "N006", "N007", "N008", "N009", "N010",
        "N011", "N012", "N013", "N014", "N015", "N016", "N017", "N018", "N019", "N020",
        "N021", "N022", "N023", "N024", "N025", "N026", "N027", "N028", "N029", "N030",
        "N031", "N032", "N033", "N034", "N035", "N036", "P001"
    };

    for (const auto & node : nodes) {
        sim.add_node(node);
    }

    // Resistors - stamp conductances
    sim.add_resistor("N004", "N015", 100e3);     // R5
    sim.add_resistor("N008", "N007", 100e3);     // R5A
    sim.add_resistor("N003", "N004", 150e3);     // R4
    sim.add_resistor("N019", "0", 1e6);          // R1
    sim.add_resistor("N033", "0", 1.5e3);        // R2
    sim.add_resistor("N028", "0", 1.5e3);        // R7
    sim.add_resistor("N003", "N009", 100e3);     // R8
    sim.add_resistor("N001", "0", 100e3);        // R9
    sim.add_resistor("N011", "N025", 680e3);     // R21
    sim.add_resistor("N029", "0", 475e3);        // R22
    sim.add_resistor("N035", "0", 1.5e3);        // R23
    sim.add_resistor("N002", "N001", 3.3e6);     // R10
    sim.add_resistor("N002", "0", 680e3);        // R11
    sim.add_resistor("N017", "N010", 82e3);      // R26
    sim.add_resistor("N030", "0", 68e3);         // R24
    sim.add_resistor("N030", "N018", 270e3);     // R25
    sim.add_resistor("N034", "0", 3.3e3);        // R30
    sim.add_resistor("N026", "N010", 274e3);     // R27
    sim.add_resistor("N002", "N027", 220e3);     // R31
    sim.add_resistor("N036", "0", 1.5e3);        // R16
    sim.add_resistor("N006", "N021", 100e3);     // R13
    sim.add_resistor("N022", "P001", 47e3);      // R105
    // sim.add_resistor("N022", "0", 47e3);         // R46 = lead_master
    sim.add_resistor("N022", "N032", 150e3);     // R102
    sim.add_resistor("N032", "0", 4.7e3);        // R101
    sim.add_resistor("N006", "N012", 120e3);     // R19
    sim.add_resistor("N023", "0", 47e3);         // R103
    sim.add_resistor("N031", "0", 1e3);          // R104
    sim.add_resistor("N023", "N032", 2.2e3);     // R12
    sim.add_resistor("N014", "N013", 15e3);      // R106
    sim.add_resistor("N005", "N005", 10e6);      // R6
    sim.add_resistor("N002", "0", 100e3);        // R32

    // Capacitors
    sim.add_capacitor("N005", "N004", 750e-12);  // C6
    sim.add_capacitor("N005", "N004", 250e-12);  // C5 (parallel with C6)
    sim.add_capacitor("N016", "N015", 0.1e-6);   // C4
    sim.add_capacitor("N024", "N015", 0.047e-6); // C3
    sim.add_capacitor("N020", "N008", 180e-12);  // C13B
    sim.add_capacitor("N033", "0", 0.47e-6);     // C1
    sim.add_capacitor("N033", "0", 22e-6);       // C2
    sim.add_capacitor("N028", "0", 22e-6);       // C13
    sim.add_capacitor("N001", "N009", 0.1e-6);   // C7
    sim.add_capacitor("N035", "N029", 120e-12);  // C22
    sim.add_capacitor("N035", "0", 2.2e-6);      // C23
    sim.add_capacitor("N002", "N001", 20e-12);   // C10
    sim.add_capacitor("N030", "0", 1000e-12);    // C24
    sim.add_capacitor("N034", "0", 0.22e-6);     // C29
    sim.add_capacitor("N027", "N026", 0.047e-6); // C30
    sim.add_capacitor("N002", "N027", 250e-12);  // C31
    sim.add_capacitor("P001", "N021", 0.047e-6); // C9
    sim.add_capacitor("N013", "N012", 0.047e-6); // C12
    sim.add_capacitor("N002", "0", 47e-12);      // C11
    sim.add_capacitor("N031", "0", 0.47e-6);     // C16
    sim.add_capacitor("N031", "0", 15e-6);       // C15
    sim.add_capacitor("N001", "N011", 0.02e-6);  // C21
    sim.add_capacitor("N002", "0", 500e-12);     // C32
    sim.add_capacitor("N018", "N017", 0.022e-6); // C25

    // input:
    sim.add_voltage_source("N019", "0", 1.0, true);  // Vin

    // Voltage sources (DC bias)
    sim.add_voltage_source("N003", "0", 405, false); // VE
    sim.add_voltage_source("N010", "0", 410, false); // VC
    sim.add_voltage_source("N006", "0", 410, false); // VC2

    // tone stack:
    sim.add_potentiometer("treble", "N005", "N007", "N016", 250e3, 'A');
    sim.add_variable_resistor("bass", "N016", "0", 250e3, 'A');
    sim.add_variable_resistor("mid", "N024", "0", 10e3, 'A');

    // volume1:
    sim.add_potentiometer("volume1", "N008", "N020", "0", 1e6, 'A');
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
    sim.add_triode("N004", "N019", "N033");

    // XV1B: 12AX7 N009 N020 N028
    sim.add_triode("N009", "N020", "N028");

    // XV3B: 12AX7 N017 N029 N035
    sim.add_triode("N017", "N029", "N035");

    // XV4A: 12AX7 N026 N030 N034
    sim.add_triode("N026", "N030", "N034");

    // XV2B: 12AX7 N021 N002 N036
    sim.add_triode("N021", "N002", "N036");

    // XV2A: 12AX7 N012 N023 N031
    sim.add_triode("N012", "N023", "N031");
#endif

    sim.prepare_to_play();

    // Set controls
    sim.set_parameter("volume1", 0.75);
    sim.set_parameter("treble", 0.8);
    sim.set_parameter("mid", 0.33);
    sim.set_parameter("bass", 0.05);
    sim.set_parameter("lead_master", 0.75);
    sim.set_parameter("master", 0.5);

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
        double phaseIncrement = 2.0f * M_PI * instantaneousFreq / sampleRate;

        // Add to accumulated phase and generate sample
        phase += phaseIncrement;
        double x = std::sin(phase);

        x = sim.process_sample(x);

        std::cout
            << t << ","
            << x
            << std::endl;
    }

    return 0;
}
