
#include <iostream>

#include "sim.hpp"

int main(int argc, char *argv[]) {
    const double sampleRate = 48000.0;
    RealtimeTubeSim sim { sampleRate };

    // Define nodes
    sim.add_node("in"); sim.add_node("g"); sim.add_node("p"); sim.add_node("k");
    sim.add_node("V_P"); sim.add_node("ts1"); sim.add_node("ts2"); sim.add_node("ts3");
    sim.add_node("out"); sim.add_node("out_final");

    sim.add_resistor("in", "gnd", 1e6);
    sim.add_resistor("g", "gnd", 1e6);
    sim.add_resistor("V_P", "p", 100e3);
    sim.add_resistor("k", "gnd", 1.5e3);
    sim.add_resistor("out_final", "gnd", 1e6);
    sim.add_capacitor("in", "g", 22e-9);
    sim.add_capacitor("k", "gnd", 22e-6);
    sim.add_capacitor("p", "ts1", 22e-9);
    sim.add_capacitor("ts1", "ts3", 250e-12);
    sim.add_capacitor("ts2", "ts3", 22e-9);
    sim.add_capacitor("out", "out_final", 100e-9);

    sim.add_voltage_source("V_P", "gnd", 300.0, false);
    sim.add_voltage_source("in", "gnd", 0.0, true);

    m_pot_treble = {m_node("ts1"), m_node("ts2"), -1, 250e3};
    m_pot_bass = {m_node("ts2"), m_node("gnd"), -1, 1e6};
    m_pot_mid = {m_node("ts3"), m_node("gnd"), -1, 25e3};
    m_pot_master = {m_node("ts3"), m_node("gnd"), m_node("out"), 1e6};

    // Add components
    sim.add_resistor("in", "gnd", 1e6);
    sim.add_resistor("g", "gnd", 1e6);
    // ... and so on for all R and C components ...

    // Add the triode with specified model parameters
    sim.add_triode("p", "g", "k", 96.20, 1.437, 613.4, 740.3, 1672.0, 2000.0);

    // Add sources
    sim.add_voltage_source("V_P", "gnd", 300.0, false);
    sim.add_voltage_source("in", "gnd", 0.0, true); // Input signal

    sim.prepare_to_play(); // Finalizes matrices and solves DC point

    // // Set tone controls
    // sim.set_volume(0.75);
    // sim.set_treble(0.8);
    // sim.set_mid(0.33);
    // sim.set_bass(0.05);
    // sim.set_master(0.5);

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
