// triode_preamp_sim.cpp
// Refactored nonlinear 12AX7 simulation with tone stack and full gain structuring from LTspice netlist

#include <cmath>
#include <algorithm>
#include <iostream>

constexpr float sampleRate = 48000.0f;

struct RCFilter {
    float R, C;
    float v = 0.0f;
    float process(float input) {
        float alpha = 1.0f / (1.0f + R * C * sampleRate);
        v = alpha * input + (1.0f - alpha) * v;
        return v;
    }
};

struct Triode12AX7 {
    constexpr static float MU = 96.2f, EX = 1.437f, KG1 = 613.4f, KP = 740.3f, KVB = 1672.0f;
    float B_plus = 300.0f, R_p = 100000.0f, R_k = 1500.0f, C_k = 0.0f;
    float V_k = 1.0f, V_ck = 0.0f;
    int max_iter = 20;
    float tol = 1e-4f;
    float process(float V_g) {
        float V_a = B_plus, V_gk = V_g - V_k;
        for (int i = 0; i < max_iter; ++i) {
            float f = (B_plus - V_a) / R_p - Ia(V_a, V_gk);
            float dfdVa = -1.0f / R_p - dIa_dVa(V_a, V_gk);
            float delta = f / dfdVa;
            V_a -= delta;
            if (std::abs(delta) < tol) break;
        }
        float current = Ia(V_a, V_gk);
        if (C_k > 0.0f) {
            float i_c = current - V_ck / R_k;
            V_ck += (i_c / C_k) / sampleRate;
            V_k = V_ck + current * R_k;
        } else V_k = current * R_k;
        return V_a;
    }
private:
    float Ia(float V_ak, float V_gk) const {
        float sqrt_arg = std::sqrt(KVB + V_ak * V_ak);
        float exponent = KP * (1.0f / MU + V_gk / sqrt_arg);
        float log_arg = 1.0f + std::exp(exponent);
        float inner = V_ak / KP * std::log(log_arg);
        float Ia = std::pow(inner, EX) / KG1;
        return std::max(0.0f, Ia);
    }
    float dIa_dVa(float V_ak, float V_gk) const {
        float h = 1e-3f;
        return (Ia(V_ak + h, V_gk) - Ia(V_ak - h, V_gk)) / (2.0f * h);
    }
};

struct FullToneStack {
    float treble = 0.8f, mid = 0.33f, bass = 0.05f;
    const float C1 = (750.0f + 250.0f) * 1e-12f;
    const float C2 = 0.1e-6f;
    const float C3 = 0.047e-6f;
    const float R1 = 250000.0f;
    const float R2 = 250000.0f;
    const float R3 = 10000.0f;
};

struct PreampSimulator {
    Triode12AX7 v1a, v1b, v3b, v4a, v2b, v2a;
    RCFilter inputR1C10, inputR10C10, c13bR5a, toneBass, toneMid, toneTreble;
    RCFilter c21r21, gainRloadC, c25r25, masterVolRC;
    FullToneStack tone;
    float gain = 0.5f, master = 0.5f, volume1 = 0.75f;

    PreampSimulator() {
        v1a.R_p = 150000.0f;  v1a.R_k = 1500.0f; v1a.C_k = 22e-6f;   v1a.B_plus = 405.0f;
        v1b.R_p = 100000.0f;  v1b.R_k = 1500.0f; v1b.C_k = 22e-6f;   v1b.B_plus = 405.0f;
        v3b.R_p = 82000.0f;   v3b.R_k = 1500.0f; v3b.C_k = 2.2e-6f;  v3b.B_plus = 410.0f;
        v4a.R_p = 274000.0f;  v4a.R_k = 3300.0f; v4a.C_k = 0.22e-6f; v4a.B_plus = 410.0f;
        v2b.R_p = 100000.0f;  v2b.R_k = 1500.0f; v2b.C_k = 0.0f;     v2b.B_plus = 410.0f;
        v2a.R_p = 120000.0f;  v2a.R_k = 1000.0f; v2a.C_k = 15.47f;   v2a.B_plus = 410.0f;

        inputR1C10.R = 1000000.0f; inputR1C10.C = 20e-12f;
        inputR10C10.R = 3300000.0f; inputR10C10.C = 20e-12f;
        c13bR5a.R = 100000.0f; c13bR5a.C = 180e-12f;

        toneBass.R = 250000.0f * (1 - tone.bass); toneBass.C = tone.C2;
        toneMid.R = 10000.0f * (1 - tone.mid); toneMid.C = tone.C3;
        toneTreble.R = (250000.0f * (1 - tone.treble)) + 100000.0f; toneTreble.C = tone.C1;

        c21r21.R = 680000.0f; c21r21.C = 0.02e-6f;
        gainRloadC.R = 475000.0f; gainRloadC.C = 120e-12f;
        c25r25.R = 270000.0f; c25r25.C = 0.022e-6f;
        masterVolRC.R = 1000000.0f * master; masterVolRC.C = 47e-12f;
    }

    float process(float input) {
        //float x = inputR1C10.process(input);
        float x = v1a.process(input);

        float toneIn = x;
        static float vb = 0.0f, vm = 0.0f, vt = 0.0f;
        float aB = 1.0f / (1.0f + toneBass.R * toneBass.C * sampleRate);
        float aM = 1.0f / (1.0f + toneMid.R * toneMid.C * sampleRate);
        float aT = 1.0f / (1.0f + toneTreble.R * toneTreble.C * sampleRate);
        vb = aB * toneIn + (1 - aB) * vb;
        vm = aM * toneIn + (1 - aM) * vm;
        vt = aT * toneIn + (1 - aT) * vt;
        float toneOut = (vb + vm + vt) / 3.0f * volume1;

        float x2 = v1b.process(c13bR5a.process(toneOut));
        float gainIn = c7r9.process(x2);

        // TODO
        //x = inputR10C10.process(gainIn);

        // float top = 1000000.0f * (1 - gain), bot = 1000000.0f * gain;
        // float gainOut = bot / (top + bot + 1e-9f) * gainIn;
        // x = v3b.process(gainRloadC.process(gainOut));
        // x = v4a.process(c25r25.process(x));
        //x = v2b.process(masterVolRC.process(x));
        return gainIn;
    }
};

int main() {
    PreampSimulator amplifier;
    amplifier.tone.treble = 0.8f;
    amplifier.tone.mid = 0.33f;
    amplifier.tone.bass = 0.05f;

    float startFreq = 20.0f;        // Start frequency in Hz
    float endFreq = 20000.0f;       // End frequency in Hz
    float duration = 5.0f;          // Sweep duration in seconds

    float k = (endFreq / startFreq);
    // float L = duration / std::log(k);

    double phase = 0.0;
    for (long i = 0; i < static_cast<long>(sampleRate * duration); ++i) {
        float t = static_cast<float>(i) / sampleRate;
        float instantaneousFreq = startFreq * std::pow(k, t / duration);

        // Calculate phase increment for this sample
        float phaseIncrement = 2.0f * static_cast<float>(M_PI) * instantaneousFreq / sampleRate;

        // Add to accumulated phase and generate sample
        phase += phaseIncrement;
        float x = std::sin(phase);

        x = amplifier.process(x);

        std::cout
            << t << ","
            << x
            << std::endl;
    }

    return 0;
}
