// triode_preamp_sim.cpp
// Full nonlinear triode simulation with full analog tone stack emulation for PCM processing

#include <cmath>
#include <algorithm>
#include <complex>
#include <iostream>

struct RCFilter {
    explicit RCFilter(float r, float c) : R(r), C(c) {}

    float R, C;
    float v = 0.0f;
    float sampleRate = 48000.0f;

    float process(float input) {
        float alpha = 1.0f / (1.0f + R * C * sampleRate);
        v = alpha * input + (1.0f - alpha) * v;
        return v;
    }
};

struct Triode12AX7 {
    Triode12AX7(float r_p, float r_k, float c_k, float b_plus, float fs = 48000.0f)
      : B_plus(b_plus),
        R_p(r_p),
        R_k(r_k),
        C_k(c_k),
        sampleRate(fs)
    {}

    constexpr static float MU = 96.2f;
    constexpr static float EX = 1.437f;
    constexpr static float KG1 = 613.4f;
    constexpr static float KP = 740.3f;
    constexpr static float KVB = 1672.0f;

    float B_plus = 300.0f;
    float R_p = 100000.0f;
    float R_k = 1500.0f;
    float C_k = 0.0f; // Cathode bypass capacitor in Farads
    float sampleRate = 48000.0f;

    float V_k = 1.0f;
    float V_ck = 0.0f; // Capacitor voltage (for bypass dynamics)

    constexpr static int max_iter = 20;
    constexpr static float tol = 1e-4f;

    float process(float V_g) {
        float V_a = B_plus;
        float V_gk = V_g - V_k;

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
            float dv = i_c / C_k / sampleRate;
            V_ck += dv;
            V_k = V_ck + current * R_k;
        } else {
            V_k = current * R_k;
        }

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
        float Ia1 = Ia(V_ak + h, V_gk);
        float Ia0 = Ia(V_ak - h, V_gk);
        return (Ia1 - Ia0) / (2.0f * h);
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
    Triode12AX7 stage1, stage2, stage3, stage4, stage5, stage6;
    FullToneStack tone;
    RCFilter inputGridFilter1, inputGridFilter2, gridDropFilter3;
    RCFilter gainPotFilter;
    float gain = 0.5f;
    float master = 0.5f;
    float volume1 = 0.75f;

    PreampSimulator() :
        stage1(150000.0f, 1500.0f, 22.47e-6f, 405.0f),
        stage2(100000.0f, 1500.0f, 22e-6f,    405.0f),
        stage3(82000.0f,  1500.0f,  2.2e-6f,  410.0f),
        stage4(270000.0f, 3300.0f,  0.22e-6f, 410.0f),
        stage5(47000.0f,  1000.0f, 15e-6f,    410.0f),
        stage6(47000.0f,  2200.0f,  0.0f,     410.0f)
    {
        inputGridFilter1.R = 1000000.0f; inputGridFilter1.C = 20e-12f;
        inputGridFilter2.R = 3300000.0f; inputGridFilter2.C = 20e-12f;
        gridDropFilter3.R = 220000.0f;   gridDropFilter3.C = 250e-12f;
        gainPotFilter.R = 100000.0f * gain + 100000.0f; gainPotFilter.C = 180e-12f;
    }

    float process(float input) {
        float x = input;
        x = inputGridFilter1.process(x);
        x = inputGridFilter2.process(x);
        x = stage1.process(x);

        float tone_in = x;
        static float v_bass = 0.0f;
        static float v_mid = 0.0f;
        static float v_treble = 0.0f;

        float Rb = tone.R2 * (1.0f - tone.bass);
        float Rm = tone.R3 * (1.0f - tone.mid);
        float Rt = tone.R1 * (1.0f - tone.treble);

        float alpha_bass = 1.0f / (1.0f + Rb * tone.C2 * 48000.0f);
        float alpha_mid  = 1.0f / (1.0f + Rm * tone.C3 * 48000.0f);
        float alpha_treb = 1.0f / (1.0f + Rt * tone.C1 * 48000.0f);

        v_bass = alpha_bass * tone_in + (1.0f - alpha_bass) * v_bass;
        v_mid  = alpha_mid  * tone_in + (1.0f - alpha_mid)  * v_mid;
        v_treble = alpha_treb * tone_in + (1.0f - alpha_treb) * v_treble;

        float tone_out = (v_bass + v_mid + v_treble) / 3.0f;

        float v_after_vol = tone_out * volume1;
        float v_lp = gainPotFilter.process(v_after_vol);

        x = stage2.process(v_lp);
        x = stage3.process(x);
        x = stage4.process(gridDropFilter3.process(x));
        x = stage5.process(x);
        x = stage6.process(x);
        return x * master;
    }
};

int main() {
    // FullToneStack tone;
    // tone.treble = 0.8f;
    // tone.mid = 0.33f;
    // tone.bass = 0.05f;
    // tone.plotFrequencyResponse();

    float sampleRate = 48000.0f;    // Audio sample rate in Hz

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
    for (long i = 0; i < 48000 * 5; ++i) {
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
