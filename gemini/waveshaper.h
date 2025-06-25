#include <iostream>
#include <vector>
#include <iomanip> // For std::fixed and std::setprecision

/**
 * @class TriodeStageWaveShaper
 * @brief Simulates the static non-linear gain of a 12AX7 common-cathode stage.
 *
 * This class implements a 5th-order polynomial approximation of the
 * Vout vs. Vin transfer curve derived from a SPICE DC simulation.
 * It models the instantaneous (memoryless) behavior of the tube stage.
 *
 * Circuit Parameters Used for Simulation:
 * - Tube: 12AX7
 * - Plate Resistor: 150k Ohm
 * - Cathode Resistor: 1.5k Ohm
 * - B+ Voltage: 405V
 * - Input Signal Range: -1.0V to +1.0V centered around 0V bias.
 */
class TriodeStageWaveShaper {
public:
    TriodeStageWaveShaper() {
        // Polynomial coefficients for:
        // Vout = p5*Vin^5 + p4*Vin^4 + p3*Vin^3 + p2*Vin^2 + p1*Vin + p0
        // Derived from SPICE .DC sweep from Vin = -1.0V to +1.0V
        p5 =  2.990;
        p4 = -2.964;
        p3 = -11.97;
        p2 =  12.04;
        p1 = -76.10;
        p0 =  236.2; // DC Bias offset
    }

    /**
     * @brief Processes a single input sample (voltage) to produce an output sample.
     * @param vin The input voltage sample. Should be in the range [-1.0, 1.0].
     * @return The resulting output voltage from the simulated tube stage.
     */
    double processSample(double vin) const {
        // To accurately model the gain structure, we first get the non-linear output
        // which includes the large DC offset. This is the "absolute" plate voltage.
        // We use Horner's method for efficient polynomial evaluation:
        // p0 + x*(p1 + x*(p2 + x*(p3 + x*(p4 + x*p5))))
        double v_plate = p0 + vin * (p1 + vin * (p2 + vin * (p3 + vin * (p4 + vin * p5))));

        // In a real audio plugin, you would likely want to remove the DC offset
        // to get an audio signal centered around 0. This is done by subtracting
        // the quiescent voltage (the output when input is 0).
        // double v_audio_out = v_plate - p0;
        //
        // However, the request is to simulate the plate output, so we return the
        // full plate voltage.
        return v_plate;
    }
    
    /**
     * @brief Returns the DC bias voltage (output when input is 0V).
     */
    double getDCOffset() const {
        return p0;
    }

private:
    // Polynomial coefficients
    double p5, p4, p3, p2, p1, p0;
};
