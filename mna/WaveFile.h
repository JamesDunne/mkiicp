#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <functional>  // For std::function
#include <cstdint>     // For fixed-width integers
#include <stdexcept>
#include <algorithm>   // For std::min/max

// Use pragma pack to ensure the struct is packed without padding.
// This is essential for reading binary data directly into the struct.
#pragma pack(push, 1)
struct WavHeader {
    // RIFF Chunk
    char riff_id[4];         // "RIFF"
    uint32_t chunk_size;     // File size - 8
    char wave_id[4];         // "WAVE"
    // Fmt Sub-Chunk
    char fmt_id[4];          // "fmt "
    uint32_t fmt_size;       // Size of this sub-chunk (16 for PCM)
    uint16_t audio_format;   // 1 for PCM
    uint16_t num_channels;   // 1=mono, 2=stereo
    uint32_t sample_rate;
    uint32_t byte_rate;      // sample_rate * num_channels * bits_per_sample/8
    uint16_t block_align;    // num_channels * bits_per_sample/8
    uint16_t bits_per_sample;
};
#pragma pack(pop)

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
    const std::function<double(double, long)>& processor
) {
    // --- Open Input File ---
    std::ifstream ifs(inputFilename, std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("Could not open input file: " + inputFilename);
    }

    // --- Read and Validate Header ---
    WavHeader header;
    ifs.read(reinterpret_cast<char*>(&header), sizeof(header));

    if (std::string(header.riff_id, 4) != "RIFF" || std::string(header.wave_id, 4) != "WAVE") {
        throw std::runtime_error("Invalid RIFF/WAVE file");
    }
    if (std::string(header.fmt_id, 4) != "fmt ") {
         throw std::runtime_error("Missing 'fmt ' chunk");
    }
    if (header.audio_format != 1) { // 1 = PCM
        throw std::runtime_error("Unsupported audio format (only PCM is supported)");
    }

    // --- Find Data Chunk ---
    char data_id[4];
    uint32_t data_size;
    while (ifs.read(data_id, 4) && ifs.read(reinterpret_cast<char*>(&data_size), 4)) {
        if (std::string(data_id, 4) == "data") {
            break; // Found it
        }
        // Skip unknown chunks by seeking past them
        ifs.seekg(data_size, std::ios::cur);
    }
    if (std::string(data_id, 4) != "data") {
        throw std::runtime_error("Missing 'data' chunk");
    }

    // --- Read Sample Data ---
    std::vector<char> raw_data(data_size);
    ifs.read(raw_data.data(), data_size);
    ifs.close();

    // --- Process Samples ---
    const long num_samples_total = data_size / (header.bits_per_sample / 8);
    const int bytes_per_sample = header.bits_per_sample / 8;
    std::vector<char> processed_raw_data(data_size);

    for (long i = 0; i < num_samples_total; ++i) {
        const long sample_offset = i * bytes_per_sample;
        double sample = 0.0;

        // 1. Extract sample as double in range [-1.0, 1.0]
        switch (header.bits_per_sample) {
            case 8: // Unsigned 8-bit
                sample = (reinterpret_cast<uint8_t*>(raw_data.data())[i] - 128) / 128.0;
                break;
            case 16: // Signed 16-bit
                sample = reinterpret_cast<int16_t*>(raw_data.data())[i] / 32768.0;
                break;
            case 24: { // Signed 24-bit
                const uint8_t* p = reinterpret_cast<uint8_t*>(raw_data.data()) + sample_offset;
                int32_t val = (p[2] << 24 | p[1] << 16 | p[0] << 8) >> 8; // Sign extend
                sample = val / 8388608.0;
                break;
            }
            case 32: // Signed 32-bit
                sample = reinterpret_cast<int32_t*>(raw_data.data())[i] / 2147483648.0;
                break;
            default: throw std::runtime_error("Unsupported bits per sample");
        }

        // 2. Process the sample using the provided function
        double processed_sample = processor(sample, i);

        // Clamp to prevent audio artifacts from out-of-range values
        processed_sample = std::max(-1.0, std::min(1.0, processed_sample));

        // 3. Write processed sample back to a new raw data buffer
        switch (header.bits_per_sample) {
            case 8:
                reinterpret_cast<uint8_t*>(processed_raw_data.data())[i] = static_cast<uint8_t>(processed_sample * 127.0 + 128.0);
                break;
            case 16:
                reinterpret_cast<int16_t*>(processed_raw_data.data())[i] = static_cast<int16_t>(processed_sample * 32767.0);
                break;
            case 24: {
                int32_t val = static_cast<int32_t>(processed_sample * 8388607.0);
                uint8_t* p = reinterpret_cast<uint8_t*>(processed_raw_data.data()) + sample_offset;
                p[0] = val & 0xFF;
                p[1] = (val >> 8) & 0xFF;
                p[2] = (val >> 16) & 0xFF;
                break;
            }
            case 32:
                reinterpret_cast<int32_t*>(processed_raw_data.data())[i] = static_cast<int32_t>(processed_sample * 2147483647.0);
                break;
        }
    }

    // --- Write Output File ---
    std::ofstream ofs(outputFilename, std::ios::binary);
    if (!ofs) {
        throw std::runtime_error("Could not create output file: " + outputFilename);
    }

    ofs.write(reinterpret_cast<const char*>(&header), sizeof(header));
    ofs.write(data_id, 4);
    ofs.write(reinterpret_cast<const char*>(&data_size), 4);
    ofs.write(processed_raw_data.data(), data_size);

    ofs.close();
}
