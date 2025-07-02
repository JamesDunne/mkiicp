#pragma once
#include <functional>
#include <string>

void processWavFile(
    const std::string& inputFilename,
    const std::string& outputFilename,
    const std::function<double(double)>& processor
);
