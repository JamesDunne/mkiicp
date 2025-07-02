#include <iostream>

#include "MarkIICPreamp.h"

int main() {
    MarkIICPreamp preamp;
    preamp.setup(48000.0);

    std::cout << "Hello, World!" << std::endl;

    double x = 0.0;
    for (int i = 0; i < 100000; i++) {
        x = preamp.processSample(0.0);
    }

    std::cout << x << std::endl;
    for (int i = 0; i < 10; i++) {
        x = preamp.processSample(0.0);
        std::cout << x << std::endl;
    }

    return 0;
}