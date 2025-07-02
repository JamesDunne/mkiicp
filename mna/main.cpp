#include <iostream>

#include "MarkIICPreamp.h"

int main() {
    MarkIICPreamp preamp;
    preamp.setup(48000.0);

    std::cout << "Hello, World!" << std::endl;

    for (int i = 0; i < 10; i++) {
        double x = preamp.processSample(0.0);
        std::cout << x << std::endl;
    }

    return 0;
}