cmake_minimum_required(VERSION 3.31)
project(preamp)

set(CMAKE_CXX_STANDARD 20)

add_executable(preamp main.cpp
        IIRBiquad.h
        Preamp.h
        Preamp.cpp
        wav.cpp
        ToneStack.cpp
        ToneStack.h
        IIRFilter.h
        TubeStage.cpp
        TubeStage.h
)

target_compile_definitions(preamp PUBLIC USE_TONESTACK_MNA)
