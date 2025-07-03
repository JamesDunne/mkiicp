cmake_minimum_required(VERSION 3.31)
project(mna)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3 -fno-strict-aliasing -ffast-math")

add_executable(mna main.cpp
        Triode.h
        MNASolver.h
        V1A_ToneStack.h
        MarkIICPreamp.h
        V1B_Stage.h
        Coupling1.h
        V3B_Stage.h
        V4A_Stage.h
        WaveFile.h
        MixerAndV2B.h
        V2A_OutputStage.h)

target_compile_definitions(mna PUBLIC USE_TRIODE_CACHE)
