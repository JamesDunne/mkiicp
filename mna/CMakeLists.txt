cmake_minimum_required(VERSION 3.31)
project(mna)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -O3 -fno-strict-aliasing -ffast-math")

add_executable(mna main.cpp
        WaveFile.h
        Triode.h
        MNASolver.h
        MarkIICPreamp.h
        V1A_ToneStack.h
        V1B_Stage.h
        LeadAndMixerStage.h
        V2B_Stage.h
        V2A_OutputStage.h)

target_compile_definitions(mna PUBLIC USE_TRIODE_CACHE)
