cmake_minimum_required(VERSION 3.31)
project(mna)

set(CMAKE_CXX_STANDARD 20)

add_executable(mna main.cpp
        Triode.h
        MNASolver.h
        V1A_ToneStack.h
        MarkIICPreamp.h
        V1B_Stage.h
        Coupling1.h
        V3B_Stage.h
        LeadPathCoupling2.h
        V4A_Stage.h
        WaveFile.h)
