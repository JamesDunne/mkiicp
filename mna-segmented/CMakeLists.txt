cmake_minimum_required(VERSION 3.31)
project(mna_segmented)

set(CMAKE_CXX_STANDARD 20)

add_executable(mna_segmented main.cpp
        Matrix.h
        MNASolver.h
        InputAndToneStage.h
        InputAndToneStage.cpp
        TubeParams.h
        TubeStageV1B.h
        TubeStageV1B.cpp
        Interstage_C7_R9.h
        Interstage_C7_R9.cpp
        LeadDriveNetwork.h
        LeadDriveNetwork.cpp
        MainMixer.h
        MainMixer.cpp
        TubeModel.h
        TubeStageV3B.h
        TubeStageV3B.cpp
        LeadChannelFinisher.h
        LeadChannelFinisher.cpp
        TubeStageV2B.h
        TubeStageV2B.cpp
        FinalStageAndMaster.cpp
        FinalStageAndMaster.h
        MkiicPreamp.cpp
        MkiicPreamp.h)
