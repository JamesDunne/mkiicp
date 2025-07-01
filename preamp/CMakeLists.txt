cmake_minimum_required(VERSION 3.31)
project(preamp)

set(CMAKE_CXX_STANDARD 20)

add_executable(preamp main.cpp
        IIRBiquad.h
        Preamp.h
        Preamp.cpp)
