cmake_minimum_required(VERSION 3.31)
project(mkiicp)

set(CMAKE_CXX_STANDARD 20)

add_executable(mkiicp
        mkiicp.hpp
        main.cpp
        PreampSimulator.h)
