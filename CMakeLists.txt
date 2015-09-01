
cmake_minimum_required (VERSION 2.6)

project (NumericTester)

set(CXX_COMPILE_FLAGS "-g -std=c++14")
set(LD_LIBRARIES "-lmpfr")

set(CMAKE_CXX_FLAGS "${CXX_COMPILE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${LD_LIBRARIES}")

add_executable(numerictester dotprod.cpp)