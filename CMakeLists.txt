
cmake_minimum_required (VERSION 2.6)

project (NumericTester)

set(CXX_COMPILE_FLAGS "-g -std=c++14 -Wall -finline-functions -findirect-inlining")
set(LD_LIBRARIES "-lmpfr -lgtest")

set(CMAKE_CXX_FLAGS "${CXX_COMPILE_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${LD_LIBRARIES}")

add_executable(dptest dotprod.cpp numerictester.cpp)
add_executable(quadtest quad.cpp numerictester.cpp)
add_executable(tests test.cpp numerictester.cpp)
