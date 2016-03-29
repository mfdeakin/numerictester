
cmake_minimum_required (VERSION 2.6)

project (NumericTester)

set(CXX_COMPILE_FLAGS "-O3 -std=c++14 -Wall --save-temps")

set(CMAKE_CXX_FLAGS "${CXX_COMPILE_FLAGS}")

add_executable(dptest dotprod.cpp numerictester.cpp)
add_executable(quadtest quad.cpp numerictester.cpp)
add_executable(tests test.cpp numerictester.cpp)

target_link_libraries(dptest mpfr)
target_link_libraries(quadtest mpfr)
target_link_libraries(tests gtest mpfr pthread)
