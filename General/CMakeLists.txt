cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(General VERSION 0.1.0)
include_directories(include)
add_library(
    General 
    src/XYZ.cpp
    src/RBL.cpp
    src/DateTime.cpp
    src/StateTypes.cpp
    src/MathFunctions.cpp 
    src/LegendrePolynom.cpp  )
add_executable(
    TestGeneral
    test/Test.cpp)
target_link_libraries(
    TestGeneral
    General)