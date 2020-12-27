cmake_minimum_required (VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(Ballistic)
add_subdirectory(General)
add_subdirectory(Tasks)
add_subdirectory(Formats)

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${CMAKE_CURRENT_SOURCE_DIR}/General/include
	${CMAKE_CURRENT_SOURCE_DIR}/Tasks/include
)

add_library(
	Ballistic 
	src/Ballistic.cpp
	src/PlahovIntegrator.cpp
	src/RungeKuttaIntegrator.cpp
	src/AdamsBashforthIntegrator.cpp)

add_executable(
	TestBallistic
	test/Test.cpp
)

target_link_libraries(
	Ballistic
	General
	Tasks
)

target_link_libraries(
	TestBallistic
	Ballistic
)