cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(Tasks VERSION 0.1.0)
include_directories(
	include
	${CMAKE_SOURCE_DIR}/General/include)
add_library(
	Tasks
	src/AstroValues.cpp
	src/GeoPotential.cpp
)
add_executable(
	TestTasks
	test/Test.cpp
)
target_link_libraries(
	Tasks
	General
)
target_link_libraries(
	TestTasks
	Tasks
)