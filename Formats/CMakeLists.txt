cmake_minimum_required(VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(Formats VERSION 0.1.0)
include_directories(
	${PROJECT_SOURCE_DIR}/include
	${CMAKE_SOURCE_DIR}/General/include
	${CMAKE_SOURCE_DIR}/Tasks/include
)
add_library(
	Formats
	src/TleFormat.cpp
)
target_link_libraries(
	Formats
	Tasks
)
add_executable(
	TestFormats
	test/Test.cpp
)
target_link_libraries(
	TestFormats
	Formats
)