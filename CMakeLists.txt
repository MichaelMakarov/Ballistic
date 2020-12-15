cmake_minimum_required (VERSION 3.8)

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
	src/PlahovIntegrator.cpp)

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