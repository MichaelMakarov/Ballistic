cmake_minimum_required (VERSION 3.8)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(Ballistic)

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${CMAKE_CURRENT_SOURCE_DIR}/General/include
	${CMAKE_CURRENT_SOURCE_DIR}/Space/include
)

add_library(
	Ballistic 

	General/src/Geometry.cpp
	General/src/Times.cpp
	General/src/Mathematics.cpp
	General/src/Structures.cpp

	Space/src/GeoPotential.cpp
	Space/src/StaticAtmosphere.cpp

	src/Ballistic.cpp
)

add_executable(
	TestBallistic
	test/TestBallistic.cpp
)
add_executable(
	TestGeneralFunctions
	test/TestGeneral.cpp
)
add_executable(
	TestSpaceFunctions
	test/TestSpace.cpp
)

target_link_libraries(
	TestBallistic
	Ballistic
)
target_link_libraries(
	TestGeneralFunctions
	Ballistic
)
target_link_libraries(
	TestSpaceFunctions
	Ballistic
)