cmake_minimum_required (VERSION 3.8)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(Ballistic)

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/include/ballistic/atmosphere
	${CMAKE_CURRENT_SOURCE_DIR}/include/ballistic/gravity
	${CMAKE_CURRENT_SOURCE_DIR}/include/ballistic/solar
	${CMAKE_CURRENT_SOURCE_DIR}/include/ballistic/integrators
	${CMAKE_CURRENT_SOURCE_DIR}/include/ballistic/models
	${CMAKE_CURRENT_SOURCE_DIR}/include/ballistic/utils3d
	${CMAKE_CURRENT_SOURCE_DIR}/include/ballistic
	${CMAKE_CURRENT_SOURCE_DIR}/libs/General/include
)

add_library(
	Ballistic 

	src/Structures.cpp
	src/Ballistic.cpp

	src/GeoPotential.cpp

	src/StaticAtmosphere.cpp

	src/EGM96.cpp
	src/JGM3.cpp 
	src/PZ90.cpp
	
	
	src/SolarModel.cpp

	src/Parsers.cpp
	src/Object.cpp
	"src/TrajectoryProxy.cpp"
 )

add_subdirectory(
	${CMAKE_CURRENT_SOURCE_DIR}/libs/General
)

target_link_libraries(
	Ballistic
	General
)

add_executable(
	TestBallistic
	test/Test.cpp
)

add_executable(
	TestAdd
	test/add.cpp
)

target_link_libraries(
	TestBallistic
	Ballistic
)

