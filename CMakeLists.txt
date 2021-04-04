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
	${CMAKE_CURRENT_SOURCE_DIR}/include/ballistic/celestrack
	${CMAKE_CURRENT_SOURCE_DIR}/include/ballistic
	${CMAKE_CURRENT_SOURCE_DIR}/libs/General/include
)

add_library(
	Ballistic 

	src/Structures.cpp
	src/Ballistic.cpp

	src/GeoPotential.cpp

	src/Atmosphere1981.cpp
	src/Atmosphere2004.cpp

	src/EGM96.cpp
	src/JGM3.cpp 
	src/PZ90.cpp
	
	src/SolarModel.cpp

	src/Parsers.cpp
	src/Object.cpp

	src/AdjustTrajectory.cpp

	src/DataProviders.cpp
	src/Conversions.cpp)

add_subdirectory(
	${CMAKE_CURRENT_SOURCE_DIR}/libs/General
)

target_link_libraries(Ballistic	General)


add_executable(TestConversions test/TestConversions.cpp)
target_link_libraries(TestConversions Ballistic)

add_executable(TestForecast test/TestForecast.cpp)
target_link_libraries(TestForecast Ballistic)

add_executable(TestAtmosphere test/TestAtmosphere.cpp)
target_link_libraries(TestAtmosphere Ballistic)