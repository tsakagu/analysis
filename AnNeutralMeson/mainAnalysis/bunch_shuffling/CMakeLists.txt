cmake_minimum_required(VERSION 3.10)
project(AsymmetryCalc)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find ROOT package
find_package(ROOT REQUIRED COMPONENTS RIO)

# Include ROOT directories
include(${ROOT_USE_FILE})

# Source files
file(GLOB SOURCES src/*.cc)

# Include directories
include_directories(include)

# Create shared library
add_library(AsymmetryCalc SHARED ${SOURCES})

# Create executables
add_subdirectory(run)

install(DIRECTORY include/AsymmetryCalc DESTINATION include)
install(TARGETS AsymmetryCalc DESTINATION lib)
