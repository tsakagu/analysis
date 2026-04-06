cmake_minimum_required(VERSION 3.10)
project(CustomHistogramFit)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find ROOT package
find_package(ROOT REQUIRED)

# Include ROOT directories
include(${ROOT_USE_FILE})

# Source files
file(GLOB SOURCES src/*.cc)

# Include directories
include_directories(include)

# Create shared library
add_library(CustomHistogramFit SHARED ${SOURCES})

install(DIRECTORY include/CustomHistogramFit DESTINATION include)
install(TARGETS CustomHistogramFit DESTINATION lib)
