# CMakeLists for GCF_Generator_Suite

# Set up
cmake_minimum_required (VERSION 3.6)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)
# genCCQE and genPhoto call TH1::GetRandom(TRandom*), added in ROOT 6.24.
# FindROOT.cmake ignores a version argument, so check root-config's here.
string(REPLACE "/" "." ROOT_VERSION_DOTTED "${ROOTVERSION}")
if(ROOT_VERSION_DOTTED VERSION_LESS 6.24)
    message(FATAL_ERROR "ROOT ${ROOTVERSION} found; ROOT >= 6.24 is required")
endif()

project (GCF_Generator_Suite)

# Default to an optimized build: with no CMAKE_BUILD_TYPE, CMake compiles
# without any optimization flags (-O0).
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

if(NOT CMAKE_CXX_FLAGS)
       set(CMAKE_CXX_FLAGS "-std=c++17")
endif()

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include)
set(LIBRARIES ${CMAKE_SOURCE_DIR}/libraries)

include_directories(${INCLUDES} ${LIBRARIES})

add_subdirectory(programs)