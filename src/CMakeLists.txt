# CMakeLists for GCF_Generator_Suite

# Set up
cmake_minimum_required (VERSION 3.6)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)

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