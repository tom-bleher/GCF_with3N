cmake_minimum_required(VERSION 3.6)
project(genQE_FSI_project)

# Default to an optimized build: with no CMAKE_BUILD_TYPE, CMake compiles
# without any optimization flags (-O0).
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../../cmake_modules")
find_package(ROOT REQUIRED)

set(CMAKE_CXX_STANDARD 17)

# Single switch: ON -> FSI enabled via GENIE, OFF -> no FSI at all
option(USE_FSI "Enable FSI via GENIE intranuclear cascade" ON)
set(GENIE_SOURCE_DIR "${CMAKE_SOURCE_DIR}/Generator-R-3_06_02" CACHE PATH
    "Path to built GENIE source tree (contains src/ and lib/)")

set(SRC_ROOT "${CMAKE_SOURCE_DIR}/../..")
set(LIBRARIES "${SRC_ROOT}/libraries")

set(COMMON_INCLUDE_DIRS
    ${ROOT_INCLUDE_DIR}
    ${SRC_ROOT}/include
    ${LIBRARIES}
    ${CMAKE_SOURCE_DIR}
)
include_directories(${COMMON_INCLUDE_DIRS})

set(GCF_LIBS
    ${LIBRARIES}/cross_sections/eNCrossSection.cc
    ${LIBRARIES}/generator/gcfGenerator.cc
    ${LIBRARIES}/nucleus/gcfNucleus.cc
    ${LIBRARIES}/nucleus/gcfSRC.cc
)

# Shared FSI helpers (compiled only when USE_FSI is on; needs GENIE headers).
set(FSI_LIBS
    ${LIBRARIES}/fsi/GenieFSIHelpers.cc
)

if(USE_FSI)
    add_compile_definitions(USE_FSI)
    add_compile_definitions(GENIE_OVERRIDE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/config/genie_override")

    if(NOT EXISTS "${GENIE_SOURCE_DIR}/src/Framework/GHEP/GHepRecord.h")
        message(FATAL_ERROR "USE_FSI=ON but GENIE headers not found in ${GENIE_SOURCE_DIR}/src")
    endif()
    if(NOT EXISTS "${GENIE_SOURCE_DIR}/lib")
        message(FATAL_ERROR "USE_FSI=ON but GENIE libs directory not found: ${GENIE_SOURCE_DIR}/lib")
    endif()

    include_directories("${GENIE_SOURCE_DIR}/src")

    find_path(LOG4CPP_INCLUDE_DIR
        NAMES log4cpp/Category.hh
        HINTS /opt/homebrew/include /usr/local/include)
    find_library(LOG4CPP_LIBRARY
        NAMES log4cpp
        HINTS /opt/homebrew/lib /usr/local/lib)

    find_path(LIBXML2_INCLUDE_DIR
        NAMES libxml/tree.h
        HINTS
            /opt/homebrew/opt/libxml2/include/libxml2
            /opt/homebrew/include/libxml2
            /usr/local/opt/libxml2/include/libxml2
            /usr/include/libxml2)
    find_library(LIBXML2_LIBRARY
        NAMES xml2
        HINTS
            /opt/homebrew/opt/libxml2/lib
            /opt/homebrew/lib
            /usr/local/opt/libxml2/lib
            /usr/local/lib)

    if(NOT LOG4CPP_INCLUDE_DIR OR NOT LOG4CPP_LIBRARY)
        message(FATAL_ERROR "USE_FSI=ON but log4cpp headers/libs not found. Install log4cpp (e.g. brew install log4cpp)")
    endif()
    if(NOT LIBXML2_INCLUDE_DIR OR NOT LIBXML2_LIBRARY)
        message(FATAL_ERROR "USE_FSI=ON but libxml2 headers/libs not found. Install libxml2 (e.g. brew install libxml2)")
    endif()

    include_directories(${LOG4CPP_INCLUDE_DIR} ${LIBXML2_INCLUDE_DIR})

    file(GLOB GENIE_SHARED_LIBS "${GENIE_SOURCE_DIR}/lib/libG*.dylib")
    if(NOT GENIE_SHARED_LIBS)
        file(GLOB GENIE_SHARED_LIBS "${GENIE_SOURCE_DIR}/lib/libG*.so")
    endif()
    if(NOT GENIE_SHARED_LIBS)
        message(FATAL_ERROR "No GENIE shared libs found in ${GENIE_SOURCE_DIR}/lib")
    endif()

    # Embed GENIE lib path in RPATH so DYLD_LIBRARY_PATH is not needed at runtime
    set(CMAKE_BUILD_RPATH "${GENIE_SOURCE_DIR}/lib")
endif()

add_executable(genQE
    genQE.cc
    QEGenerator.cc
    ${GCF_LIBS}
)
target_link_libraries(genQE ${ROOT_LIBRARIES})

add_executable(genQE_FSI
    genQE_FSI.cc
    QEGenerator.cc
    QEGeneratorFSI.cc
    ${GCF_LIBS}
    ${FSI_LIBS}
)
target_link_libraries(genQE_FSI ${ROOT_LIBRARIES})

add_executable(SRC_analysis_2N
    analysis/SRC_analysis_2N.cpp
    QEGenerator.cc
    QEGeneratorFSI.cc
    ${GCF_LIBS}
    ${FSI_LIBS}
)
target_link_libraries(SRC_analysis_2N ${ROOT_LIBRARIES})

if(USE_FSI)
    target_link_libraries(genQE_FSI ${GENIE_SHARED_LIBS} ${LOG4CPP_LIBRARY} ${LIBXML2_LIBRARY})
    target_link_libraries(SRC_analysis_2N ${GENIE_SHARED_LIBS} ${LOG4CPP_LIBRARY} ${LIBXML2_LIBRARY})
endif()
