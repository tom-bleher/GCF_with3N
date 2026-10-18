# Pipeline stages: gcf_generate (PWIA event generation) and gcf_fsi (GENIE
# FSI on an event file). See docs/PIPELINE.md.
#
#   cmake -S src/programs/gcf_stages -B src/programs/gcf_stages/build \
#         -DLOG4CPP_INCLUDE_DIR=... -DLOG4CPP_LIBRARY=...
#   cmake --build src/programs/gcf_stages/build
#
# Both tools need the vendored GENIE built in place (install.sh does that):
# gcf_fsi runs the cascade, and gcf_generate draws the interaction vertex
# from GENIE's nuclear density.

cmake_minimum_required(VERSION 3.6)
project(gcf_stages)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../../cmake_modules")
find_package(ROOT REQUIRED)
set(CMAKE_CXX_STANDARD 17)

get_filename_component(SRC_ROOT "${CMAKE_SOURCE_DIR}/../.." REALPATH)
set(LIBRARIES "${SRC_ROOT}/libraries")
set(GEN3N_DIR "${SRC_ROOT}/programs/genQE_3N")
set(GEN2N_DIR "${SRC_ROOT}/programs/genQE_FSI")
set(GENIE_SOURCE_DIR "${GEN2N_DIR}/Generator-R-3_06_02" CACHE PATH
    "Path to built GENIE source tree (contains src/ and lib/)")

if(NOT EXISTS "${GENIE_SOURCE_DIR}/src/Framework/GHEP/GHepRecord.h")
    message(FATAL_ERROR "GENIE headers not found in ${GENIE_SOURCE_DIR}/src")
endif()
file(GLOB GENIE_SHARED_LIBS "${GENIE_SOURCE_DIR}/lib/libG*.so")
if(NOT GENIE_SHARED_LIBS)
    file(GLOB GENIE_SHARED_LIBS "${GENIE_SOURCE_DIR}/lib/libG*.dylib")
endif()
if(NOT GENIE_SHARED_LIBS)
    message(FATAL_ERROR "No GENIE shared libs in ${GENIE_SOURCE_DIR}/lib; build GENIE first (install.sh)")
endif()

find_path(LOG4CPP_INCLUDE_DIR NAMES log4cpp/Category.hh
          HINTS $ENV{HOME}/local/include /opt/homebrew/include /usr/local/include)
find_library(LOG4CPP_LIBRARY NAMES log4cpp
             HINTS $ENV{HOME}/local/lib /opt/homebrew/lib /usr/local/lib)
find_path(LIBXML2_INCLUDE_DIR NAMES libxml/tree.h
          HINTS /usr/include/libxml2 /opt/homebrew/opt/libxml2/include/libxml2
                /usr/local/include/libxml2)
find_library(LIBXML2_LIBRARY NAMES xml2
             HINTS /opt/homebrew/opt/libxml2/lib /usr/local/lib)
if(NOT LOG4CPP_INCLUDE_DIR OR NOT LOG4CPP_LIBRARY)
    message(FATAL_ERROR "log4cpp not found; pass -DLOG4CPP_INCLUDE_DIR=... -DLOG4CPP_LIBRARY=...")
endif()
if(NOT LIBXML2_INCLUDE_DIR OR NOT LIBXML2_LIBRARY)
    message(FATAL_ERROR "libxml2 not found; pass -DLIBXML2_INCLUDE_DIR=... -DLIBXML2_LIBRARY=...")
endif()

# Code version recorded in every output file.
execute_process(COMMAND git -C "${SRC_ROOT}" describe --always --dirty --tags
                OUTPUT_VARIABLE GCF_GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET)
if(NOT GCF_GIT_VERSION)
    set(GCF_GIT_VERSION "unknown")
endif()

add_compile_definitions(
    USE_FSI
    GENIE_OVERRIDE_DIR="${GEN2N_DIR}/config/genie_override"
    GCF_3N_SOURCE_DIR="${GEN3N_DIR}"
    GCF_GIT_VERSION="${GCF_GIT_VERSION}")
include_directories(${ROOT_INCLUDE_DIR} ${SRC_ROOT}/include ${LIBRARIES}
                    ${GEN3N_DIR} ${GEN2N_DIR} "${GENIE_SOURCE_DIR}/src"
                    ${LOG4CPP_INCLUDE_DIR} ${LIBXML2_INCLUDE_DIR})
set(CMAKE_BUILD_RPATH "${GENIE_SOURCE_DIR}/lib")

set(GCF_LIBS
    ${LIBRARIES}/cross_sections/eNCrossSection.cc
    ${LIBRARIES}/generator/gcfGenerator.cc
    ${LIBRARIES}/nucleus/gcfNucleus.cc
    ${LIBRARIES}/nucleus/gcfSRC.cc
    ${LIBRARIES}/fsi/GenieFSIHelpers.cc)
set(DEP_LIBS ${ROOT_LIBRARIES} ${GENIE_SHARED_LIBS} ${LOG4CPP_LIBRARY} ${LIBXML2_LIBRARY})

add_executable(gcf_generate gcf_generate.cc
    ${GEN3N_DIR}/QEGenerator_3N.cc ${GEN2N_DIR}/QEGeneratorFSI.cc ${GCF_LIBS})
target_link_libraries(gcf_generate ${DEP_LIBS})

add_executable(gcf_fsi gcf_fsi.cc ${LIBRARIES}/fsi/FSIStage.cc ${GCF_LIBS})
target_link_libraries(gcf_fsi ${DEP_LIBS})
