# genQE_3N — supports both standalone and suite (programs/) builds.
#
# Standalone:   cd src/programs/genQE_3N && cmake -S . -B build && cmake --build build
# Suite build:  cd src && cmake -S . -B build && cmake --build build
#
# In suite-build mode the parent project already supplies LIBRARIES,
# ROOT_LIBRARIES and the global include dirs.

if(NOT DEFINED LIBRARIES)
    cmake_minimum_required(VERSION 3.6)
    project(genQE_3N_project)
    set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/../../cmake_modules")
    find_package(ROOT REQUIRED)
    set(CMAKE_CXX_STANDARD 17)
    set(SRC_ROOT "${CMAKE_SOURCE_DIR}/../..")
    set(LIBRARIES "${SRC_ROOT}/libraries")
    include_directories(
        ${ROOT_INCLUDE_DIR}
        ${SRC_ROOT}/include
        ${LIBRARIES}
    )
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Absolute location of the 3N density tables, so the binaries do not have to
# be run from this directory (overridable at run time via $GCF_3N_ARRAY_DIR).
add_compile_definitions(GCF_3N_ARRAY_DIR="${CMAKE_CURRENT_SOURCE_DIR}/array")

# ---------------------------------------------------------------------------
# Optional GENIE-FSI integration (matches genQE_FSI's USE_FSI block).
# Defaults ON for standalone builds; in suite builds the user can opt in via
# -DUSE_FSI=ON.

option(USE_FSI "Enable FSI via GENIE intranuclear cascade" OFF)
set(GENIE_SOURCE_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/../genQE_FSI/Generator-R-3_06_02"
    CACHE PATH "Path to built GENIE source tree (contains src/ and lib/)")

if(USE_FSI)
    if(NOT EXISTS "${GENIE_SOURCE_DIR}/src/Framework/GHEP/GHepRecord.h")
        message(WARNING "USE_FSI=ON but GENIE headers not found in ${GENIE_SOURCE_DIR}/src - disabling FSI build.")
        set(USE_FSI OFF)
    endif()
endif()

if(USE_FSI)
    if(NOT EXISTS "${GENIE_SOURCE_DIR}/lib")
        message(WARNING "USE_FSI=ON but GENIE libs directory not found: ${GENIE_SOURCE_DIR}/lib - disabling FSI build.")
        set(USE_FSI OFF)
    endif()
endif()

if(USE_FSI)
    add_compile_definitions(USE_FSI)
    # Reuse the 2N project's GENIE XML overrides so the 3N FSI side picks up
    # the same NuclearModel = LocalFGM choice without sourcing run_SRC.sh.
    add_compile_definitions(GENIE_OVERRIDE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/../genQE_FSI/config/genie_override")
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

    if(NOT LOG4CPP_INCLUDE_DIR OR NOT LOG4CPP_LIBRARY OR NOT LIBXML2_INCLUDE_DIR OR NOT LIBXML2_LIBRARY)
        message(WARNING "USE_FSI=ON but log4cpp/libxml2 not found - disabling FSI build.")
        set(USE_FSI OFF)
    else()
        include_directories(${LOG4CPP_INCLUDE_DIR} ${LIBXML2_INCLUDE_DIR})

        file(GLOB GENIE_SHARED_LIBS "${GENIE_SOURCE_DIR}/lib/libG*.dylib")
        if(NOT GENIE_SHARED_LIBS)
            file(GLOB GENIE_SHARED_LIBS "${GENIE_SOURCE_DIR}/lib/libG*.so")
        endif()
        if(NOT GENIE_SHARED_LIBS)
            message(WARNING "No GENIE shared libs found in ${GENIE_SOURCE_DIR}/lib - disabling FSI build.")
            set(USE_FSI OFF)
        else()
            set(CMAKE_BUILD_RPATH "${GENIE_SOURCE_DIR}/lib")
        endif()
    endif()
endif()

# ---------------------------------------------------------------------------
# Targets
set(GCF_LIBS_3N
    ${LIBRARIES}/cross_sections/eNCrossSection.cc
    ${LIBRARIES}/generator/gcfGenerator.cc
    ${LIBRARIES}/nucleus/gcfNucleus.cc
    ${LIBRARIES}/nucleus/gcfSRC.cc
)

add_executable(genQE_3N genQE_3N.cc ${GCF_LIBS_3N} QEGenerator_3N.cc)
target_link_libraries(genQE_3N ${ROOT_LIBRARIES})

add_executable(theta_plot analysis/theta_plot.cpp ${GCF_LIBS_3N} QEGenerator_3N.cc)
target_link_libraries(theta_plot ${ROOT_LIBRARIES})

add_executable(theta_plot_fast analysis/theta_plot_fast.cpp ${GCF_LIBS_3N} QEGenerator_3N.cc)
target_link_libraries(theta_plot_fast ${ROOT_LIBRARIES})

add_executable(theta_plot_sinelaw analysis/theta_plot_sinelaw.cpp ${GCF_LIBS_3N} QEGenerator_3N.cc)
target_link_libraries(theta_plot_sinelaw ${ROOT_LIBRARIES})

add_executable(SRC_analysis_3N analysis/SRC_analysis_3N.cpp ${GCF_LIBS_3N} QEGenerator_3N.cc)
target_link_libraries(SRC_analysis_3N ${ROOT_LIBRARIES})

add_executable(theta_scan analysis/theta_scan.cpp ${GCF_LIBS_3N} QEGenerator_3N.cc)
target_link_libraries(theta_scan ${ROOT_LIBRARIES})

if(USE_FSI)
    add_executable(genQE_3N_FSI
        genQE_3N_FSI.cc
        QEGenerator_3N.cc
        QEGeneratorFSI_3N.cc
        ${GCF_LIBS_3N}
        ${LIBRARIES}/fsi/GenieFSIHelpers.cc
    )
    target_link_libraries(genQE_3N_FSI
        ${ROOT_LIBRARIES}
        ${GENIE_SHARED_LIBS}
        ${LOG4CPP_LIBRARY}
        ${LIBXML2_LIBRARY}
    )
endif()
