# Geant4 de-excitation of the GCF A-2 remnant (tools/run_eic_sim.sh --deexcite). Build inside eic-shell:
#   cmake -S tools/deexcite -B tools/deexcite/build && cmake --build tools/deexcite/build
cmake_minimum_required(VERSION 3.16)
project(deexcite CXX)
set(CMAKE_CXX_STANDARD 17)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
find_package(Geant4 REQUIRED)
find_package(HepMC3 REQUIRED)
add_executable(deexcite deexcite.cc)
target_include_directories(deexcite PRIVATE ${HEPMC3_INCLUDE_DIR})
target_link_libraries(deexcite PRIVATE ${Geant4_LIBRARIES} ${HEPMC3_LIBRARIES})
