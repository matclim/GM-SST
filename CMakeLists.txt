cmake_minimum_required(VERSION 3.16)
project(StrawTracker VERSION 1.0 LANGUAGES CXX)

# ── C++ standard ──────────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Dependencies ──────────────────────────────────────────────────────────────
find_package(Geant4 REQUIRED ui_all vis_all)
find_package(GeoModelCore REQUIRED COMPONENTS GeoModelKernel)
find_package(GeoModelG4   REQUIRED COMPONENTS GeoModel2G4)
find_package(ROOT         REQUIRED COMPONENTS Core RIO Tree)

# GeoModel installs everything under one prefix.
# GeoModelCore_DIR is the cmake config dir, typically:
#   <prefix>/lib/cmake/GeoModelCore   ->  three levels up = <prefix>
get_filename_component(GEOMODEL_PREFIX "${GeoModelCore_DIR}/../../.." ABSOLUTE)
set(GEOMODEL_INCLUDE_DIR "${GEOMODEL_PREFIX}/include")
set(GEOMODEL_LIB_DIR     "${GEOMODEL_PREFIX}/lib")
message(STATUS "GeoModel prefix  : ${GEOMODEL_PREFIX}")
message(STATUS "GeoModel includes: ${GEOMODEL_INCLUDE_DIR}")

# Eigen3 -- GeoModelKernel/GeoDefinitions.h includes <Eigen/Dense>.
# Try the CMake config package first; fall back to a path search.
find_package(Eigen3 QUIET)
if(TARGET Eigen3::Eigen)
    get_target_property(EIGEN3_INCLUDE_DIR Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
    message(STATUS "Found Eigen3 (imported target): ${EIGEN3_INCLUDE_DIR}")
else()
    find_path(EIGEN3_INCLUDE_DIR
        NAMES Eigen/Dense
        HINTS
            ${GEOMODEL_PREFIX}/include
            ${GEOMODEL_PREFIX}/include/eigen3
            /usr/include/eigen3
            /usr/local/include/eigen3
            /opt/local/include/eigen3
        PATH_SUFFIXES eigen3
    )
    if(NOT EIGEN3_INCLUDE_DIR)
        message(FATAL_ERROR
            "Eigen3 not found.\n"
            "Install it (e.g. 'sudo apt install libeigen3-dev') or pass its location:\n"
            "  cmake .. -DEIGEN3_INCLUDE_DIR=/path/to/include")
    endif()
    message(STATUS "Found Eigen3 (path search): ${EIGEN3_INCLUDE_DIR}")
endif()

include(${Geant4_USE_FILE})
include(${ROOT_USE_FILE})

# ── Sources ───────────────────────────────────────────────────────────────────
set(SOURCES
    src/MaterialManager.cpp
    src/MagneticFieldRegion.cpp
    src/StrawTrackerBuilder.cpp
    src/TrackerDetectorConstruction.cpp
    src/TrackerSD.cpp
    src/TrackerRunAction.cpp
    src/TrackerEventAction.cpp
    src/TrackerActionInitialization.cpp
    src/TrackerPrimaryGeneratorAction.cpp
)

# ── Executable ────────────────────────────────────────────────────────────────
add_executable(run_StrawTracker exec/main.cpp ${SOURCES})

# Enforce C++20 on the target -- takes precedence over anything Geant4's
# USE_FILE injects.
target_compile_features(run_StrawTracker PRIVATE cxx_std_20)

target_include_directories(run_StrawTracker
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${GEOMODEL_INCLUDE_DIR}
        ${EIGEN3_INCLUDE_DIR}
)

target_link_libraries(run_StrawTracker
    PRIVATE
        ${Geant4_LIBRARIES}
        -L${GEOMODEL_LIB_DIR} -lGeoModelKernel -lGeoModel2G4
        ${ROOT_LIBRARIES}
)

# ── Copy vis macro to build dir ───────────────────────────────────────────────
configure_file(straw_vis.mac ${CMAKE_CURRENT_BINARY_DIR}/straw_vis.mac COPYONLY)

# ── Summary ───────────────────────────────────────────────────────────────────
message(STATUS "")
message(STATUS "=== StrawTracker build configuration ===")
message(STATUS "  Geant4 version : ${Geant4_VERSION}")
message(STATUS "  ROOT version   : ${ROOT_VERSION}")
message(STATUS "  GeoModel prefix: ${GEOMODEL_PREFIX}")
message(STATUS "  Eigen3 include : ${EIGEN3_INCLUDE_DIR}")
message(STATUS "  C++ standard   : ${CMAKE_CXX_STANDARD}")
message(STATUS "========================================")
message(STATUS "")
