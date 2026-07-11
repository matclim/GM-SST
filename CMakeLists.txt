cmake_minimum_required(VERSION 3.16)
project(StrawTracker VERSION 1.0 LANGUAGES CXX)

# ── Default build type ────────────────────────────────────────────────────────
# A bare `cmake ..` otherwise leaves CMAKE_BUILD_TYPE empty -> no NDEBUG, which
# leaves ACTS' internal debug asserts (e.g. the field-interpolation bounds
# check) live. Default to Release so runs are optimized and assert-free.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    message(STATUS "No build type specified -- defaulting to Release")
endif()

# ── C++ standard ──────────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ── Dependencies ──────────────────────────────────────────────────────────────
find_package(Geant4 REQUIRED ui_all vis_all)
find_package(GeoModelCore REQUIRED COMPONENTS GeoModelKernel)
find_package(GeoModelG4   REQUIRED COMPONENTS GeoModel2G4)
find_package(ROOT         REQUIRED COMPONENTS Core RIO Tree)
find_package(GeoModelIO REQUIRED COMPONENTS GeoModelWrite)

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

# ── Geant4 simulation sources ─────────────────────────────────────────────────
set(SOURCES
    src/MaterialManager.cpp
    src/MagneticFieldRegion.cpp
    src/ShipFieldMap.cpp
    src/StrawTrackerBuilder.cpp
    src/TrackerDetectorConstruction.cpp
    src/TrackerSD.cpp
    src/TrackerRunAction.cpp
    src/TrackerEventAction.cpp
    src/TrackerActionInitialization.cpp
    src/TrackerPrimaryGeneratorAction.cpp
)

# ── Simulation executable ─────────────────────────────────────────────────────
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
        -L${GEOMODEL_LIB_DIR} -lGeoModelKernel -lGeoModel2G4 -lGeoModelWrite -lGeoModelDBManager
        ${ROOT_LIBRARIES}
)

# ── Copy vis macro to build dir ───────────────────────────────────────────────
configure_file(straw_vis.mac ${CMAKE_CURRENT_BINARY_DIR}/straw_vis.mac COPYONLY)

# ══════════════════════════════════════════════════════════════════════════════
#  ACTS-based reconstruction (tracking + vertexing)
#
#  Standalone library `ship_reco` + app `run_reco`, built ONLY when ACTS is
#  found so the simulation above always configures. To point CMake at ACTS:
#     - source an LCG view on CVMFS that contains it, OR
#     - configure with -DCMAKE_PREFIX_PATH=<acts-install>  (or -DActs_DIR=...)
# ══════════════════════════════════════════════════════════════════════════════
option(STRAWTRACKER_BUILD_RECO "Build the ACTS-based reconstruction" ON)

if(STRAWTRACKER_BUILD_RECO)
    # ACTS ships a CMake config package. Core is enough for geometry, field,
    # propagation, CKF/KF fitting, and vertexing (all in Acts::Core).
    find_package(Acts CONFIG QUIET COMPONENTS Core)

    if(Acts_FOUND)
        # Namespaced target in modern ACTS; fall back to the old name.
        if(TARGET Acts::Core)
            set(ACTS_CORE_TARGET Acts::Core)
        elseif(TARGET ActsCore)
            set(ACTS_CORE_TARGET ActsCore)
        else()
            message(FATAL_ERROR "ACTS found but neither Acts::Core nor ActsCore target exists")
        endif()

        # Reco modules. Add the others as you build them:
        #   geometry/ShipTrackingGeometry.cpp
        #   edm/ShipMeasurement.cpp
        #   seeding/ShipSeeder.cpp
        #   finding/ShipTrackFinding.cpp
        #   fitting/ShipTrackFitting.cpp
        #   vertexing/ShipVertexing.cpp
        set(RECO_SOURCES
            src/ActsShipField.cpp
            src/ShipGeometry.cpp
            src/ShipStrawGeometry.cpp
            src/ShipHitReader.cpp
            src/ShipSeeder.cpp
            src/ShipTrackFitter.cpp
            src/ShipVertexer.cpp
        )

        add_library(ship_reco STATIC ${RECO_SOURCES})
        target_compile_features(ship_reco PUBLIC cxx_std_20)
        target_include_directories(ship_reco
            PUBLIC
                ${CMAKE_CURRENT_SOURCE_DIR}/include
        )
        target_link_libraries(ship_reco
            PUBLIC
                ${ACTS_CORE_TARGET}
                ROOT::Core ROOT::RIO ROOT::Tree
        )

        # Reconstruction driver.
        add_executable(run_reco src/run_reco.cpp)
        target_link_libraries(run_reco PRIVATE ship_reco)
        target_compile_features(run_reco PRIVATE cxx_std_20)

        message(STATUS "ACTS found (${Acts_VERSION}) -> building ship_reco [${ACTS_CORE_TARGET}]")
    else()
        message(STATUS "ACTS not found -> skipping reconstruction. "
                       "Set CMAKE_PREFIX_PATH or Acts_DIR, or source an LCG view.")
    endif()
endif()

# ── Summary ───────────────────────────────────────────────────────────────────
message(STATUS "")
message(STATUS "=== StrawTracker build configuration ===")
message(STATUS "  Geant4 version : ${Geant4_VERSION}")
message(STATUS "  ROOT version   : ${ROOT_VERSION}")
message(STATUS "  GeoModel prefix: ${GEOMODEL_PREFIX}")
message(STATUS "  Eigen3 include : ${EIGEN3_INCLUDE_DIR}")
message(STATUS "  C++ standard   : ${CMAKE_CXX_STANDARD}")
message(STATUS "  Build reco     : ${STRAWTRACKER_BUILD_RECO} (ACTS found: ${Acts_FOUND})")
message(STATUS "========================================")
message(STATUS "")
