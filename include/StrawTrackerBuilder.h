#pragma once
// StrawTrackerBuilder.h
// Builds the full straw tracker geometry using GeoModel primitives.
//
// Geometry summary
// ─────────────────
//   4 stations, each containing 4 straw layers.
//   Each layer is a double sub-layer of staggered straws.
//   Straw: diameter 2 cm (radius 1 cm), length 4 m, horizontal (along X).
//   Wall:  30 µm Mylar.
//   Fill:  Ar/CO2 70/30.
//   Station active area: 4 m (X) × 6 m (Y).
//
//   Stereo angles per layer (relative to Y axis, rotation around Z):
//     layer 0: +2.3°   (u)
//     layer 1: -2.3°   (v)
//     layer 2: +2.3°   (u)
//     layer 3: -2.3°   (v)
//
//   Station z-centres (mm):
//     station 0: z = 26500 mm
//     station 1: z = 29000 mm
//     station 2: z = 34000 mm
//     station 3: z = 35500 mm

#include "GeoModelKernel/GeoPhysVol.h"
#include "GeoModelKernel/GeoFullPhysVol.h"
#include "GeoModelKernel/GeoDefinitions.h"
#include <array>

class StrawTrackerBuilder {
public:
    // Returns the world physical volume; call once.
    GeoPhysVol* buildWorld();

    // ── Geometry constants ────────────────────────────────────────────────
    static constexpr int    kNStations    = 4;
    static constexpr int    kNLayers      = 4;   // per station
    static constexpr int    kNSubLayers   = 2;   // per layer (staggered pair)

    // Straw dimensions (mm)
    static constexpr double kStrawRadius  = 10.0;   // 1 cm radius → 2 cm diam
    static constexpr double kStrawLength  = 4000.0; // 4 m along X
    static constexpr double kWallThick    = 0.030;  // 30 µm Mylar

    // Station active face (mm)
    static constexpr double kStationX     = 4000.0; // width  (straw length)
    static constexpr double kStationY     = 6000.0; // height

    // Straws per sub-layer
    static constexpr int kNStraws = static_cast<int>(kStationY / (2.0 * kStrawRadius)); // 300

    // Stereo angle magnitude (degrees)
    static constexpr double kStereoAngle  = 2.3;

    // Z positions of station centres (mm) — converted from cm
    static constexpr std::array<double,4> kStationZ = {26500., 29000., 34000., 35500.};

private:
    // ── Internal builders ─────────────────────────────────────────────────

    /// Build one station volume and place its 4 layers inside it.
    GeoPhysVol* buildStation(int stationID);

    /// Build one stereo layer (two sub-layers + all straws).
    /// @param stereoAngleDeg  signed rotation about the straw axis direction
    GeoPhysVol* buildLayer(int layerID, double stereoAngleDeg);

    /// Build one sub-layer of N parallel straws.
    /// @param offsetY  half-pitch offset applied to the second sub-layer
    GeoPhysVol* buildSubLayer(bool shifted);

    /// Build a single straw (wall tube + gas tube).
    GeoPhysVol* buildStraw();
};
