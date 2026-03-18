// StrawTrackerBuilder.cpp
// ─────────────────────────────────────────────────────────────────────────────
// Full GeoModel geometry description of the straw tracker.
//
// Coordinate convention (GeoModel / Geant4 units are mm):
//   X  — along the straw axis (horizontal, 4 m)
//   Y  — vertical, along the straw pitch direction (unrotated)
//   Z  — beam direction
//
// Hierarchy
// ─────────
//   World  (air box)
//    └── Station_i  (air box, i = 0..3)
//         └── Layer_j  (air box, rotated by stereo angle, j = 0..3)
//              └── SubLayer_k  (air box, k = 0..1; k=1 shifted ½ pitch in Y)
//                   └── Straw_n  (wall tube + gas tube, n = 0..N-1)
//
// Stereo angles (rotation of the layer volume around Z):
//   j=0:  +2.3°   j=1: -2.3°   j=2: +2.3°   j=3: -2.3°
//
// The second sub-layer (k=1) is staggered:
//   - shifted by one straw radius in Z  (+kStrawRadius)
//   - shifted by half a pitch in Y      (+kStrawRadius)
// ─────────────────────────────────────────────────────────────────────────────

#include "StrawTrackerBuilder.h"
#include "MaterialManager.h"

// GeoModel kernel headers
#include "GeoModelKernel/GeoBox.h"
#include "GeoModelKernel/GeoTube.h"
#include "GeoModelKernel/GeoLogVol.h"
#include "GeoModelKernel/GeoPhysVol.h"
#include "GeoModelKernel/GeoFullPhysVol.h"
#include "GeoModelKernel/GeoTransform.h"
#include "GeoModelKernel/GeoNameTag.h"
#include "GeoModelKernel/GeoIdentifierTag.h"
#include "GeoModelKernel/Units.h"
#include "GeoModelKernel/GeoDefinitions.h"
#include "GeoModelKernel/GeoXF.h"  

//G4 visualisations
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

// Standard transforms
#include <cmath>
#include <string>
#include <iostream>
// ── Unit aliases ──────────────────────────────────────────────────────────────
namespace GU = GeoModelKernelUnits;

// ── Helper: sign of stereo angle for a given layer index ─────────────────────
static double stereoSign(int layerID) {
    // Alternating: layer 0,2 → +1;  layer 1,3 → -1
    return (layerID % 2 == 0) ? +1.0 : -1.0;
}

// =============================================================================
// buildWorld
// =============================================================================
GeoPhysVol* StrawTrackerBuilder::buildWorld() {
    auto& Mm = MaterialManager::instance();

    // World box: large enough to contain all stations.
    // X: straw length + margin  → 4500 mm
    // Y: station height + margin → 6500 mm
    // Z: from first to last station + margin:
    //    last station centre = 35500 mm, first = 26500 mm → span 9000 mm, add 2000
    const double worldX = 2500.0  * GU::mm;   // half-lengths for GeoBox
    const double worldY = 3500.0  * GU::mm;
    const double worldZ = 7500.0  * GU::mm;   // half of 15000 mm span (centred at ~31000)

    auto* worldBox  = new GeoBox(worldX, worldY, worldZ);
    auto* worldLog  = new GeoLogVol("World", worldBox, Mm.Air());
    auto* worldPhys = new GeoPhysVol(worldLog);

    // ── Place 4 stations ──────────────────────────────────────────────────────
    // World is centred at origin in X,Y but we shift Z to place the world box
    // symmetrically around the midpoint of the station span.
    // Station z-centres: 26500, 29000, 34000, 35500 mm
    // Midpoint: (26500+35500)/2 = 31000 mm  → world origin at z=31000 in the
    // lab frame. We place everything relative to world origin = (0,0,31000)mm
    // but GeoModel works in its own local frame, so station positions below
    // are relative to the world centre.

    constexpr double worldZOrigin = 31000.0 * GU::mm;  // world centre in lab Z

    for (int iStation = 0; iStation < kNStations; ++iStation) {
        GeoPhysVol* stationPhys = buildStation(iStation);

        const double zPos = kStationZ[iStation] * GU::mm - worldZOrigin;

        auto* nameTag = new GeoNameTag("Station_" + std::to_string(iStation));
        auto* idTag   = new GeoIdentifierTag(iStation);
        auto* xf      = new GeoTransform(GeoTrf::TranslateZ3D(zPos));

        worldPhys->add(nameTag);
        worldPhys->add(idTag);
        worldPhys->add(xf);
        worldPhys->add(stationPhys);
    }

    return worldPhys;
}

// =============================================================================
// buildStation
// =============================================================================
GeoPhysVol* StrawTrackerBuilder::buildStation(int stationID) {
    auto& Mm = MaterialManager::instance();

    // Station envelope: just large enough to hold 4 straw layers.
    // Each layer has thickness = 2 sub-layers × 2×kStrawRadius = 4×kStrawRadius.
    // Layer envelope also needs a thin margin in Z for the rotated placement.
    // With stereo angle α = 2.3° the Y-extent of a rotated layer is:
    //   Y_rot = kStationY * cos(α) + kStrawLength * sin(α)
    //         ≈ 6000 * 0.9992 + 4000 * 0.040 ≈ 6155 mm  (use 6500 mm to be safe)
    //
    // Stacking 4 layers in Z:
    //   each layer envelope half-thickness in Z ≈ 2 * kStrawRadius = 20 mm
    //   plus gaps → use 5 mm gap between layers
    //   total Z half-extent = 4*(20 + 5) = 100 mm  → use 120 mm (generous)

    const double stHalfX = (kStrawLength / 2.0 + 50.0)  * GU::mm;
    const double stHalfY = (kStationY    / 2.0 + 300.0) * GU::mm;
    const double stHalfZ = 120.0 * GU::mm;

    auto* stationBox  = new GeoBox(stHalfX, stHalfY, stHalfZ);
    auto* stationLog  = new GeoLogVol("Station", stationBox, Mm.Air());
    auto* stationPhys = new GeoPhysVol(stationLog);

    // Layer z-positions inside station: stack them with 5 mm gap.
    // Layer half-thickness in Z = 2 * kStrawRadius + 1 mm clearance = 21 mm
    const double layerHalfZ = (2.0 * kStrawRadius + 1.0) * GU::mm;  // 21 mm
    const double layerGap   = 5.0 * GU::mm;
    const double layerPitch = 2.0 * layerHalfZ + layerGap;  // 47 mm

    // Centre the stack of 4 layers at z=0 in the station frame.
    const double stackHalfZ = 0.5 * (kNLayers - 1) * layerPitch;

    for (int iLayer = 0; iLayer < kNLayers; ++iLayer) {
        const double signedAngle = stereoSign(iLayer) * kStereoAngle;
        GeoPhysVol* layerPhys = buildLayer(iLayer, signedAngle);

        // Position in Z
        const double zLay = -stackHalfZ + iLayer * layerPitch;

        // Rotation of the layer around Z-axis by the stereo angle
        const double angleRad = signedAngle * M_PI / 180.0;
        // Note: zLay is already in mm and we applied GU::mm above; correct:
        const double zLayMM = -stackHalfZ + iLayer * layerPitch; // pure mm value
        GeoTrf::Transform3D xfLayer = GeoTrf::TranslateZ3D(zLayMM) * GeoTrf::RotateZ3D(angleRad);

        auto* nameTag = new GeoNameTag("Layer_" + std::to_string(iLayer));
        auto* idTag   = new GeoIdentifierTag(iLayer);
        auto* xf      = new GeoTransform(xfLayer);

        stationPhys->add(nameTag);
        stationPhys->add(idTag);
        stationPhys->add(xf);
        stationPhys->add(layerPhys);
    }

    return stationPhys;
}

// =============================================================================
// buildLayer
// =============================================================================
GeoPhysVol* StrawTrackerBuilder::buildLayer(int layerID, double /*stereoAngleDeg*/) {
    // Layer envelope holds 2 sub-layers.
    // Sub-layer 0: centred at dz = -kStrawRadius
    // Sub-layer 1: centred at dz = +kStrawRadius, shifted +½ pitch in Y.
    //
    // The layer envelope is unrotated here; the stereo rotation is applied
    // by the parent station when placing this volume.

    auto& Mm = MaterialManager::instance();

    const double layHalfX = (kStrawLength / 2.0 + 10.0) * GU::mm;
    const double layHalfY = (kStationY    / 2.0 + 20.0) * GU::mm;
    const double layHalfZ = (2.0 * kStrawRadius + 1.0)  * GU::mm;

    auto* layerBox  = new GeoBox(layHalfX, layHalfY, layHalfZ);
    auto* layerLog  = new GeoLogVol("StrawLayer", layerBox, Mm.Air());
    auto* layerPhys = new GeoPhysVol(layerLog);

    const double dz0 = -kStrawRadius * GU::mm;  // sub-layer 0 z-offset
    const double dz1 = +kStrawRadius * GU::mm;  // sub-layer 1 z-offset

    // Sub-layer 0 (not shifted in Y)
    {
        GeoPhysVol* sl0 = buildSubLayer(false);
        auto* nameTag = new GeoNameTag("SubLayer_0");
        auto* idTag   = new GeoIdentifierTag(0);
        auto* xf      = new GeoTransform(GeoTrf::TranslateZ3D(dz0));
        layerPhys->add(nameTag);
        layerPhys->add(idTag);
        layerPhys->add(xf);
        layerPhys->add(sl0);
    }

    // Sub-layer 1 (shifted half-pitch in Y for staggering)
    {
        GeoPhysVol* sl1 = buildSubLayer(true);
        auto* nameTag = new GeoNameTag("SubLayer_1");
        auto* idTag   = new GeoIdentifierTag(1);
        // Translate: +kStrawRadius in Y (half pitch), +kStrawRadius in Z
        auto* xf = new GeoTransform(
            GeoTrf::TranslateZ3D(dz1) * GeoTrf::TranslateY3D(kStrawRadius * GU::mm)
        );
        layerPhys->add(nameTag);
        layerPhys->add(idTag);
        layerPhys->add(xf);
        layerPhys->add(sl1);
    }

    return layerPhys;
}

// =============================================================================
// buildSubLayer
// =============================================================================
GeoPhysVol* StrawTrackerBuilder::buildSubLayer(bool shifted) {
    // A sub-layer is a thin air slab containing kNStraws straws.
    // Straws are placed along Y, centred at y = (n + 0.5) * pitch - kStationY/2
    // where pitch = 2 * kStrawRadius.
    // The "shifted" flag moves the slab by half a pitch (handled by the parent
    // layer transform); the slab itself is always centred at y=0.
    std::cout << "buildSubLayer called!" << std::endl;
    auto& Mm = MaterialManager::instance();

    const double pitch    = 2.0 * kStrawRadius;           // 20 mm
    const double slHalfX  = (kStrawLength / 2.0 + 5.0) * GU::mm;
    const double slHalfY  = (kStationY   / 2.0 + pitch) * GU::mm; // slight extra
    const double slHalfZ  = (kStrawRadius + 0.5)        * GU::mm;

    const std::string name = shifted ? "SubLayer_shifted" : "SubLayer_nominal";
    auto* slBox  = new GeoBox(slHalfX, slHalfY, slHalfZ);
    auto* slLog  = new GeoLogVol(name, slBox, Mm.Air());
    auto* slPhys = new GeoPhysVol(slLog);

    const double yStart = -(kNStraws - 1) * 0.5 * pitch;

    for (int iStraw = 0; iStraw < kNStraws; ++iStraw) {
        const double yStraw = (yStart + iStraw * pitch) * GU::mm;

        auto* nameTag = new GeoNameTag("Straw_" + std::to_string(iStraw));
        auto* idTag   = new GeoIdentifierTag(iStraw);
        auto* xf = new GeoTransform(
            GeoTrf::TranslateY3D(yStraw) * GeoTrf::RotateY3D(M_PI / 2.0)
        );

        slPhys->add(nameTag);
        slPhys->add(idTag);
        slPhys->add(xf);
        slPhys->add(buildStraw());   // fresh GeoPhysVol per straw
    } 

    return slPhys;
}

// =============================================================================
// buildStraw
// =============================================================================
GeoPhysVol* StrawTrackerBuilder::buildStraw() {
    // A straw is a coaxial pair of tubes:
    //   Outer tube (wall):  rMin = kStrawRadius - kWallThick,
    //                       rMax = kStrawRadius,
    //                       half-length = kStrawLength/2
    //                       material: Mylar
    //   Inner tube (gas):   rMin = 0,
    //                       rMax = kStrawRadius - kWallThick,
    //                       half-length = kStrawLength/2
    //                       material: ArCO2
    //
    // GeoTube axis is along local Z; the parent sub-layer rotates this to X.

    auto& Mm = MaterialManager::instance();

    const double rGas  = (kStrawRadius - kWallThick) * GU::mm;
    const double rWall = kStrawRadius                 * GU::mm;
    const double half  = (kStrawLength / 2.0)         * GU::mm;

    // ── Outer (wall) tube ─────────────────────────────────────────────────────
    auto* wallTube = new GeoTube(rGas, rWall, half);
    auto* wallLog  = new GeoLogVol("StrawWall", wallTube, Mm.Mylar());
    auto* wallPhys = new GeoPhysVol(wallLog);

    // ── Inner (gas) tube ──────────────────────────────────────────────────────
    auto* gasTube  = new GeoTube(0.0, rGas, half);
    auto* gasLog   = new GeoLogVol("StrawGas", gasTube, Mm.ArCO2());
    auto* gasPhys  = new GeoPhysVol(gasLog);

    // Place gas inside wall (both share the same axis, no transform needed)
    wallPhys->add(new GeoNameTag("StrawGas"));
    wallPhys->add(gasPhys);

    return wallPhys;
}
