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
//              ├── StrawFrame_j   (hollow rectangle, material = Aluminum by
//              │                   default; GeoShapeSubtraction = outer − hole)
//              ├── SubLayer_0     (air slab with N straws, at z = -kStrawRadius)
//              └── SubLayer_1     (air slab, staggered, at z = +kStrawRadius)
//
// Stereo angles (rotation of the layer volume around Z):
//   j=0:  +2.3°   j=1: -2.3°   j=2: +2.3°   j=3: -2.3°
//
// The second sub-layer (k=1) is staggered:
//   - shifted by one straw radius in Z  (+kStrawRadius)
//   - shifted by half a pitch in Y      (+kStrawRadius)
//
// Frame (FairShip-style, one per view):
//   - outer half-sizes: aperture + kFrameWidth* material bar
//   - inner aperture:   slightly larger than the straw pattern so the sub-layer
//                       envelope fits cleanly inside with no geometry overlap
//   - built as (outer_box − inner_box) via GeoShapeSubtraction
//   - placed at z = 0 inside the layer envelope, so it rotates with the view
// ─────────────────────────────────────────────────────────────────────────────

#include "StrawTrackerBuilder.h"
#include "MaterialManager.h"

// GeoModel kernel headers
#include "GeoModelKernel/GeoBox.h"
#include "GeoModelKernel/GeoTube.h"
#include "GeoModelKernel/GeoShapeSubtraction.h"
#include "GeoModelKernel/GeoLogVol.h"
#include "GeoModelKernel/GeoPhysVol.h"
#include "GeoModelKernel/GeoFullPhysVol.h"
#include "GeoModelKernel/GeoTransform.h"
#include "GeoModelKernel/GeoNameTag.h"
#include "GeoModelKernel/GeoIdentifierTag.h"
#include "GeoModelKernel/Units.h"
#include "GeoModelKernel/GeoDefinitions.h"
#include "GeoModelKernel/GeoXF.h"
#include "GeoModelDBManager/GMDBManager.h"
#include "GeoModelWrite/WriteGeoModel.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <iostream>

namespace GU = GeoModelKernelUnits;

// ── Geometry constants that are reused across builders ───────────────────────
namespace {
    // Aperture clearance beyond the nominal active area. Needs to be a bit
    // larger than the nominal straw pattern so the staggered sub-layer (which
    // pushes an extra kStrawRadius in Y) and the straw outer radius fit inside.
    constexpr double kApClearX = 5.0;   // mm
    constexpr double kApClearY = 15.0;  // mm

    // Small clearance between layer envelope and frame outer, and between the
    // sub-layer envelope and the frame aperture (to keep Geant4 happy).
    constexpr double kEnvClearance = 5.0; // mm

    // ── Derived half-sizes (pure mm, multiply by GU::mm at use site) ─────────
    // Frame aperture (inner hole).
    inline constexpr double apHalfX_mm() {
        return StrawTrackerBuilder::kStationX / 2.0 + kApClearX;   // 2005
    }
    inline constexpr double apHalfY_mm() {
        return StrawTrackerBuilder::kStationY / 2.0 + kApClearY;   // 3015
    }
    // Frame outer.
    inline constexpr double frHalfX_mm() {
        return apHalfX_mm() + StrawTrackerBuilder::kFrameWidthX;   // 2105
    }
    inline constexpr double frHalfY_mm() {
        return apHalfY_mm() + StrawTrackerBuilder::kFrameWidthY;   // 3115
    }
    // Layer envelope.
    inline constexpr double layHalfX_mm() {
        return frHalfX_mm() + kEnvClearance;                       // 2110
    }
    inline constexpr double layHalfY_mm() {
        return frHalfY_mm() + kEnvClearance;                       // 3120
    }
    inline constexpr double layHalfZ_mm() {
        return StrawTrackerBuilder::kFrameHalfZ + kEnvClearance;   // 27
    }
}

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

    // World box sized to contain all 4 stations with generous XY margin to
    // accommodate the frames and the magnetic field volume (placed between
    // stations 1 and 2).
    const double worldX = 2500.0  * GU::mm;   // half-lengths for GeoBox
    const double worldY = 3500.0  * GU::mm;
    const double worldZ = 7500.0  * GU::mm;

    auto* worldBox  = new GeoBox(worldX, worldY, worldZ);
    auto* worldLog  = new GeoLogVol("World", worldBox, Mm.Air());
    auto* worldPhys = new GeoPhysVol(worldLog);

    // World centre in lab frame Z (see README).
    constexpr double worldZOrigin = 31000.0 * GU::mm;

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

void StrawTrackerBuilder::writeDB(GeoPhysVol* world, const std::string& filename) {
    // Remove any stale DB from a previous run so GMDBManager starts with a
    // clean schema. Without this, SQLite prints "table already exists" and
    // "UNIQUE constraint failed" errors on every run. `std::remove` silently
    // succeeds if the file is missing, so this is always safe.
    std::remove(filename.c_str());

    GMDBManager db(filename);
    GeoModelIO::WriteGeoModel writer(db);
    world->exec(&writer);
    writer.saveToDB();
    std::cout << "[StrawTrackerBuilder] Geometry written to " << filename << "\n";
}

// =============================================================================
// buildStation
// =============================================================================
GeoPhysVol* StrawTrackerBuilder::buildStation(int stationID) {
    auto& Mm = MaterialManager::instance();

    // Layer envelope after stereo rotation expands in XY. Compute the
    // bounding box of the rotated rectangle with the worst-case angle.
    const double angleRad = kStereoAngle * M_PI / 180.0;
    const double cosA = std::cos(angleRad);
    const double sinA = std::sin(angleRad);
    const double rotHalfX = layHalfX_mm() * cosA + layHalfY_mm() * sinA;
    const double rotHalfY = layHalfY_mm() * cosA + layHalfX_mm() * sinA;

    const double stHalfX = (rotHalfX + 30.0) * GU::mm;
    const double stHalfY = (rotHalfY + 30.0) * GU::mm;

    // Z stack: 4 layers spaced by layerPitch, centred on the station centre.
    const double layerGap   = 5.0;
    const double layerPitch = 2.0 * layHalfZ_mm() + layerGap;
    const double stackHalfZ = 0.5 * (kNLayers - 1) * layerPitch + layHalfZ_mm();
    const double stHalfZ    = (stackHalfZ + 10.0) * GU::mm;

    auto* stationBox  = new GeoBox(stHalfX, stHalfY, stHalfZ);
    auto* stationLog  = new GeoLogVol("Station", stationBox, Mm.Air());
    auto* stationPhys = new GeoPhysVol(stationLog);

    for (int iLayer = 0; iLayer < kNLayers; ++iLayer) {
        const double signedAngle = stereoSign(iLayer) * kStereoAngle;
        GeoPhysVol* layerPhys = buildLayer(iLayer, signedAngle);

        // Z position of this layer inside the station.
        const double zLay = -0.5 * (kNLayers - 1) * layerPitch + iLayer * layerPitch;
        const double signedAngleRad = signedAngle * M_PI / 180.0;
        GeoTrf::Transform3D xfLayer =
            GeoTrf::TranslateZ3D(zLay * GU::mm) * GeoTrf::RotateZ3D(signedAngleRad);

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
// buildLayer  (= one "view" in FairShip terminology)
// =============================================================================
GeoPhysVol* StrawTrackerBuilder::buildLayer(int layerID, double /*stereoAngleDeg*/) {
    // Layer envelope holds a material frame plus 2 sub-layers of straws.
    // The stereo rotation is applied by the parent station when this volume
    // is placed.

    auto& Mm = MaterialManager::instance();

    const double layHalfX = layHalfX_mm() * GU::mm;
    const double layHalfY = layHalfY_mm() * GU::mm;
    const double layHalfZ = layHalfZ_mm() * GU::mm;

    auto* layerBox  = new GeoBox(layHalfX, layHalfY, layHalfZ);
    auto* layerLog  = new GeoLogVol("StrawLayer", layerBox, Mm.Air());
    auto* layerPhys = new GeoPhysVol(layerLog);

    // ── Material frame (FairShip-style, one per view) ────────────────────
    {
        GeoPhysVol* framePhys = buildFrame(layerID);
        auto* nameTag = new GeoNameTag("StrawFrame_" + std::to_string(layerID));
        auto* idTag   = new GeoIdentifierTag(100 + layerID);
        auto* xf      = new GeoTransform(GeoTrf::Transform3D::Identity());
        layerPhys->add(nameTag);
        layerPhys->add(idTag);
        layerPhys->add(xf);
        layerPhys->add(framePhys);
    }

    // ── Two sub-layers of straws ─────────────────────────────────────────
    // Centres are at z = ±(kStrawRadius + 0.05) mm in the layer frame. The
    // extra 0.05 mm separation — on top of the 0.5 mm z-margin inside the
    // sub-layer envelope — is what keeps the two sub-layer envelopes from
    // overlapping each other at z = 0 (they used to overlap by 1 mm in air,
    // which CheckOverlaps flags once pSurfChk is enabled).
    const double dz0 = -(kStrawRadius + 0.55) * GU::mm;  // -10.55 mm
    const double dz1 = +(kStrawRadius + 0.55) * GU::mm;  // +10.55 mm

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

    {
        GeoPhysVol* sl1 = buildSubLayer(true);
        auto* nameTag = new GeoNameTag("SubLayer_1");
        auto* idTag   = new GeoIdentifierTag(1);
        // Only a Z offset here — the half-pitch Y stagger is applied INSIDE
        // buildSubLayer via the `shifted` flag, so both sub-layer envelopes
        // share the same symmetric XY footprint and fit cleanly in the frame.
        auto* xf = new GeoTransform(GeoTrf::TranslateZ3D(dz1));
        layerPhys->add(nameTag);
        layerPhys->add(idTag);
        layerPhys->add(xf);
        layerPhys->add(sl1);
    }

    return layerPhys;
}

// =============================================================================
// buildFrame  (FairShip-style material frame around the view)
// =============================================================================
GeoPhysVol* StrawTrackerBuilder::buildFrame(int layerID) {
    static int frameUID = 0;
    auto& Mm = MaterialManager::instance();

    // Inner hole dimensions: slightly larger than the active straw pattern so
    // the sub-layer envelopes fit inside with clearance and no overlap.
    const double apHalfX = apHalfX_mm() * GU::mm;
    const double apHalfY = apHalfY_mm() * GU::mm;

    // Outer dimensions: aperture + frame bar width.
    const double frHalfX = frHalfX_mm() * GU::mm;
    const double frHalfY = frHalfY_mm() * GU::mm;
    const double frHalfZ = kFrameHalfZ  * GU::mm;

    // Build outer box and (slightly thicker) inner box to punch through cleanly.
    auto* outerBox = new GeoBox(frHalfX, frHalfY, frHalfZ);
    auto* innerBox = new GeoBox(apHalfX, apHalfY, frHalfZ + 1.0 * GU::mm);

    // Subtract: outer − inner = hollow rectangular frame.
    auto* frameShape = new GeoShapeSubtraction(outerBox, innerBox);

    auto* frameMat = Mm.frameMaterialByName(m_frameMaterialName);
    // Globally-unique log-volume name so visualisation can target each frame.
    const std::string logName = "StrawFrameLV_" + std::to_string(frameUID++)
                              + "_layer" + std::to_string(layerID);
    auto* frameLog  = new GeoLogVol(logName, frameShape, frameMat);
    auto* framePhys = new GeoPhysVol(frameLog);

    return framePhys;
}

// =============================================================================
// buildSubLayer
// =============================================================================
GeoPhysVol* StrawTrackerBuilder::buildSubLayer(bool shifted) {
    // A sub-layer is a thin air slab containing kNStraws straws.
    // Both variants (nominal + shifted) have the same symmetric XY footprint
    // so they fit cleanly inside the frame aperture. The half-pitch Y stagger
    // between the two is applied *here* (via the `shifted` flag) rather than
    // as a translation in buildLayer, which would push one envelope into the
    // frame material on the +Y side and cause a sibling-volume overlap.

    static int strawUID = 0;
    auto& Mm = MaterialManager::instance();

    const double pitch = 2.0 * kStrawRadius;  // 20 mm
    // Half-pitch stagger applied per-straw in the shifted sub-layer.
    const double yStagger = shifted ? kStrawRadius : 0.0;   // +10 mm if shifted

    // Sub-layer envelope: symmetric, large enough to contain the staggered
    // straw pattern (which reaches y = ±3000 + kStrawRadius after the
    // stagger), and small enough to fit inside the frame aperture with a
    // kFrameClearance-mm gap on each side.
    const double slHalfX = (apHalfX_mm() - kFrameClearance) * GU::mm;
    const double slHalfY = (apHalfY_mm() - kFrameClearance) * GU::mm;
    const double slHalfZ = (kStrawRadius + 0.5)             * GU::mm;

    const std::string name = shifted ? "SubLayer_shifted" : "SubLayer_nominal";
    auto* slBox  = new GeoBox(slHalfX, slHalfY, slHalfZ);
    auto* slLog  = new GeoLogVol(name, slBox, Mm.Air());
    auto* slPhys = new GeoPhysVol(slLog);

    const double yStart = -(kNStraws - 1) * 0.5 * pitch;

    for (int iStraw = 0; iStraw < kNStraws; ++iStraw) {
        // Nominal sub-layer: straw centers at y = -2990, -2970, ..., +2990
        // Shifted sub-layer: straw centers at y = -2980, -2960, ..., +3000
        const double yStraw = (yStart + iStraw * pitch + yStagger) * GU::mm;

        auto* nameTag = new GeoNameTag("Straw_" + std::to_string(iStraw));
        auto* idTag   = new GeoIdentifierTag(iStraw);
        auto* xf = new GeoTransform(
            GeoTrf::TranslateY3D(yStraw) * GeoTrf::RotateY3D(M_PI / 2.0)
        );

        slPhys->add(nameTag);
        slPhys->add(idTag);
        slPhys->add(xf);
        slPhys->add(buildStraw(strawUID++));
    }

    return slPhys;
}

// =============================================================================
// buildStraw
// =============================================================================
GeoPhysVol* StrawTrackerBuilder::buildStraw(int uid) {
    // A straw is built as a SOLID cylinder (the "wall" LV) filled with Mylar,
    // with a gas daughter that occupies the interior:
    //
    //   Outer tube (wall):  rMin = 0, rMax = kStrawRadius,  Mylar
    //   Inner tube (gas):   rMin = 0, rMax = rGas,          ArCO2
    //
    // Geant4 uses the gas's material inside the gas, and the wall's material
    // everywhere else in the wall, so the physics is identical to the "wall
    // is a hollow shell" model — but now the gas is properly contained
    // inside the wall's shape (no mother-daughter overlap, which Geant4's
    // overlap checker catches on every single straw if the wall is hollow).
    //
    // GeoTube axis is along local Z; the parent sub-layer rotates this to X.

    auto& Mm = MaterialManager::instance();

    const double rGas  = (kStrawRadius - kWallThick) * GU::mm;
    const double rWall = kStrawRadius                 * GU::mm;
    const double half  = (kStrawLength / 2.0)         * GU::mm;

    auto* wallTube = new GeoTube(0.0, rWall, half);
    auto* wallLog  = new GeoLogVol("StrawWall_" + std::to_string(uid), wallTube, Mm.Mylar());
    auto* wallPhys = new GeoPhysVol(wallLog);

    auto* gasTube = new GeoTube(0.0, rGas, half);
    auto* gasLog  = new GeoLogVol("StrawGas_" + std::to_string(uid), gasTube, Mm.ArCO2());
    auto* gasPhys = new GeoPhysVol(gasLog);

    wallPhys->add(new GeoNameTag("StrawGas"));
    wallPhys->add(gasPhys);

    return wallPhys;
}
