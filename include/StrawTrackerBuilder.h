#pragma once
// StrawTrackerBuilder.h
// Builds the full straw tracker geometry using GeoModel primitives.
//
// Geometry summary
// ─────────────────
//   4 stations, each containing 4 straw layers (= 4 views).
//   Each layer is a double sub-layer of staggered straws, surrounded by a
//   material frame (FairShip-style "frame around the view").
//   Straw: diameter 2 cm (radius 1 cm), length 4 m, horizontal (along X).
//   Wall:  30 µm Mylar.
//   Fill:  Ar/CO2 70/30.
//   Station active (aperture) area: 4 m (X) × 6 m (Y).
//
//   Stereo angles per layer (rotation of the layer volume around Z):
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
//
//   Frame (per layer, FairShip-style):
//     Shape: outer rectangular box minus inner aperture rectangle
//            (a GeoShapeSubtraction, yielding a hollow rectangle).
//     Outer half-sizes (X,Y): aperture + kFrameWidthX / kFrameWidthY margin.
//     Z thickness: kFrameHalfZ (covers both sub-layers of the view).
//     Material: configurable via setFrameMaterial(), default Aluminum.
//     Rotation: the frame is a daughter of the (already stereo-rotated)
//               layer envelope, so it rotates with the view.

#include "GeoModelKernel/GeoPhysVol.h"
#include "GeoModelKernel/GeoFullPhysVol.h"
#include <array>
#include <string>

class StrawTrackerBuilder {
public:
    // Returns the world physical volume; call once.
    GeoPhysVol* buildWorld();

    // ── Geometry constants ────────────────────────────────────────────────
    static constexpr int    kNStations    = 4;
    static constexpr int    kNLayers      = 4;   // per station (= views)
    static constexpr int    kNSubLayers   = 2;   // per layer (staggered pair)

    // Straw dimensions (mm)
    static constexpr double kStrawRadius  = 10.0;   // 1 cm radius → 2 cm diam
    static constexpr double kStrawLength  = 4000.0; // 4 m along X
    static constexpr double kWallThick    = 0.030;  // 30 µm Mylar

    // Station active face / aperture (mm). This is the INNER hole of the frame.
    static constexpr double kStationX     = 4000.0; // aperture width  (X, straw length)
    static constexpr double kStationY     = 6000.0; // aperture height (Y)

    // Straws per sub-layer
    static constexpr int kNStraws = static_cast<int>(kStationY / (2.0 * kStrawRadius)); // 300

    // Stereo angle magnitude (degrees)
    static constexpr double kStereoAngle  = 2.3;

    // Z positions of station centres (mm)
    // SHiP frame (mm). Station 0 at 84320; spacings preserved (2500/5000/1500).
    // Magnet centre sits at 89220, inside the 5000 mm gap between stations 1&2.
    static constexpr std::array<double,4> kStationZ = {84320., 86820., 91820., 93320.};

    /// SHiP-frame z of the Geant4/GeoModel world centre. World spans SHiP
    /// z in [0, 120000] mm, so the centre is at 60000. All internal (world)
    /// coordinates are SHiP - kWorldZOrigin; inputs/outputs are SHiP.
    static constexpr double kWorldZOrigin = 60000.0;
    static constexpr double kWorldHalfX   = 6000.0;    // 12 m x 12 m: must exceed
    static constexpr double kWorldHalfY   = 6000.0;    // the spectrometer volume
    static constexpr double kWorldHalfZ   = 60000.0;   // SHiP z 0 .. 120 m
    static constexpr double kMagnetZ      = 89220.0;   // SHiP frame

    // ---- Spectrometer volume ------------------------------------------------
    // A single 20 m-long volume, centred on the magnet, that contains BOTH the
    // four straw stations AND the whole field map. The field is attached here
    // (propagated to daughters), and the map's own bounds return B = 0 outside
    // itself -- so within this volume the field IS the map, and zero elsewhere.
    // x/y match the field map exactly (+-5000 mm), which is why the world had
    // to grow to +-6000 mm.
    //
    // Rationale: an undersized field volume silently TRUNCATES the field
    // integral (a 3.5 m box against a 15 m map cost us a 50% momentum bias).
    // Sizing the volume to the map, and containing the stations inside it,
    // means a station moved outside the mapped region fails LOUDLY in Geant4
    // rather than quietly reconstructing at the wrong momentum.
    static constexpr double kSpecHalfX    = 5000.0;    // = field map extent
    static constexpr double kSpecHalfY    = 5000.0;    // = field map extent
    static constexpr double kSpecHalfZ    = 10000.0;   // 20 m long
    static constexpr double kSpecZ        = kMagnetZ;  // centred on the magnet

    /// Dump every straw wire (WORLD frame = SHiP - kWorldZOrigin, mm) to a
    /// "Straws" TTree, for the ACTS reconstruction to build StrawSurfaces from.
    /// Standalone: computed from the constants above, no built detector needed.
    static void dumpStrawTable(const std::string& outFile,
                               double worldZOriginMM = kWorldZOrigin);


    // ── Frame geometry (per view) ─────────────────────────────────────────
    // Frame material width beyond the aperture, in X and Y (mm).
    // These are the half-thicknesses of the frame "bars".
    static constexpr double kFrameWidthX  = 100.0;  // 10 cm frame in X
    static constexpr double kFrameWidthY  = 100.0;  // 10 cm frame in Y

    // Frame thickness (half-length) along Z, must cover both sub-layers.
    // Sub-layers sit at z = ±kStrawRadius with half-Z ≈ kStrawRadius, so the
    // full envelope of straw material is roughly ±2*kStrawRadius = ±20 mm.
    static constexpr double kFrameHalfZ   = 22.0;   // covers [-22, +22] mm

    // Small clearance so daughters never touch the frame surface.
    static constexpr double kFrameClearance = 0.5;  // mm

    // Frame material: configurable, default Aluminum. Must match one of
    // MaterialManager::frameMaterialByName() supported names.
    void setFrameMaterial(const std::string& name) { m_frameMaterialName = name; }
    const std::string& frameMaterial() const { return m_frameMaterialName; }

    void writeDB(GeoPhysVol* world, const std::string& filename = "StrawTracker.db");

private:
    // ── Internal builders ─────────────────────────────────────────────────

    /// Build one station volume and place its 4 layers inside it.
    GeoPhysVol* buildStation(int stationID);

    /// Build one stereo layer (frame + two sub-layers + all straws).
    /// @param stereoAngleDeg  signed rotation about the beam axis Z
    GeoPhysVol* buildLayer(int layerID, double stereoAngleDeg);

    /// Build the material frame for one view (hollow rectangle).
    GeoPhysVol* buildFrame(int layerID);

    /// Build one sub-layer of N parallel straws.
    /// @param shifted  if true, this sub-layer is the staggered one
    GeoPhysVol* buildSubLayer(bool shifted);

    /// Build a single straw (wall tube + gas tube).
    GeoPhysVol* buildStraw(int uid);

    // ── State ─────────────────────────────────────────────────────────────
    std::string m_frameMaterialName = "Aluminum";
};
