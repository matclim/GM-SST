#include "MagneticFieldRegion.h"
#include "ShipFieldMap.h"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4TransportationManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"

#include <iostream>

void MagneticFieldRegion::build(G4LogicalVolume* worldLog,
                                double worldZOrigin_mm,
                                const std::string& fieldMapFile) {
    // ── Geometry ──────────────────────────────────────────────────────────────
    // Field volume sits between tracker stations 1 (z=29000) and 2 (z=34000).
    // Its lab-frame centre is at z = 31500 mm. Dimensions match the tracker
    // active area (4 m × 6 m) and the magnet gap (3.5 m).
    const double halfX =  2000.0 * mm;   // 4 m
    const double halfY =  3000.0 * mm;   // 6 m
    const double halfZ =  1750.0 * mm;   // 3.5 m

    const double labZ   = 31500.0 * mm;
    const double localZ = labZ - worldZOrigin_mm * mm;

    auto* airMat = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");

    auto* fieldBox = new G4Box("MagnetVolume", halfX, halfY, halfZ);
    auto* fieldLog = new G4LogicalVolume(fieldBox, airMat, "MagnetVolume");

    // Visualise as a semi-transparent blue box
    auto* va = new G4VisAttributes(G4Colour(0.2, 0.2, 1.0, 0.15));
    va->SetForceSolid(true);
    fieldLog->SetVisAttributes(va);

    new G4PVPlacement(nullptr,
                      G4ThreeVector(0., 0., localZ),
                      fieldLog,
                      "MagnetVolume",
                      worldLog,
                      false, 0, false);

    // ── Pick a field implementation ───────────────────────────────────────────
    G4MagneticField* magField = nullptr;

    if (!fieldMapFile.empty()) {
        // Origin of the field map in the Geant4 world frame. The world is
        // centred at lab z = worldZOrigin_mm, and by the user's convention
        // the field-map z = 0 corresponds to lab z = 31500 mm. So in the
        // Geant4 frame the map origin sits at z = (31500 − worldZOrigin_mm).
        const G4ThreeVector mapOriginG4(
            0.0,
            0.0,
            (31500.0 - worldZOrigin_mm) * mm
        );

        auto* mapField = new ShipFieldMap(fieldMapFile, mapOriginG4);
        if (mapField->valid()) {
            magField = mapField;
            std::cout << "[MagneticFieldRegion] Using field map from '"
                      << fieldMapFile << "'\n";
        } else {
            std::cerr << "[MagneticFieldRegion] Field map invalid, falling "
                      << "back to uniform dipole.\n";
            delete mapField;
        }
    }

    if (!magField) {
        // Original fallback: uniform -0.15 T along X
        magField = new G4UniformMagField(G4ThreeVector(-0.15 * tesla, 0., 0.));
        std::cout << "[MagneticFieldRegion] Using uniform dipole "
                  << "(-0.15 T along X)\n";
    }

    // ── Field manager plumbing ────────────────────────────────────────────────
    auto* eqRhs       = new G4Mag_UsualEqRhs(magField);
    auto* stepper     = new G4ClassicalRK4(eqRhs);
    auto* chordFinder = new G4ChordFinder(magField, 1.0 * mm, stepper);

    auto* fieldMgr = new G4FieldManager(magField, chordFinder);
    fieldMgr->SetDetectorField(magField);

    // Attach to this volume (not propagated to daughters — we have none here).
    fieldLog->SetFieldManager(fieldMgr, true);
}
