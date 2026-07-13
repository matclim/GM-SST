#include "MagneticFieldRegion.h"
#include "ShipFieldMap.h"

#include "G4LogicalVolume.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4TransportationManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UniformMagField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4ChordFinder.hh"
#include "G4FieldManager.hh"

#include <iostream>

// =============================================================================
// The field is attached to the SPECTROMETER volume (built in StrawTrackerBuilder),
// which is 20 m long and matches the field map in x/y (+-5000 mm), and which
// CONTAINS the four straw stations as daughters. The field is propagated to
// those daughters, so it acts inside the straws -- which is where the map says
// it is: the map spans +-7500 mm in z and therefore overlaps the whole tracker.
//
// ShipFieldMap::GetFieldValue returns B = 0 outside the map's own bounds, so
// within the spectrometer the field IS the map, and zero elsewhere.
//
// WHY: an undersized field box silently TRUNCATES the field integral. A 3.5 m
// box against a 15 m map produced a ~50% momentum over-estimate that took a
// long time to find, because every component was individually correct. Sizing
// the volume to the map, with the stations inside it, means a station moved
// beyond the mapped region fails LOUDLY in Geant4 rather than quietly
// reconstructing at the wrong momentum.
// =============================================================================
void MagneticFieldRegion::build(G4LogicalVolume* worldLog,
                                double worldZOrigin_mm,
                                const std::string& fieldMapFile) {
    // SHiP-frame z of the magnet (= the field map's own origin).
    constexpr double kMagnetZ_SHiP_mm = 89220.0;

    G4MagneticField* magField = nullptr;

    if (!fieldMapFile.empty()) {
        // Map origin in the Geant4 world frame: the world centre sits at SHiP
        // z = worldZOrigin_mm, and the map's z = 0 is the magnet centre.
        const G4ThreeVector mapOriginG4(
            0.0,
            0.0,
            (kMagnetZ_SHiP_mm - worldZOrigin_mm) * mm
        );

        auto* mapField = new ShipFieldMap(fieldMapFile, mapOriginG4);
        if (mapField->valid()) {
            magField = mapField;
            std::cout << "[MagneticFieldRegion] Using field map from '"
                      << fieldMapFile << "'\n"
                      << "[MagneticFieldRegion] Field attached to the WORLD; "
                      << "the map's own bounds define its extent.\n";
        } else {
            std::cerr << "[MagneticFieldRegion] Field map invalid, falling "
                      << "back to uniform dipole.\n";
            delete mapField;
        }
    }

    if (!magField) {
        // Fallback: uniform -0.15 T along X, over the whole world.
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

    // Attach to the SPECTROMETER volume, propagated to its daughters (the four
    // straw stations). Fall back to the world if it isn't there, so a geometry
    // without the spectrometer still gets a field -- but say so loudly.
    G4LogicalVolume* target =
        G4LogicalVolumeStore::GetInstance()->GetVolume("Spectrometer", false);
    if (target) {
        std::cout << "[MagneticFieldRegion] Field attached to 'Spectrometer' "
                  << "(20 m, matches the map in x/y); propagated to the stations.\n";
    } else {
        target = worldLog;
        std::cerr << "[MagneticFieldRegion] WARNING: no 'Spectrometer' volume "
                  << "found -- attaching the field to the World instead.\n";
    }
    target->SetFieldManager(fieldMgr, true);
}
