#pragma once
// MagneticFieldRegion.h
// Dipole magnetic field volume placed between tracker stations 1 and 2.
// Can use either a uniform field (fallback) or a 3-D field map from file.

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4SystemOfUnits.hh"

#include <string>

class G4LogicalVolume;

class MagneticFieldRegion {
public:
    /// Build the field volume and attach a field manager to it.
    ///
    /// @param worldLog         the world logical volume (field volume is
    ///                         placed as a daughter of this)
    /// @param worldZOrigin_mm  lab-frame Z of the world centre (31000 mm here)
    /// @param fieldMapFile     optional path to a text field-map file. If the
    ///                         string is empty or the file cannot be loaded,
    ///                         falls back to the original uniform dipole
    ///                         (-0.15 T along X).
    static void build(G4LogicalVolume* worldLog,
                      double worldZOrigin_mm,
                      const std::string& fieldMapFile = "");
};
