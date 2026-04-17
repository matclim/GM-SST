#pragma once
// TrackerDetectorConstruction.h
// Bridges the GeoModel tree into Geant4 via GeoModel2G4, and attaches the
// magnetic field region.

#include "G4VUserDetectorConstruction.hh"

#include <string>

class TrackerDetectorConstruction : public G4VUserDetectorConstruction {
public:
    /// @param fieldMapFile    optional path to a 3-D field-map file. Empty
    ///                        string uses the uniform -0.15 T fallback.
    /// @param frameMaterial   material name for the view frames. Supported
    ///                        values: Aluminum / Aluminium / Al / Mylar / Air.
    ///                        Default: Aluminum.
    TrackerDetectorConstruction(std::string fieldMapFile = "",
                                std::string frameMaterial = "Aluminum");

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;

private:
    std::string m_fieldMapFile;
    std::string m_frameMaterial;
};
