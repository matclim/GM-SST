#pragma once
// TrackerDetectorConstruction.h
// Geant4 detector construction: converts the GeoModel tree into G4 volumes.

#include "G4VUserDetectorConstruction.hh"
#include "GeoModel2G4/ExtParameterisedVolumeBuilder.h"

class StrawTrackerBuilder;

class TrackerDetectorConstruction : public G4VUserDetectorConstruction {
public:
    TrackerDetectorConstruction();
    ~TrackerDetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;
    void ConstructSDandField() override;
};
