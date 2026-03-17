// TrackerDetectorConstruction.cpp
// Bridges the GeoModel tree into Geant4 via GeoModel2G4.

#include "TrackerDetectorConstruction.h"
#include "StrawTrackerBuilder.h"
#include "TrackerSD.h"
#include "MaterialManager.h"

#include "GeoModel2G4/ExtParameterisedVolumeBuilder.h"
#include "GeoModelKernel/GeoPhysVol.h"
#include "G4PVPlacement.hh"

#include "G4SDManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"

// Global pointer to the SD so TrackerEventAction can reach it.
TrackerSD* g_strawSD = nullptr;

TrackerDetectorConstruction::TrackerDetectorConstruction()
    : G4VUserDetectorConstruction()
{}

G4VPhysicalVolume* TrackerDetectorConstruction::Construct() {
    // 1. Build geometry using GeoModel
    StrawTrackerBuilder builder;
    GeoPhysVol* geoWorld = builder.buildWorld();

    // 2. Convert GeoModel tree → Geant4 volume tree
    ExtParameterisedVolumeBuilder g4Builder("StrawTracker");

    G4LogicalVolume* g4WorldLog = g4Builder.Build(geoWorld);

    G4VPhysicalVolume* g4World = new G4PVPlacement(
        nullptr,             // no rotation
        G4ThreeVector(),     // at origin
        g4WorldLog,
        "World",
        nullptr,             // no mother volume (this is the world)
        false, 0, false);

    return g4World;


}

void TrackerDetectorConstruction::ConstructSDandField() {
    // Attach the sensitive detector to every logical volume named "StrawGas".
    auto* sdManager = G4SDManager::GetSDMpointer();

    auto* sd = new TrackerSD("StrawTrackerSD");
    sdManager->AddNewDetector(sd);
    g_strawSD = sd;

    // Walk the logical volume store and attach SD to all "StrawGas" volumes.
    auto* lvStore = G4LogicalVolumeStore::GetInstance();
    for (auto* lv : *lvStore) {
        if (lv->GetName() == "StrawGas") {
            lv->SetSensitiveDetector(sd);
        }
    }
}
