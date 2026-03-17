// TrackerDetectorConstruction.cpp
// Bridges the GeoModel tree into Geant4 via GeoModel2G4.

#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"


#include "GeoModel2G4/ExtParameterisedVolumeBuilder.h"
#include "GeoModelKernel/GeoPhysVol.h"

#include "TrackerDetectorConstruction.h"
#include "StrawTrackerBuilder.h"
#include "TrackerSD.h"

TrackerDetectorConstruction::TrackerDetectorConstruction()
    : G4VUserDetectorConstruction()
{}

G4VPhysicalVolume* TrackerDetectorConstruction::Construct() {
    StrawTrackerBuilder builder;
    GeoPhysVol* geoWorld = builder.buildWorld();

    ExtParameterisedVolumeBuilder g4Builder("StrawTracker");
    G4LogicalVolume* g4WorldLog = g4Builder.Build(geoWorld);

    // ── Visibility attributes ─────────────────────────────────────────────
    auto* lvStore = G4LogicalVolumeStore::GetInstance();
    for (auto* lv : *lvStore) {
        const G4String& name = lv->GetName();

        if (name == "World") {
            lv->SetVisAttributes(G4VisAttributes::GetInvisible());

        } else if (name == "Station") {
            auto* va = new G4VisAttributes(false); // invisible envelope
            lv->SetVisAttributes(va);

        } else if (name == "StrawLayer") {
            auto* va = new G4VisAttributes(false); // invisible envelope
            lv->SetVisAttributes(va);

        } else if (name == "SubLayer_nominal" || name == "SubLayer_shifted") {
            auto* va = new G4VisAttributes(false); // invisible envelope
            lv->SetVisAttributes(va);

        } else if (name == "StrawWall") {
            auto* va = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.3)); // grey, translucent
            va->SetForceSolid(true);
            lv->SetVisAttributes(va);

        } else if (name == "StrawGas") {
            auto* va = new G4VisAttributes(G4Colour(0.0, 0.8, 1.0, 0.4)); // cyan, translucent
            va->SetForceSolid(true);
            lv->SetVisAttributes(va);
        }
    }

    G4VPhysicalVolume* g4World = new G4PVPlacement(
        nullptr, G4ThreeVector(), g4WorldLog,
        "World", nullptr, false, 0, false);

    return g4World;
}


void TrackerDetectorConstruction::ConstructSDandField() {
    // Attach the sensitive detector to every logical volume named "StrawGas".
    // This is called once per worker thread in MT mode; each worker gets its
    // own TrackerSD instance, which is what we want.
    auto* sdManager = G4SDManager::GetSDMpointer();

    auto* sd = new TrackerSD("StrawTrackerSD");
    sdManager->AddNewDetector(sd);

    auto* lvStore = G4LogicalVolumeStore::GetInstance();
    for (auto* lv : *lvStore) {
        if (lv->GetName() == "StrawGas")
            lv->SetSensitiveDetector(sd);
    }
}
