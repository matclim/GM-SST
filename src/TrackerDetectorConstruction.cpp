// TrackerDetectorConstruction.cpp
// Bridges the GeoModel tree into Geant4 via GeoModel2G4.

#include "G4PVPlacement.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4SDManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SystemOfUnits.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "MagneticFieldRegion.h"

#include "GeoModel2G4/ExtParameterisedVolumeBuilder.h"
#include "GeoModelKernel/GeoPhysVol.h"

#include "TrackerDetectorConstruction.h"
#include "StrawTrackerBuilder.h"
#include "TrackerSD.h"

#include <iostream>
#include <utility>

TrackerDetectorConstruction::TrackerDetectorConstruction(std::string fieldMapFile,
                                                         std::string frameMaterial)
    : G4VUserDetectorConstruction(),
      m_fieldMapFile(std::move(fieldMapFile)),
      m_frameMaterial(std::move(frameMaterial))
{}

G4VPhysicalVolume* TrackerDetectorConstruction::Construct() {
    StrawTrackerBuilder builder;
    builder.setFrameMaterial(m_frameMaterial);

    GeoPhysVol* geoWorld = builder.buildWorld();
    builder.writeDB(geoWorld, "StrawTracker.db");

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

        } else if (name.substr(0, 12) == "StrawFrameLV") {
            // Frame visualisation: solid light-gray, clearly visible.
            auto* va = new G4VisAttributes(G4Colour(0.7, 0.7, 0.75, 0.5));
            va->SetForceSolid(true);
            lv->SetVisAttributes(va);

        } else if (name.substr(0, 9) == "StrawWall") {
            auto* va = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.3));
            va->SetForceSolid(true);
            lv->SetVisAttributes(va);

        } else if (name.substr(0, 8) == "StrawGas") {
            auto* va = new G4VisAttributes(G4Colour(0.0, 0.8, 1.0, 0.4));
            va->SetForceSolid(true);
            lv->SetVisAttributes(va);
        }
    }

    // ── Magnetic field region ─────────────────────────────────────────────
    constexpr double worldZOriginMM = 31000.;
    MagneticFieldRegion::build(g4WorldLog, worldZOriginMM, m_fieldMapFile);

    G4VPhysicalVolume* g4World = new G4PVPlacement(
        nullptr, G4ThreeVector(), g4WorldLog,
        "World", nullptr, false, 0, /*pSurfChk=*/true);

    // ── Sanity-check the geometry for sibling overlaps ────────────────────
    // This walks every placed physical volume and runs Geant4's point-sampling
    // overlap test. It's a one-off cost at startup and prints a clear warning
    // for each overlap detected, which is invaluable for diagnosing "no hits"
    // symptoms caused by undefined navigation through overlapping siblings.
    {
        auto* pvStore = G4PhysicalVolumeStore::GetInstance();
        std::cout << "[TrackerDetectorConstruction] Checking "
                  << pvStore->size() << " placed volumes for overlaps...\n";
        int nWithOverlaps = 0;
        for (auto* pv : *pvStore) {
            if (pv->CheckOverlaps(200, 0.0, /*verbose=*/false))
                ++nWithOverlaps;
        }
        if (nWithOverlaps == 0)
            std::cout << "[TrackerDetectorConstruction] "
                      << "No overlaps detected.\n";
        else
            std::cout << "[TrackerDetectorConstruction] WARNING: "
                      << nWithOverlaps << " placements have overlaps. "
                      << "Re-run with G4's verbose=true for details.\n";
    }

    return g4World;
}


void TrackerDetectorConstruction::ConstructSDandField() {
    // Attach the sensitive detector to every logical volume named "StrawGas*".
    auto* sdManager = G4SDManager::GetSDMpointer();

    auto* sd = new TrackerSD("StrawTrackerSD");
    sdManager->AddNewDetector(sd);

    auto* lvStore = G4LogicalVolumeStore::GetInstance();
    for (auto* lv : *lvStore) {
        if (lv->GetName().substr(0, 8) == "StrawGas")
            lv->SetSensitiveDetector(sd);
    }
}
