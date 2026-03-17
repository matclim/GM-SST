// exec/main.cpp
// Entry point for the straw tracker Geant4 simulation.
//
// Usage:
//   ./run_StrawTracker [options]
//
// Options:
//   --n-events   <N>          Events to simulate            (default: 0)
//   --output     <file>       ROOT output file name         (default: StrawTracker_hits.root)
//   --seed       <N>          Random seed (0 = auto)        (default: 0)
//   --particle   <name>       Geant4 particle name          (default: mu-)
//   --energy-MeV <E>          Kinetic energy [MeV]          (default: 10000)
//   --pos-mm     <x> <y> <z>  Gun position [mm] lab frame   (default: 0 0 24000)
//   --dir        <x> <y> <z>  Direction unit vector         (default: 0 0 1)
//   --sigma-xy-mm <s>         Gaussian beam spread in X,Y   (default: 0)
//   --visualize               Open interactive viewer
//   --vis-macro  <file>       Vis macro                     (default: straw_vis.mac)
//   --write-gdml              Export GDML
//   --gdml-out   <file>       GDML file name                (default: StrawTracker_geometry.gdml)

#include "TrackerDetectorConstruction.h"
#include "TrackerRunAction.h"
#include "TrackerEventAction.h"
#include "TrackerPrimaryGeneratorAction.h"
#include "TrackerActionInitialization.h"

// Geant4
#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "FTFP_BERT.hh"
#include "G4GDMLParser.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RanecuEngine.h"

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <stdexcept>

// Extern SD pointer set in TrackerDetectorConstruction.cpp
extern TrackerSD* g_strawSD;

// ── Simple argument parser ────────────────────────────────────────────────────
struct Config {
    int         nEvents      {0};
    std::string outFile      {"StrawTracker_hits.root"};
    long        seed         {0};
    std::string particle     {"mu-"};
    double      energyMeV    {10000.};
    double      posX{0.}, posY{0.}, posZ{24000.};   // mm, lab frame
    double      dirX{0.}, dirY{0.}, dirZ{1.};
    double      sigmaXY      {0.};
    bool        visualize    {false};
    std::string visMacro     {"straw_vis.mac"};
    bool        writeGDML    {false};
    std::string gdmlOut      {"StrawTracker_geometry.gdml"};
};

static Config parseArgs(int argc, char** argv) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&](const char* flag) -> std::string {
            if (++i >= argc) throw std::runtime_error(std::string(flag) + " requires argument");
            return argv[i];
        };
        if      (a == "--n-events")    cfg.nEvents      = std::stoi(next("--n-events"));
        else if (a == "--output")      cfg.outFile       = next("--output");
        else if (a == "--seed")        cfg.seed          = std::stol(next("--seed"));
        else if (a == "--particle")    cfg.particle      = next("--particle");
        else if (a == "--energy-MeV")  cfg.energyMeV     = std::stod(next("--energy-MeV"));
        else if (a == "--sigma-xy-mm") cfg.sigmaXY       = std::stod(next("--sigma-xy-mm"));
        else if (a == "--vis-macro")   cfg.visMacro      = next("--vis-macro");
        else if (a == "--gdml-out")    cfg.gdmlOut       = next("--gdml-out");
        else if (a == "--visualize")   cfg.visualize     = true;
        else if (a == "--write-gdml")  cfg.writeGDML     = true;
        else if (a == "--pos-mm") {
            cfg.posX = std::stod(next("--pos-mm"));
            cfg.posY = std::stod(next("--pos-mm"));
            cfg.posZ = std::stod(next("--pos-mm"));
        } else if (a == "--dir") {
            cfg.dirX = std::stod(next("--dir"));
            cfg.dirY = std::stod(next("--dir"));
            cfg.dirZ = std::stod(next("--dir"));
        } else {
            std::cerr << "[main] Unknown option: " << a << "\n";
        }
    }
    return cfg;
}

// ── main ──────────────────────────────────────────────────────────────────────
int main(int argc, char** argv) {
    Config cfg;
    try { cfg = parseArgs(argc, argv); }
    catch (const std::exception& e) {
        std::cerr << "Argument error: " << e.what() << "\n";
        return 1;
    }

    // ── Random seed ───────────────────────────────────────────────────────────
    if (cfg.seed != 0) {
      
        CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
        CLHEP::HepRandom::setTheSeed(cfg.seed);
    }

    // ── Run manager ───────────────────────────────────────────────────────────
    auto* runMgr = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Default);

    // Detector
    runMgr->SetUserInitialization(new TrackerDetectorConstruction());

    // Physics list
    runMgr->SetUserInitialization(new FTFP_BERT());

    // User actions
    runMgr->SetUserInitialization(new TrackerActionInitialization(cfg.outFile));
    // Configure gun from CLI
    {
        
    auto* genAction = static_cast<const TrackerPrimaryGeneratorAction*>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
    auto* gun = genAction->gun();


    }
        gun->SetParticleDefinition(pdef);
        gun->SetParticleEnergy(cfg.energyMeV * MeV);

        // Convert lab-frame gun position to world-local frame:
        // world origin is at z=31000 mm in lab frame.
        constexpr double worldZOriginMM = 31000.;
        gun->SetParticlePosition(G4ThreeVector(
            cfg.posX * mm,
            cfg.posY * mm,
            (cfg.posZ - worldZOriginMM) * mm
        ));
        gun->SetParticleMomentumDirection(
            G4ThreeVector(cfg.dirX, cfg.dirY, cfg.dirZ).unit()
        );
    }

    runMgr->Initialize();

    // Wire SD into event action after geometry is initialised
    evtAction->setSD(g_strawSD);

    // ── Export GDML ───────────────────────────────────────────────────────────
    if (cfg.writeGDML) {
        G4GDMLParser parser;
        parser.Write(cfg.gdmlOut,
                     G4TransportationManager::GetTransportationManager()
                         ->GetNavigatorForTracking()
                         ->GetWorldVolume());
        std::cout << "[main] GDML written to " << cfg.gdmlOut << "\n";
    }

    // ── Visualisation ─────────────────────────────────────────────────────────
    if (cfg.visualize) {
        G4UIExecutive* ui = new G4UIExecutive(argc, argv);
        G4VisManager* vis = new G4VisExecutive();
        vis->Initialize();

        auto* uiMgr = G4UImanager::GetUIpointer();
        uiMgr->ApplyCommand("/control/execute " + cfg.visMacro);

        if (cfg.nEvents > 0) {
            uiMgr->ApplyCommand("/run/beamOn " + std::to_string(cfg.nEvents));
        }

        ui->SessionStart();
        delete ui;
        delete vis;
    } else if (cfg.nEvents > 0) {
        // Batch mode
        runMgr->BeamOn(cfg.nEvents);
    }

    delete runMgr;
    return 0;
}
