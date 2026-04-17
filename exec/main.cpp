// exec/main.cpp
// Entry point for the straw tracker Geant4 simulation.
//
// Options:
//   --n-events       <N>           Events to simulate            (default: 0)
//   --output         <file>        ROOT output file name         (default: StrawTracker_hits.root)
//   --seed           <N>           Random seed (0 = auto)        (default: 0)
//   --particle       <n>           Geant4 particle name          (default: mu-)
//   --energy-MeV     <E>           Kinetic energy [MeV]          (default: 10000)
//   --pos-mm         <x> <y> <z>   Gun position [mm] lab frame   (default: 0 0 24000)
//   --dir            <x> <y> <z>   Direction unit vector         (default: 0 0 1)
//   --field-map      <file>        Magnetic field-map text file
//                                  (empty = uniform -0.15 T fallback)
//   --frame-material <n>           Frame material (Aluminum, Mylar, Air)
//                                                                 (default: Aluminum)
//   --visualize                    Open interactive viewer
//   --vis-macro      <file>        Vis macro                     (default: straw_vis.mac)
//   --write-gdml                   Export GDML
//   --gdml-out       <file>        GDML file name                (default: StrawTracker_geometry.gdml)

#include "TrackerDetectorConstruction.h"
#include "TrackerActionInitialization.h"

#include "G4RunManagerFactory.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4TransportationManager.hh"
#include "FTFP_BERT.hh"
#include "G4GDMLParser.hh"
#include "G4SystemOfUnits.hh"
#include "CLHEP/Random/Random.h"
#include "CLHEP/Random/RanecuEngine.h"
#include "TApplication.h"
#include "TFile.h"

#include <iostream>
#include <string>
#include <stdexcept>

// ── Config ────────────────────────────────────────────────────────────────────
struct Config {
    int         nEvents   {0};
    std::string outFile   {"StrawTracker_hits.root"};
    long        seed      {0};
    bool        visualize {false};
    std::string visMacro  {"straw_vis.mac"};
    bool        writeGDML {false};
    std::string gdmlOut   {"StrawTracker_geometry.gdml"};
    // Gun settings
    std::string particle  {"mu-"};
    double      energyMeV {10000.};
    double      posX      {0.}, posY{0.}, posZ{24000.};
    double      dirX      {0.}, dirY{0.}, dirZ{1.};
    bool        writeDB   {false};
    std::string dbOut     {"StrawTracker.db"};
    // New — field map and frame material
    std::string fieldMap      {""};
    std::string frameMaterial {"Aluminum"};
};

static Config parseArgs(int argc, char** argv) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&](const char* flag) -> std::string {
            if (++i >= argc)
                throw std::runtime_error(std::string(flag) + " requires argument");
            return argv[i];
        };
        if      (a == "--n-events")        cfg.nEvents       = std::stoi(next(a.c_str()));
        else if (a == "--output")          cfg.outFile       = next(a.c_str());
        else if (a == "--seed")            cfg.seed          = std::stol(next(a.c_str()));
        else if (a == "--particle")        cfg.particle      = next(a.c_str());
        else if (a == "--energy-MeV")      cfg.energyMeV     = std::stod(next(a.c_str()));
        else if (a == "--vis-macro")       cfg.visMacro      = next(a.c_str());
        else if (a == "--gdml-out")        cfg.gdmlOut       = next(a.c_str());
        else if (a == "--visualize")       cfg.visualize     = true;
        else if (a == "--write-gdml")      cfg.writeGDML     = true;
        else if (a == "--write-db")        cfg.writeDB       = true;
        else if (a == "--db-out")          cfg.dbOut         = next(a.c_str());
        else if (a == "--field-map")       cfg.fieldMap      = next(a.c_str());
        else if (a == "--frame-material")  cfg.frameMaterial = next(a.c_str());
        else if (a == "--pos-mm") {
            cfg.posX = std::stod(next(a.c_str()));
            cfg.posY = std::stod(next(a.c_str()));
            cfg.posZ = std::stod(next(a.c_str()));
        } else if (a == "--dir") {
            cfg.dirX = std::stod(next(a.c_str()));
            cfg.dirY = std::stod(next(a.c_str()));
            cfg.dirZ = std::stod(next(a.c_str()));
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

    // Random seed
    if (cfg.seed != 0) {
        CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
        CLHEP::HepRandom::setTheSeed(cfg.seed);
    }

    TApplication rootApp("rootApp", nullptr, nullptr);
    {
        TFile* dummy = TFile::Open("/dev/null", "READ");
        if (dummy) delete dummy;
    }

    // Run manager
    auto* runMgr = G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly);
    runMgr->SetUserInitialization(
        new TrackerDetectorConstruction(cfg.fieldMap, cfg.frameMaterial));
    runMgr->SetUserInitialization(new FTFP_BERT());

    // Build gun config and pass to action initialisation
    TrackerActionInitialization::GunConfig gunCfg;
    gunCfg.particle  = cfg.particle;
    gunCfg.energyMeV = cfg.energyMeV;
    gunCfg.posX      = cfg.posX;
    gunCfg.posY      = cfg.posY;
    gunCfg.posZ      = cfg.posZ;
    gunCfg.dirX      = cfg.dirX;
    gunCfg.dirY      = cfg.dirY;
    gunCfg.dirZ      = cfg.dirZ;

    runMgr->SetUserInitialization(new TrackerActionInitialization(cfg.outFile, gunCfg));

    runMgr->Initialize();

    // Export GDML
    if (cfg.writeGDML) {
        G4GDMLParser parser;
        parser.Write(cfg.gdmlOut,
                     G4TransportationManager::GetTransportationManager()
                         ->GetNavigatorForTracking()->GetWorldVolume());
        std::cout << "[main] GDML written to " << cfg.gdmlOut << "\n";
    }

    // Visualisation or batch run
    if (cfg.visualize) {
        auto* ui  = new G4UIExecutive(argc, argv);
        auto* vis = new G4VisExecutive();
        vis->Initialize();

        auto* uiMgr = G4UImanager::GetUIpointer();
        uiMgr->ApplyCommand("/control/execute " + cfg.visMacro);
        if (cfg.nEvents > 0)
            uiMgr->ApplyCommand("/run/beamOn " + std::to_string(cfg.nEvents));

        ui->SessionStart();
        delete ui;
        delete vis;
    } else if (cfg.nEvents > 0) {
        runMgr->BeamOn(cfg.nEvents);
    }

    delete runMgr;
    return 0;
}
