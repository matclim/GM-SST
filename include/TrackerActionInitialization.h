#pragma once
// TrackerActionInitialization.h
// Required by Geant4 MT: all user actions must be registered here.

#include "G4VUserActionInitialization.hh"
#include <string>

class TrackerActionInitialization : public G4VUserActionInitialization {
public:
    // gunCfg holds the CLI particle/energy/position settings so each
    // worker thread can configure its own gun identically.
    struct GunConfig {
        std::string particle     {"mu-"};
        double      energyMeV    {10000.};
        double      posX         {0.};
        double      posY         {0.};
        double      posZ         {24000.};   // lab-frame mm
        double      dirX         {0.};
        double      dirY         {0.};
        double      dirZ         {1.};
    };

    TrackerActionInitialization(const std::string& outFile,
                                const GunConfig&   gunCfg)
        : m_outFile(outFile), m_gunCfg(gunCfg) {}

    // Called once on the master thread.
    void BuildForMaster() const override;

    // Called once per worker thread.
    void Build() const override;

private:
    std::string m_outFile;
    GunConfig   m_gunCfg;
};
