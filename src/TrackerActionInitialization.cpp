// TrackerActionInitialization.cpp
#include "TrackerActionInitialization.h"
#include "TrackerRunAction.h"
#include "TrackerEventAction.h"
#include "TrackerPrimaryGeneratorAction.h"

#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"

void TrackerActionInitialization::BuildForMaster() const {
    // Master thread only needs a run action to open/close the ROOT file.
    SetUserAction(new TrackerRunAction(nullptr, m_outFile));
}

void TrackerActionInitialization::Build() const {
    // Primary generator
    auto* gen = new TrackerPrimaryGeneratorAction();

    // Configure gun from stored CLI settings.
    // Particle table is safe to query here; workers initialise after master.
    auto* gun  = gen->gun();
    auto* ptbl = G4ParticleTable::GetParticleTable();
    auto* pdef = ptbl->FindParticle(m_gunCfg.particle);
    if (pdef) gun->SetParticleDefinition(pdef);

    gun->SetParticleEnergy(m_gunCfg.energyMeV * MeV);

    constexpr double worldZOriginMM = 31000.;
    gun->SetParticlePosition(G4ThreeVector(
        m_gunCfg.posX * mm,
        m_gunCfg.posY * mm,
        (m_gunCfg.posZ - worldZOriginMM) * mm
    ));
    gun->SetParticleMomentumDirection(
        G4ThreeVector(m_gunCfg.dirX, m_gunCfg.dirY, m_gunCfg.dirZ).unit()
    );

    SetUserAction(gen);

    // Event and run actions
    auto* evtAction = new TrackerEventAction();
    SetUserAction(evtAction);
    SetUserAction(new TrackerRunAction(evtAction, m_outFile));
}
