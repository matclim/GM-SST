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

    constexpr double worldZOriginMM = 60000.;   // SHiP frame: world centre
    gun->SetParticlePosition(G4ThreeVector(
        m_gunCfg.posX * mm,
        m_gunCfg.posY * mm,
        (m_gunCfg.posZ - worldZOriginMM) * mm
    ));
    gun->SetParticleMomentumDirection(
        G4ThreeVector(m_gunCfg.dirX, m_gunCfg.dirY, m_gunCfg.dirZ).unit()
    );

    // LLP mode overrides the gun entirely: every charged daughter is fired from
    // the LLP's recorded decay vertex. The gun settings above are then ignored.
    if (!m_gunCfg.llpFile.empty())
        gen->openLLPFile(m_gunCfg.llpFile, worldZOriginMM);

    SetUserAction(gen);

    // Event and run actions
    auto* evtAction = new TrackerEventAction();
    SetUserAction(evtAction);
    SetUserAction(new TrackerRunAction(evtAction, m_outFile));
}
