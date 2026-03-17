// TrackerPrimaryGeneratorAction.cpp
// Default primary generator: 10 GeV muon along +Z, aimed at the tracker centre.

#include "TrackerPrimaryGeneratorAction.h"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Event.hh"

TrackerPrimaryGeneratorAction::TrackerPrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction()
    , m_gun(std::make_unique<G4ParticleGun>(1))
{
    auto* mu = G4ParticleTable::GetParticleTable()->FindParticle("mu-");
    m_gun->SetParticleDefinition(mu);
    m_gun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    m_gun->SetParticleEnergy(10000. * MeV);  // 10 GeV

    // Place gun upstream of the first station (z = 26500 mm → gun at 24000 mm).
    // World is centred at z=31000 mm in lab → local gun z = 24000-31000 = -7000 mm
    m_gun->SetParticlePosition(G4ThreeVector(0., 0., -7000. * mm));
}

void TrackerPrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    m_gun->GeneratePrimaryVertex(event);
}
