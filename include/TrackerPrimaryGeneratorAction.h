#pragma once
// TrackerPrimaryGeneratorAction.h

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include <memory>

class G4Event;

class TrackerPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    TrackerPrimaryGeneratorAction();
    ~TrackerPrimaryGeneratorAction() override = default;

    void GeneratePrimaries(G4Event*) override;

    G4ParticleGun* gun() { return m_gun.get(); }

private:
    std::unique_ptr<G4ParticleGun> m_gun;
};
