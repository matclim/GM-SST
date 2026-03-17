#pragma once
// TrackerSD.h
// Sensitive detector that records energy deposits in straw gas volumes.

#include "G4VSensitiveDetector.hh"
#include "TrackerHit.h"
#include <vector>

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class TrackerSD : public G4VSensitiveDetector {
public:
    explicit TrackerSD(const G4String& name);
    ~TrackerSD() override = default;

    void Initialize(G4HCofThisEvent*) override;
    G4bool ProcessHits(G4Step*, G4TouchableHistory*) override;
    void EndOfEvent(G4HCofThisEvent*) override;

    const std::vector<StrawHit>& hits() const { return m_hits; }
    void clearHits() { m_hits.clear(); }

private:
    std::vector<StrawHit> m_hits;
};
