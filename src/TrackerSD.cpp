// TrackerSD.cpp
// Records energy deposits in straw gas volumes.
// The touchable history is used to decode:
//   depth 0: StrawGas     (the sensitive volume itself)
//   depth 1: StrawWall    (the straw physical volume)
//   depth 2: SubLayer
//   depth 3: Layer
//   depth 4: Station
//
// The copy numbers set via GeoIdentifierTag map to straw/sub-layer/layer/station IDs.

#include "TrackerSD.h"
#include <cmath>

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SystemOfUnits.hh"

TrackerSD::TrackerSD(const G4String& name)
    : G4VSensitiveDetector(name)
{}

void TrackerSD::Initialize(G4HCofThisEvent*) {
    m_hits.clear();
}

G4bool TrackerSD::ProcessHits(G4Step* step, G4TouchableHistory*) {
    const double edep = step->GetTotalEnergyDeposit();
    if (edep <= 0.) return false;

    auto* touch = step->GetPreStepPoint()->GetTouchable();

    static bool printed = false;
    if (!printed) {
        printed = true;
        G4cout << "[TrackerSD] Touchable history depth: "
               << touch->GetHistoryDepth() << G4endl;
        for (int d = 0; d <= touch->GetHistoryDepth(); ++d) {
            G4cout << "  depth " << d << ": "
                   << touch->GetVolume(d)->GetName()
                   << "  copyNo=" << touch->GetCopyNumber(d) << G4endl;
        }
    }

    // Decode copy numbers from the touchable history.
    // Depth 0 = StrawGas, 1 = StrawWall (straw copy number),
    // 2 = SubLayer, 3 = StrawLayer, 4 = Station
    const int strawID    = (touch->GetHistoryDepth() >= 1) ? touch->GetCopyNumber(1) : -1;
    const int subLayerID = (touch->GetHistoryDepth() >= 2) ? touch->GetCopyNumber(2) : -1;
    const int layerID    = (touch->GetHistoryDepth() >= 3) ? touch->GetCopyNumber(3) : -1;
    const int stationID  = (touch->GetHistoryDepth() >= 4) ? touch->GetCopyNumber(4) : -1;

    // Step centroid (average of pre- and post-step positions)
    const G4ThreeVector pre  = step->GetPreStepPoint()->GetPosition();
    const G4ThreeVector post = step->GetPostStepPoint()->GetPosition();
    const G4ThreeVector mid  = 0.5 * (pre + post);

    const G4Track* trk = step->GetTrack();
    const G4ThreeVector vtx  = trk->GetVertexPosition();
    // truth momentum at production = |p| * direction at vertex
    const G4ThreeVector vdir = trk->GetVertexMomentumDirection();
    const double m0   = trk->GetDefinition()->GetPDGMass();
    const double eKin = trk->GetVertexKineticEnergy();
    const double pMag = std::sqrt(eKin * (eKin + 2.0 * m0));
    StrawHit hit;
    hit.trackID    = trk->GetTrackID();
    hit.parentID   = trk->GetParentID();
    hit.pdg        = trk->GetDefinition()->GetPDGEncoding();
    hit.vtxX       = vtx.x() / mm;
    hit.vtxY       = vtx.y() / mm;
    hit.vtxZ       = vtx.z() / mm;
    hit.vpx        = pMag * vdir.x() / MeV;
    hit.vpy        = pMag * vdir.y() / MeV;
    hit.vpz        = pMag * vdir.z() / MeV;
    hit.stationID  = stationID;
    hit.layerID    = layerID;
    hit.subLayerID = subLayerID;
    hit.strawID    = strawID;
    hit.edep       = edep / MeV;
    hit.x          = mid.x()  / mm;
    hit.y          = mid.y()  / mm;
    hit.z          = mid.z()  / mm;
    hit.x_entry    = pre.x()  / mm;
    hit.y_entry    = pre.y()  / mm;
    hit.z_entry    = pre.z()  / mm;
    hit.x_exit     = post.x() / mm;
    hit.y_exit     = post.y() / mm;
    hit.z_exit     = post.z() / mm;

    m_hits.push_back(hit);
    return true;
}

void TrackerSD::EndOfEvent(G4HCofThisEvent*) {
    // Hits are read by TrackerEventAction; nothing to do here.
}
