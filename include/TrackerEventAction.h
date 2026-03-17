#pragma once
// TrackerEventAction.h
// Collects hits from the SD at end of event and fills the TTree.

#include "G4UserEventAction.hh"
#include "TrackerHit.h"
#include <vector>

class TTree;
class TrackerSD;

class TrackerEventAction : public G4UserEventAction {
public:
    TrackerEventAction() = default;
    ~TrackerEventAction() override = default;

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*)   override;

    // Called by TrackerRunAction to bind branches.
    void setTree(TTree* tree);

    // Set the pointer to the SD (done after ConstructSDandField).
    void setSD(TrackerSD* sd) { m_sd = sd; }

private:
    TTree*     m_tree {nullptr};
    TrackerSD* m_sd   {nullptr};

    // Branch buffers
    std::vector<int>    m_trackID, m_stationID, m_layerID, m_subLayerID, m_strawID;
    std::vector<double> m_edep;
    std::vector<double> m_x,  m_y,  m_z;
    std::vector<double> m_xe, m_ye, m_ze;   // entry
    std::vector<double> m_xx, m_yx, m_zx;   // exit
};
