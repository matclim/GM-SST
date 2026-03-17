#pragma once
// TrackerRunAction.h
// In Geant4 MT each worker thread gets its own TrackerRunAction instance
// and writes to its own per-thread ROOT file. The master run action does
// nothing with ROOT; a post-run hadd step merges the files if desired.

#include "G4UserRunAction.hh"
#include <string>

class TFile;
class TTree;
class TrackerEventAction;

class TrackerRunAction : public G4UserRunAction {
public:
    // evtAction is nullptr on the master thread.
    TrackerRunAction(TrackerEventAction* evtAction,
                     const std::string& outFile = "StrawTracker_hits.root");
    ~TrackerRunAction() override = default;

    void BeginOfRunAction(const G4Run*) override;
    void EndOfRunAction(const G4Run*)   override;

private:
    TrackerEventAction* m_evtAction {nullptr};
    std::string         m_outFile;
    TFile*              m_rootFile  {nullptr};
    TTree*              m_tree      {nullptr};
};
