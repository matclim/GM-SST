#pragma once
// TrackerRunAction.h
// Opens / closes the ROOT output file and owns the TTree.

#include "G4UserRunAction.hh"
#include <string>

class TFile;
class TTree;
class TrackerEventAction;

class TrackerRunAction : public G4UserRunAction {
public:
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
