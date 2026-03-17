// TrackerRunAction.cpp
// Opens/closes the ROOT file and owns the TTree.

#include "TrackerRunAction.h"
#include "TrackerEventAction.h"

#include "TFile.h"
#include "TTree.h"

#include "G4Run.hh"
#include <iostream>

TrackerRunAction::TrackerRunAction(TrackerEventAction* evtAction,
                                   const std::string& outFile)
    : G4UserRunAction()
    , m_evtAction(evtAction)
    , m_outFile(outFile)
{}

void TrackerRunAction::BeginOfRunAction(const G4Run*) {
    m_rootFile = TFile::Open(m_outFile.c_str(), "RECREATE");
    if (!m_rootFile || m_rootFile->IsZombie()) {
        std::cerr << "[TrackerRunAction] ERROR: cannot open " << m_outFile << "\n";
        return;
    }

    m_tree = new TTree("Events", "Straw tracker hit tree");
    m_evtAction->setTree(m_tree);

    std::cout << "[TrackerRunAction] Output file: " << m_outFile << "\n";
}

void TrackerRunAction::EndOfRunAction(const G4Run* run) {
    if (m_rootFile && !m_rootFile->IsZombie()) {
        m_rootFile->Write();
        m_rootFile->Close();
        std::cout << "[TrackerRunAction] Wrote " << run->GetNumberOfEvent()
                  << " events to " << m_outFile << "\n";
    }
}
