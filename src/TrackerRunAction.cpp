// TrackerRunAction.cpp
// Each worker thread opens its own ROOT file named
//   <base>_t<threadID>.root
// so there are no concurrent write conflicts. The master does nothing
// with ROOT. After the run, files can be merged with:
//   hadd StrawTracker_hits.root StrawTracker_hits_t*.root

#include "TrackerRunAction.h"
#include "TrackerEventAction.h"

#include "TFile.h"
#include "TTree.h"

#include "G4Run.hh"
#include "G4Threading.hh"
#include <iostream>
#include <sstream>

TrackerRunAction::TrackerRunAction(TrackerEventAction* evtAction,
                                   const std::string& outFile)
    : G4UserRunAction()
    , m_evtAction(evtAction)
    , m_outFile(outFile)
{}

void TrackerRunAction::BeginOfRunAction(const G4Run*) {
    // Master thread has no event action — nothing to do for ROOT.
    if (!m_evtAction) return;

    // Build a per-thread filename: insert _t<id> before the extension.
    std::string fname = m_outFile;
    const std::string ext = ".root";
    std::ostringstream suffix;
    suffix << "_t" << G4Threading::G4GetThreadId();

    auto pos = fname.rfind(ext);
    if (pos != std::string::npos)
        fname.insert(pos, suffix.str());
    else
        fname += suffix.str() + ext;

    m_rootFile = TFile::Open(fname.c_str(), "RECREATE");
    if (!m_rootFile || m_rootFile->IsZombie()) {
        std::cerr << "[TrackerRunAction] ERROR: cannot open " << fname << "\n";
        return;
    }

    m_tree = new TTree("Events", "Straw tracker hit tree");
    m_evtAction->setTree(m_tree);

    std::cout << "[TrackerRunAction] Thread " << G4Threading::G4GetThreadId()
              << " writing to " << fname << "\n";
}

void TrackerRunAction::EndOfRunAction(const G4Run* run) {
    if (!m_evtAction || !m_rootFile || m_rootFile->IsZombie()) return;

    m_rootFile->Write();
    m_rootFile->Close();
    std::cout << "[TrackerRunAction] Thread " << G4Threading::G4GetThreadId()
              << " wrote " << run->GetNumberOfEvent() << " events.\n";
}
