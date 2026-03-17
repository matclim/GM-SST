// TrackerEventAction.cpp
#include "TrackerEventAction.h"
#include "TrackerSD.h"

#include "G4SDManager.hh"
#include "G4Event.hh"
#include "TTree.h"

void TrackerEventAction::setTree(TTree* tree) {
    m_tree = tree;
    if (!tree) return;

    tree->Branch("trackID",    &m_trackID);
    tree->Branch("stationID",  &m_stationID);
    tree->Branch("layerID",    &m_layerID);
    tree->Branch("subLayerID", &m_subLayerID);
    tree->Branch("strawID",    &m_strawID);
    tree->Branch("edep",       &m_edep);
    tree->Branch("x",          &m_x);
    tree->Branch("y",          &m_y);
    tree->Branch("z",          &m_z);
    tree->Branch("x_entry",    &m_xe);
    tree->Branch("y_entry",    &m_ye);
    tree->Branch("z_entry",    &m_ze);
    tree->Branch("x_exit",     &m_xx);
    tree->Branch("y_exit",     &m_yx);
    tree->Branch("z_exit",     &m_zx);
}

void TrackerEventAction::BeginOfEventAction(const G4Event*) {
    m_trackID.clear();   m_stationID.clear();
    m_layerID.clear();   m_subLayerID.clear();  m_strawID.clear();
    m_edep.clear();
    m_x.clear();  m_y.clear();  m_z.clear();
    m_xe.clear(); m_ye.clear(); m_ze.clear();
    m_xx.clear(); m_yx.clear(); m_zx.clear();
}

void TrackerEventAction::EndOfEventAction(const G4Event*) {
    if (!m_tree) return;

    // Lazy SD lookup — safe in MT because each worker has its own instance.
    if (!m_sd) {
        m_sd = static_cast<TrackerSD*>(
            G4SDManager::GetSDMpointer()->FindSensitiveDetector("StrawTrackerSD"));
    }
    if (!m_sd) return;

    for (const auto& h : m_sd->hits()) {
        m_trackID   .push_back(h.trackID);
        m_stationID .push_back(h.stationID);
        m_layerID   .push_back(h.layerID);
        m_subLayerID.push_back(h.subLayerID);
        m_strawID   .push_back(h.strawID);
        m_edep      .push_back(h.edep);
        m_x .push_back(h.x);   m_y .push_back(h.y);   m_z .push_back(h.z);
        m_xe.push_back(h.x_entry); m_ye.push_back(h.y_entry); m_ze.push_back(h.z_entry);
        m_xx.push_back(h.x_exit);  m_yx.push_back(h.y_exit);  m_zx.push_back(h.z_exit);
    }

    m_tree->Fill();
}
