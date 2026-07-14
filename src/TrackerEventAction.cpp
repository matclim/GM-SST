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
    tree->Branch("parentID",   &m_parentID);
    tree->Branch("pdg",        &m_pdg);
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
    tree->Branch("vtxX",       &m_vx);
    tree->Branch("vtxY",       &m_vy);
    tree->Branch("vtxZ",       &m_vz);
    tree->Branch("vpx",        &m_px);
    tree->Branch("vpy",        &m_py);
    tree->Branch("vpz",        &m_pz);
    tree->Branch("driftTime",  &m_driftTime);   // ns -- THE measurement
    tree->Branch("driftTrue",  &m_driftTrue);   // mm -- truth, diagnostic
    tree->Branch("driftTime",  &m_driftTime);   // ns, THE measurement
    tree->Branch("driftTrue",  &m_driftTrue);   // mm, truth (diagnostic)
}

void TrackerEventAction::BeginOfEventAction(const G4Event*) {
    m_trackID.clear();   m_parentID.clear();   m_pdg.clear();   m_stationID.clear();
    m_layerID.clear();   m_subLayerID.clear();  m_strawID.clear();
    m_edep.clear();
    m_x.clear();  m_y.clear();  m_z.clear();
    m_xe.clear(); m_ye.clear(); m_ze.clear();
    m_xx.clear(); m_yx.clear(); m_zx.clear();
    m_vx.clear(); m_vy.clear(); m_vz.clear();
    m_px.clear(); m_py.clear(); m_pz.clear();
    m_driftTime.clear(); m_driftTrue.clear();
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
        m_parentID  .push_back(h.parentID);
        m_pdg       .push_back(h.pdg);
        m_stationID .push_back(h.stationID);
        m_layerID   .push_back(h.layerID);
        m_subLayerID.push_back(h.subLayerID);
        m_strawID   .push_back(h.strawID);
        m_edep      .push_back(h.edep);
        m_x .push_back(h.x);   m_y .push_back(h.y);   m_z .push_back(h.z);
        m_xe.push_back(h.x_entry); m_ye.push_back(h.y_entry); m_ze.push_back(h.z_entry);
        m_xx.push_back(h.x_exit);  m_yx.push_back(h.y_exit);  m_zx.push_back(h.z_exit);
        m_vx.push_back(h.vtxX);    m_vy.push_back(h.vtxY);    m_vz.push_back(h.vtxZ);
        m_px.push_back(h.vpx);     m_py.push_back(h.vpy);     m_pz.push_back(h.vpz);
        m_driftTime.push_back(h.driftTime);
        m_driftTrue.push_back(h.driftTrue);
    }

    m_tree->Fill();
}
