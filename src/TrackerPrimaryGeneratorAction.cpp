// TrackerPrimaryGeneratorAction.cpp
//
// Gun mode: one particle, configured from the CLI.
// LLP mode: replay a decay sample -- fire every CHARGED daughter of the LLP from
//           its recorded decay vertex. See the header for the file format.

#include "TrackerPrimaryGeneratorAction.h"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryVertex.hh"
#include "G4PrimaryParticle.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4Event.hh"

#include <TFile.h>
#include <TTree.h>

#include <cmath>
#include <iostream>

namespace {
// The LLP file is in GeV and metres; Geant4 here is MeV and mm.
constexpr double kGeVToMeV = 1000.0;
constexpr double kMToMM    = 1000.0;
}  // namespace

TrackerPrimaryGeneratorAction::TrackerPrimaryGeneratorAction()
    : G4VUserPrimaryGeneratorAction()
    , m_gun(std::make_unique<G4ParticleGun>(1))
{
    auto* mu = G4ParticleTable::GetParticleTable()->FindParticle("mu-");
    m_gun->SetParticleDefinition(mu);
    m_gun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
    m_gun->SetParticleEnergy(10000. * MeV);
    m_gun->SetParticlePosition(G4ThreeVector(0., 0., 19320. * mm));  // SHiP 79320
}

TrackerPrimaryGeneratorAction::~TrackerPrimaryGeneratorAction() {
    if (m_file) { m_file->Close(); delete m_file; }
}

bool TrackerPrimaryGeneratorAction::openLLPFile(const std::string& path,
                                                double worldZOriginMM) {
    m_file = TFile::Open(path.c_str(), "READ");
    if (!m_file || m_file->IsZombie()) {
        std::cerr << "[PrimaryGen] ERROR: cannot open LLP file '" << path << "'\n";
        m_file = nullptr;
        return false;
    }
    m_tree = dynamic_cast<TTree*>(m_file->Get("Events"));
    if (!m_tree) {
        std::cerr << "[PrimaryGen] ERROR: no 'Events' tree in '" << path << "'\n";
        return false;
    }

    m_tree->SetBranchAddress("vtx_x", &m_vx);
    m_tree->SetBranchAddress("vtx_y", &m_vy);
    m_tree->SetBranchAddress("vtx_z", &m_vz);
    m_tree->SetBranchAddress("d_px",  &m_dpx);
    m_tree->SetBranchAddress("d_py",  &m_dpy);
    m_tree->SetBranchAddress("d_pz",  &m_dpz);
    m_tree->SetBranchAddress("d_E",   &m_dE);
    m_tree->SetBranchAddress("d_pdg", &m_dpdg);
    if (m_tree->GetBranch("LLP_weight"))
        m_tree->SetBranchAddress("LLP_weight", &m_weight);

    m_nEntries     = m_tree->GetEntries();
    m_worldZOrigin = worldZOriginMM;
    m_llpMode      = true;
    m_next         = 0;

    std::cout << "[PrimaryGen] LLP mode: '" << path << "', "
              << m_nEntries << " decays\n"
              << "[PrimaryGen] units: GeV -> MeV, m -> mm;  SHiP z -> local z - "
              << m_worldZOrigin << " mm\n";
    return true;
}

void TrackerPrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
    if (m_llpMode) generateFromLLP(event);
    else           m_gun->GeneratePrimaryVertex(event);
}

void TrackerPrimaryGeneratorAction::generateFromLLP(G4Event* event) {
    if (!m_tree || m_nEntries == 0) return;

    // Replay sequentially; wrap around if more events are requested than exist.
    const long entry = m_next % m_nEntries;
    ++m_next;
    m_tree->GetEntry(entry);

    // Decay vertex: SHiP metres -> Geant4 local mm.
    const double vxG4 =  m_vx * kMToMM * mm;
    const double vyG4 =  m_vy * kMToMM * mm;
    const double vzG4 = (m_vz * kMToMM - m_worldZOrigin) * mm;

    auto* vertex = new G4PrimaryVertex(G4ThreeVector(vxG4, vyG4, vzG4), 0.0);

    auto* table = G4ParticleTable::GetParticleTable();
    const std::size_t nd = m_dpdg ? m_dpdg->size() : 0;
    int nCharged = 0;
    m_nBodiesTrue = static_cast<int>(nd);   // total decay bodies (charged + neutral)

    for (std::size_t i = 0; i < nd; ++i) {
        const int pdg = (*m_dpdg)[i];
        auto* def = table->FindParticle(pdg);
        if (!def) {
            std::cerr << "[PrimaryGen] unknown PDG " << pdg << ", skipping\n";
            continue;
        }
        // Only CHARGED daughters leave straw hits; neutrals would just be noise
        // in the tracker, and we cannot vertex them anyway.
        if (std::abs(def->GetPDGCharge()) < 1e-6) continue;

        auto* p = new G4PrimaryParticle(def);
        p->SetMomentum((*m_dpx)[i] * kGeVToMeV * MeV,
                       (*m_dpy)[i] * kGeVToMeV * MeV,
                       (*m_dpz)[i] * kGeVToMeV * MeV);
        p->SetWeight(m_weight);
        vertex->SetPrimary(p);
        ++nCharged;
    }

    m_nBodiesChargedTrue = nCharged;   // charged bodies actually fired

    if (nCharged == 0) {
        // Nothing trackable; still add the (empty) vertex so event IDs line up.
        std::cerr << "[PrimaryGen] entry " << entry << ": no charged daughters\n";
    }
    event->AddPrimaryVertex(vertex);
}
