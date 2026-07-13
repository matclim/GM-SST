#pragma once
// TrackerPrimaryGeneratorAction.h
//
// Two modes:
//   * particle gun  (default) — one particle, set via the CLI
//   * LLP file      (--llp-file) — read an LLP decay sample and fire ALL of its
//     charged daughters from the recorded decay vertex. The LLP itself is never
//     tracked: it is neutral and has already decayed. Geant4 therefore sees N
//     primaries emitted from one point, which is exactly the topology we want
//     to vertex.
//
// The LLP file (e.g. data/DP_4pi.root) stores, per event:
//     LLP_px/py/pz/E/m [GeV], LLP_pdg, LLP_weight
//     vtx_x/y/z        [m], in the SHiP frame
//     d_px/py/pz/E/m   [GeV] (vector, one per daughter), d_pdg (vector)
// Units are converted on the way in: GeV -> MeV, m -> mm.

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"

#include <memory>
#include <string>
#include <vector>

class G4Event;
class TFile;
class TTree;

class TrackerPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    TrackerPrimaryGeneratorAction();
    ~TrackerPrimaryGeneratorAction() override;

    void GeneratePrimaries(G4Event*) override;

    G4ParticleGun* gun() { return m_gun.get(); }

    /// Switch to LLP mode. worldZOriginMM converts SHiP mm -> Geant4 local mm.
    /// Returns false (and stays in gun mode) if the file cannot be read.
    bool openLLPFile(const std::string& path, double worldZOriginMM);

    /// Number of events available in the LLP file (0 in gun mode).
    long llpEntries() const { return m_nEntries; }

private:
    void generateFromLLP(G4Event*);

    std::unique_ptr<G4ParticleGun> m_gun;

    // ---- LLP mode ----
    bool        m_llpMode      {false};
    double      m_worldZOrigin {60000.0};   // mm
    TFile*      m_file         {nullptr};
    TTree*      m_tree         {nullptr};
    long        m_nEntries     {0};
    long        m_next         {0};         // next entry to read

    // branch buffers
    float                m_vx{0}, m_vy{0}, m_vz{0};       // m
    std::vector<float>  *m_dpx{nullptr}, *m_dpy{nullptr}, *m_dpz{nullptr};
    std::vector<float>  *m_dE {nullptr};
    std::vector<int>    *m_dpdg{nullptr};
    float                m_weight{1.0f};
};
