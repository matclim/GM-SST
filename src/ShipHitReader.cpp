// reco: ShipHitReader.cpp
#include "ShipHitReader.hpp"
#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace shipreco {
std::vector<RawEvent> readEvents(const std::string& rootFile, bool primaryOnly) {
  std::unique_ptr<TFile> f(TFile::Open(rootFile.c_str(), "READ"));
  if (!f || f->IsZombie())
    throw std::runtime_error("ShipHitReader: cannot open '" + rootFile + "'");
  auto* T = dynamic_cast<TTree*>(f->Get("Events"));
  if (!T) throw std::runtime_error("ShipHitReader: no 'Events' tree");

  std::vector<int>*    trackID=nullptr,*parentID=nullptr,*pdg=nullptr,*stationID=nullptr;
  std::vector<int>*    layerID=nullptr,*subLayerID=nullptr,*strawID=nullptr;
  std::vector<double>* x=nullptr,*y=nullptr,*z=nullptr;
  std::vector<double>* xe=nullptr,*ye=nullptr,*ze=nullptr,*xx=nullptr,*yx=nullptr,*zx=nullptr;
  std::vector<double>* vx=nullptr,*vy=nullptr,*vz=nullptr;
  std::vector<double>* px=nullptr,*py=nullptr,*pz=nullptr;
  T->SetBranchAddress("trackID",&trackID);   T->SetBranchAddress("stationID",&stationID);
  T->SetBranchAddress("layerID",&layerID);   T->SetBranchAddress("subLayerID",&subLayerID);
  T->SetBranchAddress("strawID",&strawID);
  T->SetBranchAddress("x",&x); T->SetBranchAddress("y",&y); T->SetBranchAddress("z",&z);
  T->SetBranchAddress("x_entry",&xe); T->SetBranchAddress("y_entry",&ye); T->SetBranchAddress("z_entry",&ze);
  T->SetBranchAddress("x_exit",&xx);  T->SetBranchAddress("y_exit",&yx);  T->SetBranchAddress("z_exit",&zx);
  const bool hasParent = T->GetBranch("parentID"); if (hasParent) T->SetBranchAddress("parentID",&parentID);
  const bool hasPdg    = T->GetBranch("pdg");       if (hasPdg)    T->SetBranchAddress("pdg",&pdg);
  const bool hasVtx    = T->GetBranch("vtxX");
  if (hasVtx) { T->SetBranchAddress("vtxX",&vx); T->SetBranchAddress("vtxY",&vy); T->SetBranchAddress("vtxZ",&vz); }
  const bool hasMom   = T->GetBranch("vpx");
  const bool hasDrift  = T->GetBranch("driftTime");
  const bool hasWeight = T->GetBranch("weight");
  if (hasMom) { T->SetBranchAddress("vpx",&px); T->SetBranchAddress("vpy",&py); T->SetBranchAddress("vpz",&pz); }
  std::vector<double> *dt=nullptr, *dtr=nullptr;
  if (hasDrift) { T->SetBranchAddress("driftTime",&dt); T->SetBranchAddress("driftTrue",&dtr); }

  std::vector<RawEvent> events;
  const Long64_t n = T->GetEntries(); events.reserve(n);
  for (Long64_t e = 0; e < n; ++e) {
    T->GetEntry(e);
    RawEvent ev;
    double sx=0, sy=0, sz=0; int nsec=0;
    for (std::size_t h = 0; h < stationID->size(); ++h) {
      const int par = hasParent ? (*parentID)[h] : 0;
      if (primaryOnly && (*trackID)[h] != 1) continue;
      RawHit rh;
      rh.trackID=(*trackID)[h]; rh.parentID=par; rh.pdg=hasPdg?(*pdg)[h]:0;
      rh.stationID=(*stationID)[h]; rh.layerID=(*layerID)[h];
      rh.subLayerID=(*subLayerID)[h]; rh.strawID=(*strawID)[h];
      rh.x=(*x)[h]; rh.y=(*y)[h]; rh.z=(*z)[h];
      rh.xe=(*xe)[h]; rh.ye=(*ye)[h]; rh.ze=(*ze)[h];
      rh.xx=(*xx)[h]; rh.yx=(*yx)[h]; rh.zx=(*zx)[h];
      if (hasMom)   { rh.vpx=(*px)[h]; rh.vpy=(*py)[h]; rh.vpz=(*pz)[h]; }
      if (hasDrift) { rh.driftTime=(*dt)[h]; rh.driftTrue=(*dtr)[h]; }
      if (hasVtx) { rh.vtxX=(*vx)[h]; rh.vtxY=(*vy)[h]; rh.vtxZ=(*vz)[h];
                    if (par == 1) { sx+=rh.vtxX; sy+=rh.vtxY; sz+=rh.vtxZ; ++nsec; } }
      ev.hits.push_back(rh);
    }
    if (hasVtx && nsec > 0) { ev.hasTruthVtx=true;
      ev.truthVtxX=sx/nsec; ev.truthVtxY=sy/nsec; ev.truthVtxZ=sz/nsec; }
    events.push_back(std::move(ev));
  }
  return events;
}
}  // namespace shipreco
