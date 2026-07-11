#pragma once
// reco: ShipHitReader.hpp — read the Geant4 Events tree into memory.
#include <string>
#include <vector>

namespace shipreco {

struct RawHit {
  int    trackID{0};
  int    parentID{0};
  int    pdg{0};
  int    stationID{0};
  int    layerID{0};
  int    subLayerID{0};
  int    strawID{0};
  double x{0}, y{0}, z{0};        // step mid-point
  double xe{0}, ye{0}, ze{0};    // entry
  double xx{0}, yx{0}, zx{0};    // exit
  double vtxX{0}, vtxY{0}, vtxZ{0};  // truth production (decay) vertex
  double vpx{0}, vpy{0}, vpz{0};     // truth momentum at production (MeV)
};

struct RawEvent {
  std::vector<RawHit> hits;
  bool   hasTruthVtx{false};
  double truthVtxX{0}, truthVtxY{0}, truthVtxZ{0};  // event decay vertex (mean of secondaries)
};

/// primaryOnly=false keeps ALL tracks (needed for multi-track / vertexing).
std::vector<RawEvent> readEvents(const std::string& rootFile,
                                 bool primaryOnly = false);

}  // namespace shipreco
