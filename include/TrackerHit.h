#pragma once
// TrackerHit.h — plain data structs for straw hits collected per event.
#include <vector>

struct StrawHit {
    int    trackID   {-1};
    int    parentID  {-1};   // 0 = primary; 1 = direct daughter of primary
    int    pdg       {0};    // PDG code of the track
    int    stationID {-1};
    int    layerID   {-1};
    int    subLayerID{-1};
    int    strawID   {-1};
    double edep      {0.};   // MeV
    double x         {0.};   // mm (step centroid)
    double y         {0.};
    double z         {0.};
    double x_entry   {0.};
    double y_entry   {0.};
    double z_entry   {0.};
    double x_exit    {0.};
    double y_exit    {0.};
    double z_exit    {0.};
    double vtxX      {0.};   // production vertex (= decay vertex for daughters)
    double vtxY      {0.};
    double vtxZ      {0.};
    double vpx       {0.};   // TRUTH momentum at production (MeV)
    double vpy       {0.};
    double vpz       {0.};
};
