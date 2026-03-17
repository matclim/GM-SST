#pragma once
// TrackerHit.h
// Plain data structs for straw hits collected per event.

#include <vector>

struct StrawHit {
    int    trackID   {-1};
    int    stationID {-1};
    int    layerID   {-1};
    int    subLayerID{-1};
    int    strawID   {-1};
    double edep      {0.};   // MeV
    double x         {0.};   // mm  (step centroid)
    double y         {0.};
    double z         {0.};
    double x_entry   {0.};
    double y_entry   {0.};
    double z_entry   {0.};
    double x_exit    {0.};
    double y_exit    {0.};
    double z_exit    {0.};
};
