#pragma once
// reco: ShipStrawGeometry.hpp
// Builds ACTS StrawSurfaces (line surfaces along each wire) from a straw table
// dumped by the sim (which sources it from the GeoModel geometry). The reco has
// NO GeoModel dependency: it only reads (ids + wire center + wire direction).
//
// Straw table = ROOT TTree "Straws" with per-straw branches:
//   station, layer, subLayer, straw   (int)
//   cx, cy, cz                        (double, wire CENTRE, global mm)
//   ux, uy, uz                        (double, wire DIRECTION, unit)
//   [optional] radius, halfLength     (double, mm; default 10 / 2000)

#include <map>
#include <memory>
#include <string>
#include <vector>

#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Surfaces/Surface.hpp>

namespace shipreco {

class ShipStrawGeometry {
 public:
  /// Load the straw table and build a StrawSurface per wire.
  void loadTable(const std::string& rootFile);

  /// Stable GeometryIdentifier for a straw (packs the four IDs).
  static Acts::GeometryIdentifier geoId(int station, int layer,
                                        int subLayer, int straw);

  /// Surface lookups.
  const Acts::Surface* surfaceFor(int station, int layer,
                                  int subLayer, int straw) const;
  const Acts::Surface* surfaceForGeoId(Acts::GeometryIdentifier id) const;

  std::size_t size() const { return m_surfaces.size(); }

 private:
  std::vector<std::shared_ptr<Acts::Surface>>              m_surfaces;
  std::map<Acts::GeometryIdentifier, const Acts::Surface*> m_byId;
};

}  // namespace shipreco
