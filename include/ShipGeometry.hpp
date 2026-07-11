#pragma once
// reco/geometry/ShipGeometry.hpp
// Builds the tracking surfaces for the straw spectrometer.
//
// Design choice: each straw LAYER is modelled as one Acts::PlaneSurface at its
// z, with an in-plane rotation `stereoRad` that orients the measured direction
// (loc0). This measures one coordinate per layer — the natural straw model —
// and needs NO TrackingGeometry: the fit uses Acts::DirectNavigator with the
// ordered surface list returned by orderedSurfaces().
//
// Fill `defaultLayers()` (or build your own vector) from StrawTrackerBuilder
// so the z / stereo / bounds match your real geometry. Coordinates are in the
// ACTS global frame; if that frame == your Geant4 world (centre lab z=31000),
// pass z_world = z_lab - 31000.

#include <cstdint>
#include <map>
#include <memory>
#include <vector>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Surfaces/Surface.hpp>

namespace shipreco {

/// One measurement plane (a straw layer).
struct LayerDef {
  int    station{0};
  int    layer{0};
  double zWorld{0.0};      // mm, ACTS global frame
  double stereoRad{0.0};   // in-plane rotation of the measurement axis (loc0)
  double halfX{2500.0};    // plane half-size (mm)
  double halfY{3500.0};
};

class ShipGeometry {
 public:
  explicit ShipGeometry(const std::vector<LayerDef>& layers);

  /// Surfaces in beam (z) order — hand this to the KalmanFitter/DirectNavigator.
  const std::vector<const Acts::Surface*>& orderedSurfaces() const { return m_ordered; }

  /// Look up a surface by its GeometryIdentifier (used by the SourceLink accessor).
  const Acts::Surface* surfaceForGeoId(Acts::GeometryIdentifier id) const;

  /// Stable GeometryIdentifier for a (station,layer) pair.
  static Acts::GeometryIdentifier geoId(int station, int layer);

  /// A plausible 4-station layout — REPLACE the numbers with your real ones.
  static std::vector<LayerDef> defaultLayers(double worldZOriginMm = 31000.0);

 private:
  std::vector<std::shared_ptr<Acts::Surface>>              m_surfaces;
  std::vector<const Acts::Surface*>                        m_ordered;
  std::map<Acts::GeometryIdentifier, const Acts::Surface*> m_byId;
};

}  // namespace shipreco
