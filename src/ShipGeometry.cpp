// src/ShipGeometry.cpp  (real straw geometry; ACTS main API)
//
// Each straw layer is a PlaneSurface normal to z. Straws run along X, so a
// layer measures the coordinate PERPENDICULAR to the wire = Y (the bend
// coordinate). With a stereo rotation alpha of the wires about z, the measured
// (loc0) direction is (-sin a, cos a, 0), so the surface is rotated by a+90 deg
// and loc0 = -x*sin a + y*cos a. loc1 runs along the wire (no measurement).
#include "ShipGeometry.hpp"

#include <algorithm>
#include <cmath>

#include <Acts/Definitions/Units.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/RectangleBounds.hpp>

namespace shipreco {

Acts::GeometryIdentifier ShipGeometry::geoId(int station, int layer) {
  return Acts::GeometryIdentifier()
      .withVolume(static_cast<Acts::GeometryIdentifier::Value>(station + 1))
      .withLayer(static_cast<Acts::GeometryIdentifier::Value>(layer + 1));
}

ShipGeometry::ShipGeometry(const std::vector<LayerDef>& layersIn) {
  std::vector<LayerDef> layers = layersIn;
  std::sort(layers.begin(), layers.end(),
            [](const LayerDef& a, const LayerDef& b) { return a.zWorld < b.zWorld; });

  for (const auto& L : layers) {
    Acts::Transform3 tf = Acts::Transform3::Identity();
    tf.translation() = Acts::Vector3(0.0, 0.0, L.zWorld);
    // Rotate by (stereo + 90 deg): surface local-x (loc0) becomes the measured
    // direction perpendicular to the wire.
    tf.rotate(Acts::AngleAxis3(L.stereoRad + M_PI / 2.0, Acts::Vector3::UnitZ()));

    // Bounds: local-x (loc0 = measured, ~Y) half = 3000; local-y (along wire,
    // ~X) half = 2000.  (kStationY/2 = 3000, kStationX/2 = 2000.)
    auto bounds = std::make_shared<Acts::RectangleBounds>(L.halfX, L.halfY);
    auto surf = Acts::Surface::makeShared<Acts::PlaneSurface>(tf, bounds);
    surf->assignGeometryId(geoId(L.station, L.layer));

    m_byId.emplace(geoId(L.station, L.layer), surf.get());
    m_ordered.push_back(surf.get());
    m_surfaces.push_back(std::move(surf));
  }
}

const Acts::Surface* ShipGeometry::surfaceForGeoId(Acts::GeometryIdentifier id) const {
  auto it = m_byId.find(id);
  return (it == m_byId.end()) ? nullptr : it->second;
}

std::vector<LayerDef> ShipGeometry::defaultLayers(double worldZOriginMm) {
  // Real StrawTracker geometry (from StrawTrackerBuilder):
  //   4 stations at these lab-z centres; 4 layers each at within-station
  //   offsets of -88.5,-29.5,+29.5,+88.5 mm (pitch 59); stereo +/-2.3 deg.
  const double stZ[4] = {26500., 29000., 34000., 35500.};
  const double layHalfZ   = 27.0;               // kFrameHalfZ + clearance
  const double layerGap   = 5.0;
  const double layerPitch = 2.0 * layHalfZ + layerGap;   // 59 mm
  const double stereoDeg  = 2.3;
  const double halfMeasY  = 3000.0;   // loc0 half-range (Y aperture / 2)
  const double halfWireX  = 2000.0;   // loc1 half-range (X aperture / 2)

  std::vector<LayerDef> out;
  for (int s = 0; s < 4; ++s) {
    for (int l = 0; l < 4; ++l) {
      const double zLay  = -0.5 * (4 - 1) * layerPitch + l * layerPitch;
      const double labZ  = stZ[s] + zLay;
      const double alpha = ((l % 2 == 0) ? +1.0 : -1.0) * stereoDeg * M_PI / 180.0;
      out.push_back(LayerDef{s, l, labZ - worldZOriginMm, alpha, halfMeasY, halfWireX});
    }
  }
  return out;
}

}  // namespace shipreco
