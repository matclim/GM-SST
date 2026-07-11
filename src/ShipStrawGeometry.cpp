// reco: ShipStrawGeometry.cpp   (ACTS main @ b660c71)
#include "ShipStrawGeometry.hpp"

#include <stdexcept>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Surfaces/LineBounds.hpp>
#include <Acts/Surfaces/StrawSurface.hpp>

#include <TFile.h>
#include <TTree.h>

namespace shipreco {

Acts::GeometryIdentifier ShipStrawGeometry::geoId(int station, int layer,
                                                  int subLayer, int straw) {
  using V = Acts::GeometryIdentifier::Value;
  // ACTS main: immutable withXxx() setters (confirmed: withVolume/withLayer).
  // If withApproach/withSensitive are absent on your build, pack straw+subLayer
  // into withSensitive alone instead (see note in the reply).
  return Acts::GeometryIdentifier()
      .withVolume(static_cast<V>(station + 1))
      .withLayer(static_cast<V>(layer + 1))
      .withApproach(static_cast<V>(subLayer + 1))
      .withSensitive(static_cast<V>(straw + 1));
}

void ShipStrawGeometry::loadTable(const std::string& rootFile) {
  std::unique_ptr<TFile> f(TFile::Open(rootFile.c_str(), "READ"));
  if (!f || f->IsZombie())
    throw std::runtime_error("ShipStrawGeometry: cannot open '" + rootFile + "'");
  auto* t = dynamic_cast<TTree*>(f->Get("Straws"));
  if (!t) throw std::runtime_error("ShipStrawGeometry: no 'Straws' tree");

  int station = 0, layer = 0, subLayer = 0, straw = 0;
  double cx = 0, cy = 0, cz = 0, ux = 0, uy = 0, uz = 1;
  double radius = 10.0, halfLength = 2000.0;
  t->SetBranchAddress("station",  &station);
  t->SetBranchAddress("layer",    &layer);
  t->SetBranchAddress("subLayer", &subLayer);
  t->SetBranchAddress("straw",    &straw);
  t->SetBranchAddress("cx", &cx); t->SetBranchAddress("cy", &cy); t->SetBranchAddress("cz", &cz);
  t->SetBranchAddress("ux", &ux); t->SetBranchAddress("uy", &uy); t->SetBranchAddress("uz", &uz);
  if (t->GetBranch("radius"))     t->SetBranchAddress("radius", &radius);
  if (t->GetBranch("halfLength")) t->SetBranchAddress("halfLength", &halfLength);

  const Long64_t n = t->GetEntries();
  m_surfaces.reserve(n);
  for (Long64_t i = 0; i < n; ++i) {
    t->GetEntry(i);

    // Build the wire transform: local z-axis -> wire direction, origin at centre.
    Acts::Vector3 c(cx, cy, cz);
    Acts::Vector3 d(ux, uy, uz);
    if (d.norm() == 0) continue;
    d.normalize();
    Eigen::Quaternion<double> q =
        Eigen::Quaternion<double>::FromTwoVectors(Acts::Vector3::UnitZ(), d);
    Acts::Transform3 tf = Acts::Transform3::Identity();
    tf.linear() = q.toRotationMatrix();
    tf.translation() = c;

    // Line/straw bounds: drift radius + half wire length.
    auto lb = std::make_shared<Acts::LineBounds>(radius, halfLength);
    auto surf = Acts::Surface::makeShared<Acts::StrawSurface>(tf, lb);

    const auto gid = geoId(station, layer, subLayer, straw);
    surf->assignGeometryId(gid);
    m_byId.emplace(gid, surf.get());
    m_surfaces.push_back(std::move(surf));
  }
}

const Acts::Surface* ShipStrawGeometry::surfaceFor(int station, int layer,
                                                   int subLayer, int straw) const {
  return surfaceForGeoId(geoId(station, layer, subLayer, straw));
}

const Acts::Surface* ShipStrawGeometry::surfaceForGeoId(
    Acts::GeometryIdentifier id) const {
  auto it = m_byId.find(id);
  return (it == m_byId.end()) ? nullptr : it->second;
}

}  // namespace shipreco
