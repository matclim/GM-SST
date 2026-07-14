// reco: ShipStrawGeometry.cpp   (ACTS main @ b660c71)
#include "ShipStrawGeometry.hpp"

#include <stdexcept>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Material/HomogeneousSurfaceMaterial.hpp>
#include <Acts/Material/Material.hpp>
#include <Acts/Material/MaterialSlab.hpp>
#include <Acts/Surfaces/LineBounds.hpp>
#include <Acts/Surfaces/StrawSurface.hpp>

#include <TFile.h>
#include <TTree.h>

namespace shipreco {

namespace {

// ---------------------------------------------------------------------------
// STRAW MATERIAL
//
// The sim's straws are a 30 um Mylar wall around ~20 mm of Ar/CO2. A track
// crossing one traverses TWO walls plus the gas:
//
//     Mylar : 2 x 30 um / X0(287 mm)      = 2.09e-4
//     gas   : ~19.94 mm / X0(~130 m)      = 1.53e-4
//     TOTAL                                 3.62e-4  per straw
//
// ACTS attaches ONE homogeneous slab per surface, so we fold both into a single
// effective Mylar slab whose radiation-length budget matches exactly:
//
//     t_eff = X0(Mylar) x 3.62e-4 = 0.104 mm
//
// That reproduces the multiple scattering (which is what drives the angular
// resolution, and hence the vertex resolution) and approximates the energy loss
// -- the gas carries little mass, so a Mylar-dominated dE/dx is close enough.
//
// Scale: over 32 straws this is x/X0 ~ 1.2%, giving theta0 ~ 243 urad at 5 GeV
// but only ~20 urad at 60 GeV. So MS matters a great deal for the K_S pions and
// rather little for the stiff LLP daughters.
//
// Caveat: the slab is applied regardless of where the track crosses the straw,
// whereas the real path through the wall depends on the impact parameter. This
// is the standard homogeneous approximation; it is an AVERAGE, not exact.
// ---------------------------------------------------------------------------
std::shared_ptr<const Acts::ISurfaceMaterial> makeStrawMaterial() {
  namespace U = Acts::UnitConstants;

  // Mylar (C10H8O4): X0 = 287 mm, lambda_I = 569 mm, <A> = 12.88, <Z> = 6.46,
  // rho = 1.4 g/cm3.
  //
  // ACTS native units are mm and g/mm^3, so the density is
  //     1.4 g/cm3 = 1.4e-3 g/mm3.
  // Written out with explicit constants rather than UDLs: getting a density
  // unit wrong is a silent factor-1000 that no compiler will catch.
  const float rho = static_cast<float>(1.4e-3 * U::g / (U::mm * U::mm * U::mm));

  const Acts::Material mylar = Acts::Material::fromMassDensity(
      static_cast<float>(287.0 * U::mm),   // X0
      static_cast<float>(569.0 * U::mm),   // L0 (nuclear interaction length)
      12.88f,                              // Ar
      6.46f,                               // Z
      rho);

  const float kEffThickness = static_cast<float>(0.104 * U::mm);  // budget above
  const Acts::MaterialSlab slab(mylar, kEffThickness);
  return std::make_shared<const Acts::HomogeneousSurfaceMaterial>(slab);
}

}  // namespace

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

  // One shared material object for all 9600 straws.
  auto strawMaterial = makeStrawMaterial();

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
    // Without this the Kalman filter treats the tracker as vacuum: it does not
    // inflate the covariance between measurements, so the fit is over-confident
    // and mis-weights the hits. KalmanFitterOptions already defaults
    // multipleScattering = energyLoss = true -- it just had no material to act on.
    surf->assignSurfaceMaterial(strawMaterial);
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
