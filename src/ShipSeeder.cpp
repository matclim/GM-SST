// reco: ShipSeeder.cpp   (ACTS main @ b660c71)
//
// Direction: least-squares line through the upstream (pre-magnet) hits.
// Momentum:  3-point sagitta with the OUTER chord — mean y at all 4 stations,
//            chord between stations 0 & 3 (absorbs any incoming tilt), sagitta
//            = deviation of the middle stations (1,2, bracketing the magnet).
//            p ≈ 0.3·kBL / θ,  θ = 4·|sagitta| / L_lever.  Self-contained.
#include "ShipSeeder.hpp"

#include <algorithm>
#include <cmath>

#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Utilities/VectorHelpers.hpp>

namespace shipreco {

namespace {
struct LineFit { double x0, sx, y0, sy; int n; bool ok; };

LineFit fitLine(const RawEvent& ev, int sA, int sB) {
  double Sz = 0, Szz = 0, Sx = 0, Szx = 0, Sy = 0, Szy = 0; int n = 0;
  for (const auto& h : ev.hits) {
    if (h.stationID != sA && h.stationID != sB) continue;
    Sz += h.z; Szz += h.z * h.z;
    Sx += h.x; Szx += h.z * h.x;
    Sy += h.y; Szy += h.z * h.y;
    ++n;
  }
  LineFit f{0, 0, 0, 0, n, false};
  if (n < 2) return f;
  const double denom = n * Szz - Sz * Sz;
  if (std::fabs(denom) < 1e-9) return f;
  f.sx = (n * Szx - Sz * Sx) / denom; f.x0 = (Sx - f.sx * Sz) / n;
  f.sy = (n * Szy - Sz * Sy) / denom; f.y0 = (Sy - f.sy * Sz) / n;
  f.ok = true;
  return f;
}

// Mean (y,z) of all hits in one station.
bool stationMean(const RawEvent& ev, int station, double& y, double& z) {
  double Sy = 0, Sz = 0; int n = 0;
  for (const auto& h : ev.hits)
    if (h.stationID == station) { Sy += h.y; Sz += h.z; ++n; }
  if (n == 0) return false;
  y = Sy / n; z = Sz / n; return true;
}
}  // namespace

std::optional<Acts::BoundTrackParameters>
seedFromUpstream(const Contexts& ctx, const RawEvent& ev,
                 const Acts::Surface& reference,
                 double pFallbackMeV, double charge,
                 int upstreamA, int upstreamB,
                 int downstreamA, int downstreamB, double kBL_Tm) {
  using namespace Acts::UnitConstants;

  const LineFit up = fitLine(ev, upstreamA, upstreamB);
  if (!up.ok) return std::nullopt;
  Acts::Vector3 dir = Acts::Vector3(up.sx, up.sy, 1.0).normalized();

  // ---- momentum from the outer-chord sagitta -------------------------------
  double pMeV = pFallbackMeV;
  double y0, z0, y1, z1, y2, z2, y3, z3;
  if (stationMean(ev, 0, y0, z0) && stationMean(ev, 1, y1, z1) &&
      stationMean(ev, 2, y2, z2) && stationMean(ev, 3, y3, z3)) {
    const double L = z3 - z0;                        // outer lever (mm)
    if (std::fabs(L) > 1e-3) {
      const double sc = (y3 - y0) / L;               // chord slope (tilt-absorbing)
      const double zMid = 0.5 * (z1 + z2);
      const double chordY = y0 + sc * (zMid - z0);
      const double midY = 0.5 * (y1 + y2);
      const double sag = midY - chordY;              // pure bend sagitta (mm)
      const double theta = 4.0 * std::fabs(sag) / std::fabs(L);
      if (theta > 1e-6) {
        double pGeV = 0.3 * kBL_Tm / theta;
        pGeV = std::clamp(pGeV, 0.5, 60.0);
        pMeV = pGeV * 1000.0;
      }
    }
  }

  // Charge: if caller passes 0, derive the sign from the bend direction.
  double q = charge;
  if (q == 0.0) {
    double y0c,z0c,y1c,z1c,y2c,z2c,y3c,z3c;
    if (stationMean(ev,0,y0c,z0c) && stationMean(ev,1,y1c,z1c) &&
        stationMean(ev,2,y2c,z2c) && stationMean(ev,3,y3c,z3c)) {
      const double Lc = z3c - z0c;
      const double scc = (y3c - y0c) / Lc;
      const double sgc = 0.5*(y1c+y2c) - (y0c + scc*(0.5*(z1c+z2c)-z0c));
      q = (sgc >= 0.0) ? +1.0 : -1.0;   // convention; flip kChargeSign in app if needed
    } else q = -1.0;
  }
  const double p = pMeV * MeV;
  const double qOverP = q / p;

  const double zRef = reference.center(ctx.gctx).z();
  Acts::Vector3 posRef(up.x0 + up.sx * zRef, up.y0 + up.sy * zRef, zRef);
  auto locRes = reference.globalToLocal(ctx.gctx, posRef, dir);
  if (!locRes.ok()) return std::nullopt;
  const Acts::Vector2 loc = *locRes;

  Acts::BoundVector par = Acts::BoundVector::Zero();
  par[Acts::eBoundLoc0]   = loc[0];
  par[Acts::eBoundLoc1]   = loc[1];
  par[Acts::eBoundPhi]    = Acts::VectorHelpers::phi(dir);
  par[Acts::eBoundTheta]  = Acts::VectorHelpers::theta(dir);
  par[Acts::eBoundQOverP] = qOverP;
  par[Acts::eBoundTime]   = 0.0;

  Acts::SquareMatrix<Acts::eBoundSize> cov =
      Acts::SquareMatrix<Acts::eBoundSize>::Identity();
  cov(Acts::eBoundLoc0,   Acts::eBoundLoc0)   = 100.0;
  cov(Acts::eBoundLoc1,   Acts::eBoundLoc1)   = 100.0;
  cov(Acts::eBoundPhi,    Acts::eBoundPhi)    = 1e-4;
  cov(Acts::eBoundTheta,  Acts::eBoundTheta)  = 1e-4;
  cov(Acts::eBoundQOverP, Acts::eBoundQOverP) = std::pow(qOverP * 0.5, 2);
  cov(Acts::eBoundTime,   Acts::eBoundTime)   = 1.0;

  return Acts::BoundTrackParameters(reference.getSharedPtr(), par, cov,
                                    Acts::ParticleHypothesis::muon());
}

}  // namespace shipreco
