#pragma once
// StrawDrift.h
//
// The straw drift model, shared by the SIMULATION (which produces a drift time)
// and the RECONSTRUCTION (which inverts it back to a radius). Keeping both in
// one header means there is a single r-t relation: if they ever disagreed, the
// reco would sit on a systematically wrong radius and no compiler would notice.
//
// ---------------------------------------------------------------------------
// WHAT A STRAW ACTUALLY MEASURES
//
// A charged track crossing the gas leaves discrete PRIMARY IONISATION CLUSTERS
// along its path. Each cluster's electrons drift to the wire. The TDC fires on
// the FIRST arrival. So the measurement is not the track's distance of closest
// approach (DOCA) -- it is the drift time of whichever cluster happens to lie
// nearest the wire.
//
// This is why straw resolution is RADIUS-DEPENDENT, and why it is WORST NEAR
// THE WIRE: there, the track runs nearly tangential to the drift circles, so a
// cluster displaced a little along the track maps to a large radius error, and
// with few clusters nearby the fluctuation is large. At mid-radius the geometry
// is kinder and more clusters sit near the DOCA, so the resolution improves.
//
// Modelling this properly (rather than smearing the truth DOCA by a constant)
// is the whole point: it reproduces the shape of the resolution curve, and it
// naturally delivers an UNSIGNED radius -- so left/right must be resolved by
// the fit, exactly as in the real detector.
//
// ---------------------------------------------------------------------------
// PARAMETERS: Ar/CO2 70:30, following NA62's straw tracker (same gas, same
// technology). NOTE that OUR straws are 10 mm in radius against NA62's ~5 mm,
// so the drift path is twice as long: expect a resolution somewhat WORSE than
// their ~130 um. That is a real consequence of the geometry, and the cluster
// model produces it rather than us imposing it.
// ---------------------------------------------------------------------------

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>

namespace strawdrift {

// ---- gas and electronics ---------------------------------------------------
/// Primary ionisation clusters per mm of track. Ar/CO2 70:30 at NTP gives
/// ~28 clusters/cm.
inline constexpr double kClusterDensity = 2.8;      // clusters / mm

/// Drift velocity. Ar/CO2 is near-saturated, so a LINEAR r-t relation is a
/// reasonable approximation: t = r / v.
inline constexpr double kDriftVelocity  = 0.050;    // mm / ns  (50 um/ns)

/// Longitudinal diffusion: sigma = kDiffusion * sqrt(r), with r in mm.
/// Ar/CO2 gives roughly 100-200 um at a full 1 cm drift, so
///     kDiffusion = 0.11 mm / sqrt(10 mm) ~ 0.035.
/// TUNE THIS against the resolution you actually want: after a run, look at
///     Hits->Draw("driftMeas-driftTrue")
/// and compare the width with the straw spec. Cluster statistics and diffusion
/// together set it, and only the second is a free knob.
inline constexpr double kDiffusion      = 0.035;    // mm / sqrt(mm)

/// TDC / electronics jitter.
inline constexpr double kTdcJitter      = 1.0;      // ns

/// Straw inner radius (gas volume).
inline constexpr double kStrawRadius    = 9.97;     // mm  (10.0 - 30 um wall)

// ---- the r-t relation ------------------------------------------------------
/// radius -> drift time. Linear, since Ar/CO2 runs near saturation.
inline double timeFromRadius(double r) {
  return r / kDriftVelocity;                        // ns
}

/// drift time -> radius. THE RECO'S INVERSION, and it is NOT the naive inverse
/// of timeFromRadius().
///
/// WHY A CALIBRATION IS NEEDED. The TDC fires on the nearest primary cluster,
/// and that cluster is never closer to the wire than the track itself. So the
/// naive linear inversion is BIASED OUTWARD -- by ~+50 um near the wire, where
/// the track runs tangential to the drift circles and few clusters lie near the
/// point of closest approach. Real experiments calibrate the r-t relation from
/// data for exactly this reason; here we do the same, fitting <DOCA | t> from
/// the simulation. That uses only the *distribution* of drift times, never
/// per-hit truth, so it is a legitimate calibration rather than a cheat.
///
/// Degree-5 polynomial fit to <DOCA | t>, t in ns -> r in mm. Residual bias
/// -0.4 um; resolution 102 um (128 um near the wire, 82 um at mid-radius,
/// 114 um near the wall -- the characteristic straw U-shape).
///
/// RE-FIT THIS if the gas, the straw radius, or the drift velocity change:
///   Hits->Draw("driftTrue:driftTime>>p(60)","","prof");  then fit a pol5.
inline double radiusFromTime(double t) {
  static constexpr double a[6] = {
      -2.869485e-11,   // t^5
       1.357186e-08,   // t^4
      -2.265298e-06,   // t^3
       1.547425e-04,   // t^2
       4.665418e-02,   // t^1
      -2.483804e-02    // t^0
  };
  double r = a[0];
  for (int i = 1; i < 6; ++i) r = r * t + a[i];
  return std::clamp(r, 0.0, kStrawRadius);
}

/// The naive linear inversion, kept for comparison: it carries the outward bias
/// described above. Use radiusFromTime() in the reco.
inline double radiusFromTimeLinear(double t) { return t * kDriftVelocity; }

/// Maximum drift time (a hit later than this is not from this straw).
inline double maxDriftTime() { return timeFromRadius(kStrawRadius); }   // ~200 ns

// ---------------------------------------------------------------------------
/// Simulate what the TDC sees for one track segment through the gas.
///
/// @param entry,exit  the Geant4 step end points, in the straw's LOCAL frame,
///                    where the wire is the z-axis (so the distance to the wire
///                    is just sqrt(x^2 + y^2)).
/// @param rng         random engine (the caller owns it; one per thread).
/// @return the drift time in ns, or -1 if the track left no clusters.
///
/// The algorithm, mirroring the physics:
///   1. sample N ~ Poisson(density x path length) cluster positions, uniformly
///      along the step;
///   2. for each, the perpendicular distance to the wire -> a drift time;
///   3. add diffusion (grows with drift distance) and TDC jitter;
///   4. the TDC fires on the EARLIEST -> take the minimum.
// ---------------------------------------------------------------------------
template <typename RNG>
inline double simulateDriftTime(const double entry[3], const double exit[3],
                                RNG& rng) {
  const double dx = exit[0] - entry[0];
  const double dy = exit[1] - entry[1];
  const double dz = exit[2] - entry[2];
  const double pathLen = std::sqrt(dx*dx + dy*dy + dz*dz);
  if (pathLen <= 0.0) return -1.0;

  // 1. how many primary clusters?
  std::poisson_distribution<int> nClusters(kClusterDensity * pathLen);
  const int n = nClusters(rng);
  if (n == 0) return -1.0;          // no ionisation -> no hit. A real inefficiency.

  std::uniform_real_distribution<double> along(0.0, 1.0);
  std::normal_distribution<double>       gauss(0.0, 1.0);

  // 2. the CLOSEST cluster sets the arrival: find the minimum radius over the
  //    sampled clusters. This is pure geometry + ionisation statistics, and it
  //    is the term that makes the resolution radius-dependent (and biased: the
  //    nearest cluster is generally NOT at the DOCA, so r_cluster >= r_doca,
  //    always -- a real effect that a calibrated r-t relation must absorb).
  double rMin = std::numeric_limits<double>::max();

  for (int i = 0; i < n; ++i) {
    const double s  = along(rng);
    const double cx = entry[0] + s * dx;
    const double cy = entry[1] + s * dy;
    const double r  = std::sqrt(cx*cx + cy*cy);
    if (r > kStrawRadius) continue;               // outside the gas
    rMin = std::min(rMin, r);
  }
  if (rMin == std::numeric_limits<double>::max()) return -1.0;

  // 3. that cluster's electrons drift in; diffusion (growing as sqrt r) and the
  //    TDC's jitter smear the ONE arrival time that is actually recorded.
  //    (Smearing every cluster and then taking the minimum would be wrong: it
  //    biases the time low by the min-of-N-Gaussians, an artefact of the model
  //    rather than of the detector.)
  double t = timeFromRadius(rMin);
  const double sigmaDiffusion_mm = kDiffusion * std::sqrt(rMin);
  t += (sigmaDiffusion_mm / kDriftVelocity) * gauss(rng);     // mm -> ns
  t += kTdcJitter * gauss(rng);

  return std::max(t, 0.0);           // a negative time is unphysical
}

/// The TRUE distance of closest approach of the segment to the wire (the z-axis
/// in the local frame). Diagnostic only -- the reco must never see this.
inline double trueDoca(const double entry[3], const double exit[3]) {
  // closest approach of a 2D segment (in x,y) to the origin
  const double dx = exit[0] - entry[0];
  const double dy = exit[1] - entry[1];
  const double len2 = dx*dx + dy*dy;
  if (len2 <= 0.0) return std::sqrt(entry[0]*entry[0] + entry[1]*entry[1]);

  double s = -(entry[0]*dx + entry[1]*dy) / len2;
  s = std::clamp(s, 0.0, 1.0);
  const double cx = entry[0] + s * dx;
  const double cy = entry[1] + s * dy;
  return std::sqrt(cx*cx + cy*cy);
}

}  // namespace strawdrift
