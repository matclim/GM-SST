#pragma once
// reco/edm/ShipEDM.hpp
// Measurement store, SourceLink, surface accessor and 1-D calibrator.
//
// >>> This file concentrates the ACTS API most likely to differ by version. <<<
// Every such line is marked // ACTS-API. If the fit fails to compile, look here
// first and mirror the equivalent call from your ACTS' Examples measurement
// calibrator (Examples/Algorithms/.../MeasurementCalibration).

#include <cstddef>
#include <vector>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/EventData/SourceLink.hpp>
#include <Acts/Geometry/GeometryIdentifier.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

#include "ShipGeometry.hpp"

namespace shipreco {

/// One 1-D straw-layer measurement: loc0 (mm) with resolution sigma0 (mm).
/// One straw measurement.
///
/// A real straw delivers a DRIFT TIME, which the reco inverts (via the
/// calibrated r-t relation in StrawDrift.h) into an UNSIGNED radius. It carries
/// no information about WHICH SIDE of the wire the track passed -- the
/// left/right ambiguity -- so `radius` is always >= 0, and the sign is resolved
/// by the FIT (see ShipCalibrator below), never by truth.
struct Meas {
  Acts::GeometryIdentifier geoId{};
  double radius{0.0};     // UNSIGNED drift radius, mm
  double sigma0{0.10};    // ~100 um for our 10 mm straws
  double truthLoc0{0.0};  // SIGNED truth drift: DIAGNOSTICS ONLY, never fitted
};

using MeasStore = std::vector<Meas>;

/// Concrete source link: points into the MeasStore and carries the surface id.
class ShipSourceLink {
 public:
  ShipSourceLink(Acts::GeometryIdentifier id, std::size_t idx)
      : m_geoId(id), m_index(idx) {}
  Acts::GeometryIdentifier geometryId() const { return m_geoId; }
  std::size_t index() const { return m_index; }
 private:
  Acts::GeometryIdentifier m_geoId;
  std::size_t              m_index;
};

/// Maps a SourceLink -> Surface for the fitter's surfaceAccessor extension.
struct ShipSurfaceAccessor {
  const ShipGeometry* geo{nullptr};
  const Acts::Surface* operator()(const Acts::SourceLink& sl) const {
    const auto& s = sl.get<ShipSourceLink>();          // ACTS-API: SourceLink::get<>
    return geo->surfaceForGeoId(s.geometryId());
  }
};

/// Turns a SourceLink into a calibrated 1-D measurement on a track state.
struct ShipCalibrator {
  const MeasStore* store{nullptr};

  template <typename traj_t>
  void calibrate(const Acts::GeometryContext& /*gctx*/,
                 const Acts::CalibrationContext& /*cctx*/,
                 const Acts::SourceLink& sl,
                 typename traj_t::TrackStateProxy ts) const {
    const auto& ssl = sl.get<ShipSourceLink>();
    const Meas& m = (*store)[ssl.index()];

    ts.setUncalibratedSourceLink(Acts::SourceLink{sl});                  // ACTS-API

    // ACTS-API cluster: allocate a size-1 calibrated measurement and fill it.
    ts.allocateCalibrated(1);
    // LEFT/RIGHT. The measurement is an UNSIGNED radius (a drift time carries
    // no side information), but eBoundLoc0 on a line surface is SIGNED. So the
    // fit must choose the side: take the sign of the PREDICTED loc0, i.e. the
    // side the track is currently believed to have passed. No truth involved.
    // It is only as good as the prediction -- which is why ShipTrackFitter
    // ITERATES: a first pass from a rough seed, then a refit from that result,
    // by which point the predicted signs are nearly all correct.
    const double predLoc0 = ts.predicted()[Acts::eBoundLoc0];
    const double loc0     = (predLoc0 < 0.0 ? -1.0 : +1.0) * m.radius;
    ts.template calibrated<1>() = Acts::Vector<1>(loc0);
    ts.template calibratedCovariance<1>() =
        Acts::SquareMatrix<1>(m.sigma0 * m.sigma0);

    // ACTS-API: declare that we measure eBoundLoc0. In some versions this is
    // setProjectorSubspaceIndices(...); in others setBoundSubspaceIndices(...).
    ts.setProjectorSubspaceIndices(
        std::array<Acts::BoundIndices, 1>{Acts::eBoundLoc0});
  }
};

}  // namespace shipreco
