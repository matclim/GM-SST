// reco: ShipVertexer.cpp   (ACTS main @ b660c71)
//
// Single-vertex Billoir fit, for any N >= 2 tracks.
//
// Billoir expects tracks in the PERIGEE parameterisation *with respect to the
// vertex* (d0, z0 are impact parameters relative to it) and linearises around
// the vertex estimate. Parameters referenced on a surface far from the vertex
// break that linearisation (frozen fit, or non-normal chi2).
//
// We therefore re-express every track at its 3D point of closest approach to
// the seed vertex, using Acts::ImpactPointEstimator -- the supported tool for
// exactly this. An earlier hand-rolled version ("propagate to a PerigeeSurface
// placed at the seed") failed for long back-extrapolations: the propagator
// lands NEAR the perigee but not exactly at the PCA, and LineSurface::
// globalToLocal then rejects the point with GlobalPositionNotOnSurface. For LLP
// decays 25-50 m upstream of the tracker that killed ~95% of vertices.
#include "ShipVertexer.hpp"

#include <Acts/Definitions/Direction.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/VoidNavigator.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Vertexing/FullBilloirVertexFitter.hpp>
#include <Acts/Vertexing/HelicalTrackLinearizer.hpp>
#include <Acts/Vertexing/ImpactPointEstimator.hpp>
#include <Acts/Vertexing/TrackAtVertex.hpp>
#include <Acts/Vertexing/Vertex.hpp>
#include <Acts/Vertexing/VertexingOptions.hpp>

namespace shipreco {

using Stepper    = Acts::EigenStepper<>;
using Propagator = Acts::Propagator<Stepper, Acts::VoidNavigator>;
using Linearizer = Acts::HelicalTrackLinearizer;

struct ShipVertexer::Impl {
  std::shared_ptr<const Acts::MagneticFieldProvider> field;
  std::shared_ptr<Propagator>                        propagator;
  std::unique_ptr<Linearizer>                        linearizer;
  std::unique_ptr<Acts::FullBilloirVertexFitter>     fitter;
  std::unique_ptr<Acts::ImpactPointEstimator>        ipEst;

  explicit Impl(std::shared_ptr<const Acts::MagneticFieldProvider> f)
      : field(std::move(f)) {
    Stepper stepper(field);
    propagator = std::make_shared<Propagator>(
        std::move(stepper), Acts::VoidNavigator{},
        Acts::getDefaultLogger("VtxProp", Acts::Logging::WARNING));

    Linearizer::Config linCfg;
    linCfg.bField     = field;
    linCfg.propagator = propagator;
    linearizer = std::make_unique<Linearizer>(
        linCfg, Acts::getDefaultLogger("Linearizer", Acts::Logging::WARNING));

    // ImpactPointEstimator: finds each track's parameters AT the 3D point of
    // closest approach to the vertex. This is the supported way to perigee-
    // reference a track. Hand-rolling it as "propagate to a PerigeeSurface at
    // the seed" fails over long back-extrapolations: the propagator lands NEAR
    // the perigee but not exactly at the PCA, and LineSurface::globalToLocal
    // then rejects it with GlobalPositionNotOnSurface. (For a 25-50 m swim back
    // into the decay volume that killed ~95% of vertices.)
    Acts::ImpactPointEstimator::Config ipCfg(field, propagator);
    ipEst = std::make_unique<Acts::ImpactPointEstimator>(ipCfg);

    Acts::FullBilloirVertexFitter::Config fitCfg;
    fitCfg.extractParameters.connect<&Acts::InputTrack::extractParameters>();
    fitCfg.trackLinearizer.connect<&Linearizer::linearizeTrack>(linearizer.get());
    fitter = std::make_unique<Acts::FullBilloirVertexFitter>(
        fitCfg, Acts::getDefaultLogger("Billoir", Acts::Logging::WARNING));
  }
};

ShipVertexer::ShipVertexer(
    std::shared_ptr<const Acts::MagneticFieldProvider> field)
    : m_impl(std::make_unique<Impl>(std::move(field))) {}

ShipVertexer::~ShipVertexer() = default;

std::optional<VertexResult> ShipVertexer::fit(
    const Contexts& ctx,
    const std::vector<Acts::BoundTrackParameters>& tracks,
    const Acts::Vector3& seedPos) const {
  VertexResult out;
  out.fail = VertexFail::None;
  if (tracks.size() < 2) { out.fail = VertexFail::TooFewTracks; return out; }

  // ---- 1. re-express every track at its PCA to the seed vertex -------------
  // Billoir expects perigee parameters *with respect to the vertex* and
  // linearises around it; plane-referenced parameters from tens of metres away
  // are useless to it. ImpactPointEstimator does this properly.
  Acts::ImpactPointEstimator::State ipState{m_impl->field->makeCache(ctx.mctx)};

  std::vector<Acts::BoundTrackParameters> atPerigee;   // must outlive the fit
  atPerigee.reserve(tracks.size());

  for (const auto& t : tracks) {
    PropDiag diag;
    const Acts::Vector3 toVtx = seedPos - t.position(ctx.gctx);
    diag.backward  = (toVtx.dot(t.direction()) < 0);
    diag.distToVtx = toVtx.norm();

    auto r = m_impl->ipEst->estimate3DImpactParameters(
        ctx.gctx, ctx.mctx, t, seedPos, ipState);

    if (!r.ok()) {
      diag.error = r.error().message();
      out.propDiag.push_back(diag);
      continue;
    }
    const auto& p = *r;
    if (!p.parameters().allFinite()) {
      diag.error = "non-finite parameters at the PCA";
      out.propDiag.push_back(diag);
      continue;
    }
    if (p.covariance() && !p.covariance()->allFinite()) {
      diag.error = "non-finite covariance at the PCA";
      out.propDiag.push_back(diag);
      continue;
    }
    diag.ok = true;
    out.propDiag.push_back(diag);
    atPerigee.push_back(p);
  }

  out.nPropagated = static_cast<int>(atPerigee.size());
  if (atPerigee.size() < 2) { out.fail = VertexFail::PropagationFail; return out; }

  // ---- 2. Billoir fit on the perigee-parameterized tracks ------------------
  std::vector<Acts::InputTrack> inputs;
  inputs.reserve(atPerigee.size());
  for (const auto& p : atPerigee) inputs.emplace_back(&p);

  // Seed vertex only sets the linearization point; no constraint pull.
  Acts::Vertex seedVertex(seedPos);
  Acts::VertexingOptions opts(ctx.gctx, ctx.mctx, seedVertex,
                              /*useConstraintInFit=*/false);
  auto cache = m_impl->field->makeCache(ctx.mctx);

  auto res = m_impl->fitter->fit(inputs, opts, cache);
  if (!res.ok()) { out.fail = VertexFail::BilloirFail; return out; }

  const Acts::Vertex& v = *res;
  out.position   = v.position();
  out.covariance = v.covariance();
  out.nTracks    = static_cast<int>(atPerigee.size());
  out.fail       = VertexFail::None;
  return out;
}

}  // namespace shipreco
