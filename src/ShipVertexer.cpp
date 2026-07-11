// reco: ShipVertexer.cpp   (ACTS main @ b660c71)
//
// Single-vertex Billoir fit. CRITICAL: the Billoir algorithm expects tracks in
// the PERIGEE parameterization *with respect to the vertex* (d0, z0 are impact
// parameters relative to it) and linearizes around the vertex estimate. Feeding
// it parameters referenced on a surface far from the vertex breaks the
// linearization (non-convergence, or non-normal chi2). So we first propagate
// every track to a PerigeeSurface placed at the seed vertex, then fit those.
#include "ShipVertexer.hpp"

#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/VoidNavigator.hpp>
#include <Acts/Surfaces/PerigeeSurface.hpp>
#include <Acts/Surfaces/Surface.hpp>
#include <Acts/Utilities/Logger.hpp>
#include <Acts/Vertexing/FullBilloirVertexFitter.hpp>
#include <Acts/Vertexing/HelicalTrackLinearizer.hpp>
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
  if (tracks.size() < 2) return std::nullopt;

  // ---- 1. re-express every track at a perigee AT the seed vertex -----------
  auto perigee = Acts::Surface::makeShared<Acts::PerigeeSurface>(seedPos);
  Acts::PropagatorPlainOptions popts(ctx.gctx, ctx.mctx);

  std::vector<Acts::BoundTrackParameters> atPerigee;   // must outlive the fit
  atPerigee.reserve(tracks.size());
  for (const auto& t : tracks) {
    auto r = m_impl->propagator->propagateToSurface(t, *perigee, popts);
    if (!r.ok()) continue;
    const auto& p = *r;
    if (!p.parameters().allFinite()) continue;
    if (p.covariance() && !p.covariance()->allFinite()) continue;
    atPerigee.push_back(p);
  }
  if (atPerigee.size() < 2) return std::nullopt;

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
  if (!res.ok()) return std::nullopt;

  const Acts::Vertex& v = *res;
  VertexResult out;
  out.position   = v.position();
  out.covariance = v.covariance();
  out.nTracks    = static_cast<int>(atPerigee.size());
  return out;
}

}  // namespace shipreco
