// reco: ShipTrackFitter.cpp   (ACTS main @ b660c71)
#include "ShipTrackFitter.hpp"

#include <Acts/Propagator/DirectNavigator.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/TrackFitting/GainMatrixSmoother.hpp>
#include <Acts/TrackFitting/GainMatrixUpdater.hpp>
#include <Acts/TrackFitting/KalmanFitter.hpp>
#include <Acts/Utilities/Logger.hpp>

namespace shipreco {

using Stepper    = Acts::EigenStepper<>;
using Navigator  = Acts::DirectNavigator;
using Propagator = Acts::Propagator<Stepper, Navigator>;
using Fitter     = Acts::KalmanFitter<Propagator, Acts::VectorMultiTrajectory>;

namespace {
// Maps a SourceLink -> its straw surface (by geoId).
struct StrawAccessor {
  const ShipStrawGeometry* geo{nullptr};
  const Acts::Surface* operator()(const Acts::SourceLink& sl) const {
    return geo->surfaceForGeoId(sl.get<ShipSourceLink>().geometryId());
  }
};
}  // namespace

struct ShipTrackFitter::Impl {
  const ShipStrawGeometry* geom;
  const MeasStore*         store;
  std::unique_ptr<Fitter>  fitter;
  Acts::GainMatrixUpdater  updater;
  Acts::GainMatrixSmoother smoother;

  Impl(std::shared_ptr<const Acts::MagneticFieldProvider> field,
       const ShipStrawGeometry& g, const MeasStore& s)
      : geom(&g), store(&s) {
    Stepper stepper(std::move(field));
    Navigator navigator(Acts::getDefaultLogger("DirectNav", Acts::Logging::WARNING));
    Propagator prop(std::move(stepper), std::move(navigator),
                    Acts::getDefaultLogger("Propagator", Acts::Logging::WARNING));
    fitter = std::make_unique<Fitter>(
        std::move(prop),
        Acts::getDefaultLogger("KalmanFitter", Acts::Logging::WARNING));
  }
};

ShipTrackFitter::ShipTrackFitter(
    std::shared_ptr<const Acts::MagneticFieldProvider> field,
    const ShipStrawGeometry& geom, const MeasStore& store)
    : m_impl(std::make_unique<Impl>(std::move(field), geom, store)) {}

ShipTrackFitter::~ShipTrackFitter() = default;

bool ShipTrackFitter::fit(const Contexts& ctx,
                          const std::vector<Acts::SourceLink>& sourceLinks,
                          const Acts::BoundTrackParameters& start,
                          const std::vector<const Acts::Surface*>& sSequence,
                          TrackContainer& tracks) const {
  using Traj = Acts::VectorMultiTrajectory;

  StrawAccessor accessor{m_impl->geom};
  ShipCalibrator calibrator{m_impl->store};

  Acts::KalmanFitterExtensions<Traj> ext;
  ext.calibrator.connect<&ShipCalibrator::calibrate<Traj>>(&calibrator);
  ext.updater.connect<&Acts::GainMatrixUpdater::operator()<Traj>>(&m_impl->updater);
  ext.smoother.connect<&Acts::GainMatrixSmoother::operator()<Traj>>(&m_impl->smoother);
  ext.surfaceAccessor.connect<&StrawAccessor::operator()>(&accessor);

  Acts::PropagatorPlainOptions pOpts(ctx.gctx, ctx.mctx);

  Acts::KalmanFitterOptions<Traj> kfOpts(
      ctx.gctx, ctx.mctx, ctx.cctx, ext, pOpts,
      /*referenceSurface=*/&start.referenceSurface(),
      /*multipleScattering=*/true, /*energyLoss=*/true,
      /*reversedFiltering=*/false);

  auto res = m_impl->fitter->fit(sourceLinks.begin(), sourceLinks.end(),
                                 start, kfOpts, sSequence, tracks);
  return res.ok();
}

}  // namespace shipreco
