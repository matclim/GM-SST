#pragma once
// reco: ShipTrackFitter.hpp — KalmanFitter over a per-event straw-surface list.
#include <memory>
#include <vector>

#include <Acts/EventData/SourceLink.hpp>
#include <Acts/EventData/TrackContainer.hpp>
#include <Acts/EventData/BoundTrackParameters.hpp>
#include <Acts/EventData/VectorMultiTrajectory.hpp>
#include <Acts/EventData/VectorTrackContainer.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include "ShipRecoContext.hpp"
#include "ShipEDM.hpp"
#include "ShipStrawGeometry.hpp"

namespace shipreco {

using TrackContainer =
    Acts::TrackContainer<Acts::VectorTrackContainer, Acts::VectorMultiTrajectory,
                         std::shared_ptr>;

class ShipTrackFitter {
 public:
  ShipTrackFitter(std::shared_ptr<const Acts::MagneticFieldProvider> field,
                  const ShipStrawGeometry& geom, const MeasStore& store);
  ~ShipTrackFitter();

  /// Fit one track. `sSequence` = the straw surfaces the track visits, ordered
  /// along the trajectory (z). `sourceLinks` reference those surfaces by geoId.
  bool fit(const Contexts& ctx,
           const std::vector<Acts::SourceLink>& sourceLinks,
           const Acts::BoundTrackParameters& start,
           const std::vector<const Acts::Surface*>& sSequence,
           TrackContainer& tracks) const;

 private:
  struct Impl;
  std::unique_ptr<Impl> m_impl;
};

}  // namespace shipreco
