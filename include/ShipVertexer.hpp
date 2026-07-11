#pragma once
// reco: ShipVertexer.hpp — single-vertex Billoir fit over a set of tracks.
#include <memory>
#include <optional>
#include <vector>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/EventData/BoundTrackParameters.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>

#include "ShipRecoContext.hpp"

namespace shipreco {

struct VertexResult {
  Acts::Vector3       position{Acts::Vector3::Zero()};
  Acts::SquareMatrix<3> covariance{Acts::SquareMatrix<3>::Zero()};
  int nTracks{0};
};

class ShipVertexer {
 public:
  explicit ShipVertexer(std::shared_ptr<const Acts::MagneticFieldProvider> field);
  ~ShipVertexer();

  /// Fit one common vertex from >=2 tracks. Returns nullopt on failure.
  /// seedPos: starting vertex estimate (e.g. tracks' closest approach), in the
  /// ACTS/world frame. Billoir linearizes around it — essential for a vertex
  /// far from the origin.
  std::optional<VertexResult> fit(
      const Contexts& ctx,
      const std::vector<Acts::BoundTrackParameters>& tracks,
      const Acts::Vector3& seedPos) const;

 private:
  struct Impl;
  std::unique_ptr<Impl> m_impl;
};

}  // namespace shipreco
