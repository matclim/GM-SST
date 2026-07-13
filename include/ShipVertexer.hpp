#pragma once
#include <string>
#include <vector>
// reco: ShipVertexer.hpp — single-vertex Billoir fit over a set of tracks.
#include <memory>
#include <optional>
#include <vector>

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/EventData/BoundTrackParameters.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>

#include "ShipRecoContext.hpp"

namespace shipreco {

/// Why a vertex fit failed. Distinguishes "the tracks never reached the perigee"
/// (a propagation problem) from "Billoir diverged" (a fit problem) -- they need
/// completely different fixes, and lumping them together hides which is which.
enum class VertexFail {
  None            = 0,
  TooFewTracks    = 1,   // <2 tracks given
  PropagationFail = 2,   // <2 tracks reached the perigee surface
  BilloirFail     = 3    // the fit itself returned an error
};

/// Per-track propagation diagnostics. ACTS tells us *why* a propagation failed
/// (path limit, step limit, stepper failure, surface unreachable, ...); throwing
/// that away and inferring the cause from the outside wastes a lot of time.
struct PropDiag {
  bool        ok        {false};
  std::string error     {};      // ACTS error message, if it failed
  double      pathLen   {0.0};   // mm actually propagated
  bool        backward  {false}; // did we go against the momentum?
  double      distToVtx {0.0};   // |seed - track position| before propagating
};

struct VertexResult {
  Acts::Vector3       position{Acts::Vector3::Zero()};
  Acts::SquareMatrix<3> covariance{Acts::SquareMatrix<3>::Zero()};
  int nTracks{0};
  int nPropagated{0};        // how many reached the perigee
  VertexFail fail{VertexFail::None};
  std::vector<PropDiag> propDiag;   // one per input track
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
