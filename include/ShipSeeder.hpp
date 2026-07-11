#pragma once
// reco: ShipSeeder.hpp — self-contained seed for one track.
//
// Direction: least-squares line through ALL upstream (pre-magnet) hits.
// Momentum:  from the upstream->downstream y-deflection through the magnet,
//            p[GeV] = 0.3 * kBL[T·m] / |Δθ_y|.  No fixed guess, no truth.
#include <optional>

#include <Acts/EventData/BoundTrackParameters.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include "ShipRecoContext.hpp"
#include "ShipHitReader.hpp"

namespace shipreco {

/// pFallbackMeV is used only if the deflection can't be measured (e.g. missing
/// downstream hits or |Δθ|~0). charge: -1 for mu-. kBL: effective field
/// integral (T·m); default matches the SHiP spectrometer (~2.5 mrad at 10 GeV).
std::optional<Acts::BoundTrackParameters>
seedFromUpstream(const Contexts& ctx, const RawEvent& ev,
                 const Acts::Surface& reference,
                 double pFallbackMeV, double charge,
                 int upstreamA = 0, int upstreamB = 1,
                 int downstreamA = 2, int downstreamB = 3,
                 double kBL_Tm = 0.083);

}  // namespace shipreco
