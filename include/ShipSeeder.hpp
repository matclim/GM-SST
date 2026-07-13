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
/// integral [T·m]. NOTE: 0.49 T·m, re-derived after the field-map cm->mm
/// unit fix (the old 0.083 was tuned against a 10x-truncated field).
/// integral (T·m); 0.49 gives ~15 mrad at 10 GeV, matching the measured
/// muon deflection of -60 mm at station 3.
std::optional<Acts::BoundTrackParameters>
seedFromUpstream(const Contexts& ctx, const RawEvent& ev,
                 const Acts::Surface& reference,
                 double pFallbackMeV, double charge,
                 int upstreamA = 0, int upstreamB = 1,
                 int downstreamA = 2, int downstreamB = 3,
                 double kBL_Tm = 0.49);   // measured post units-fix (was 0.083)

}  // namespace shipreco
