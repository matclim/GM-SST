#pragma once
// reco/ShipRecoContext.hpp — the ACTS contexts carried through the chain.
#include <any>
#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <Acts/Utilities/CalibrationContext.hpp>

namespace shipreco {
struct Contexts {
  Acts::GeometryContext        gctx{std::any{}};
  Acts::MagneticFieldContext   mctx{std::any{}};
  Acts::CalibrationContext     cctx{std::any{}};
};
}  // namespace shipreco
