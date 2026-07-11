#pragma once
// reco/field/ActsShipField.hpp
#include <array>
#include <memory>
#include <string>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>

namespace shipreco {
/// Interpolated ACTS field from a SHiP ROOT map (Range+Data trees).
/// originMm : position (mm, ACTS global frame) of the map's local (0,0,0).
/// bScaleToTesla : 1.0 if file is Tesla, 0.1 if kilogauss.
std::shared_ptr<const Acts::MagneticFieldProvider>
makeShipFieldFromRootMap(const std::string& rootFile,
                         std::array<double, 3> originMm = {0.0, 0.0, 0.0},
                         double bScaleToTesla = 1.0);
}  // namespace shipreco
