#pragma once
// ShipFieldMap.h
// G4MagneticField that reads a 3-D field map from a text file and does
// trilinear interpolation.
//
// Expected file format (whitespace-separated, one data point per line):
//
//   X_SHiP   Y_SHiP   Z_SHiP   BX_SHiP   BY_SHiP   BZ_SHiP    ← optional header
//   -5000    -5000    -7500    0.00150...  9.09e-05  5.49e-04
//   ...
//
//   - Positions (X_SHiP, Y_SHiP, Z_SHiP) are assumed to be in mm, in the
//     MAGNET's own coordinate system (origin at the magnet centre).
//   - Field values (BX_SHiP, BY_SHiP, BZ_SHiP) are assumed to be in Tesla.
//   - Lines starting with '#' are treated as comments.
//   - A single header line that fails to parse as six doubles is silently
//     skipped, so the file above (with its text header) works out of the box.
//
// The caller supplies the position of the map's origin in the Geant4 world
// frame via `mapOriginG4`. The world frame in the current project is centred
// at lab z = 31000 mm, and the field-map origin is at lab z = 31500 mm, so the
// correct setting is G4ThreeVector(0, 0, 500*mm).

#include "G4MagneticField.hh"
#include "G4ThreeVector.hh"

#include <string>
#include <vector>
#include <cstddef>

class ShipFieldMap : public G4MagneticField {
public:
    /// @param file          path to the text field-map file
    /// @param mapOriginG4   position of the map's origin in the Geant4 frame
    ShipFieldMap(const std::string& file,
                 const G4ThreeVector& mapOriginG4);

    // Geant4 interface — B in Geant4 internal units, point in world coords.
    void GetFieldValue(const G4double point[4], G4double* Bfield) const override;

    // Diagnostics
    bool   valid() const { return m_valid; }
    int    nX()    const { return m_nX; }
    int    nY()    const { return m_nY; }
    int    nZ()    const { return m_nZ; }
    double xMin()  const { return m_xMin; }
    double xMax()  const { return m_xMax; }
    double yMin()  const { return m_yMin; }
    double yMax()  const { return m_yMax; }
    double zMin()  const { return m_zMin; }
    double zMax()  const { return m_zMax; }

private:
    // Flat index into a size-(nX*nY*nZ) array, Z-major, then Y, then X.
    inline std::size_t idx(int ix, int iy, int iz) const {
        return (static_cast<std::size_t>(iz) * m_nY + iy) * m_nX + ix;
    }

    void loadFile(const std::string& file);

    // Field components (Tesla), flat arrays of size nX*nY*nZ
    std::vector<float> m_Bx, m_By, m_Bz;

    // Grid geometry (map frame, mm)
    int    m_nX{0}, m_nY{0}, m_nZ{0};
    double m_xMin{0}, m_yMin{0}, m_zMin{0};
    double m_xMax{0}, m_yMax{0}, m_zMax{0};
    double m_dx{0},   m_dy{0},   m_dz{0};

    // Map origin in the Geant4 frame (internal units)
    G4ThreeVector m_mapOriginG4;

    bool m_valid{false};
};
