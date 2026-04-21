// ShipFieldMap.cpp
// Reads a 3-D tabulated magnetic field from a text file and evaluates it by
// trilinear interpolation.

#include "ShipFieldMap.h"

#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

ShipFieldMap::ShipFieldMap(const std::string& file,
                           const G4ThreeVector& mapOriginG4)
    : m_mapOriginG4(mapOriginG4)
{
    loadFile(file);
}

// ─────────────────────────────────────────────────────────────────────────────
// File loader
// ─────────────────────────────────────────────────────────────────────────────
void ShipFieldMap::loadFile(const std::string& file) {
    std::ifstream in(file);
    if (!in) {
        G4cerr << "[ShipFieldMap] ERROR: cannot open file '" << file << "'\n";
        return;
    }

    // Collect all valid data rows into memory first; the grid geometry is
    // inferred from the unique X/Y/Z values that appear.
    std::vector<std::array<double, 6>> rows;
    rows.reserve(1 << 16);

    std::string line;
    std::size_t lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        // Trim leading whitespace
        auto pos = line.find_first_not_of(" \t\r\n");
        if (pos == std::string::npos)         continue;  // blank
        if (line[pos] == '#')                 continue;  // comment

        std::istringstream iss(line);
        double x, y, z, bx, by, bz;
        if (!(iss >> x >> y >> z >> bx >> by >> bz)) {
            // Header line ("X_SHiP Y_SHiP ...") or malformed — skip silently.
            continue;
        }
        rows.push_back({x, y, z, bx, by, bz});
    }

    if (rows.empty()) {
        G4cerr << "[ShipFieldMap] ERROR: no usable data rows in '" << file << "'\n";
        return;
    }

    // Infer grid extent by collecting unique coordinates.
    std::set<double> xs, ys, zs;
    for (const auto& r : rows) {
        xs.insert(r[0]);
        ys.insert(r[1]);
        zs.insert(r[2]);
    }
    m_nX = static_cast<int>(xs.size());
    m_nY = static_cast<int>(ys.size());
    m_nZ = static_cast<int>(zs.size());

    m_xMin = *xs.begin();   m_xMax = *xs.rbegin();
    m_yMin = *ys.begin();   m_yMax = *ys.rbegin();
    m_zMin = *zs.begin();   m_zMax = *zs.rbegin();

    m_dx = (m_nX > 1) ? (m_xMax - m_xMin) / (m_nX - 1) : 0.0;
    m_dy = (m_nY > 1) ? (m_yMax - m_yMin) / (m_nY - 1) : 0.0;
    m_dz = (m_nZ > 1) ? (m_zMax - m_zMin) / (m_nZ - 1) : 0.0;

    const std::size_t expected =
        static_cast<std::size_t>(m_nX) * m_nY * m_nZ;

    if (rows.size() != expected) {
        G4cerr << "[ShipFieldMap] WARNING: " << rows.size() << " data rows "
               << "but inferred grid has " << expected << " points ("
               << m_nX << " x " << m_nY << " x " << m_nZ
               << "). Missing points will read as B = 0.\n";
    }

    // Allocate + fill the flat arrays. We bin each row by nearest grid index
    // to be robust against small floating-point drift in the coordinates.
    m_Bx.assign(expected, 0.f);
    m_By.assign(expected, 0.f);
    m_Bz.assign(expected, 0.f);

    for (const auto& r : rows) {
        const int ix = (m_dx > 0.0) ? static_cast<int>(std::lround((r[0] - m_xMin) / m_dx)) : 0;
        const int iy = (m_dy > 0.0) ? static_cast<int>(std::lround((r[1] - m_yMin) / m_dy)) : 0;
        const int iz = (m_dz > 0.0) ? static_cast<int>(std::lround((r[2] - m_zMin) / m_dz)) : 0;
        if (ix < 0 || ix >= m_nX) continue;
        if (iy < 0 || iy >= m_nY) continue;
        if (iz < 0 || iz >= m_nZ) continue;
        const std::size_t i = idx(ix, iy, iz);
        m_Bx[i] = static_cast<float>(r[3]);
        m_By[i] = static_cast<float>(r[4]);
        m_Bz[i] = static_cast<float>(r[5]);
    }

    m_valid = true;

    G4cout << std::fixed << std::setprecision(1);
    G4cout << "[ShipFieldMap] Loaded '" << file << "'\n"
           << "               grid "  << m_nX << " x " << m_nY << " x " << m_nZ
           << " = " << expected << " points\n"
           << "               X in [" << m_xMin << ", " << m_xMax << "] mm, "
           << "dx = " << m_dx << " mm\n"
           << "               Y in [" << m_yMin << ", " << m_yMax << "] mm, "
           << "dy = " << m_dy << " mm\n"
           << "               Z in [" << m_zMin << ", " << m_zMax << "] mm, "
           << "dz = " << m_dz << " mm\n"
           << "               Map origin in Geant4 frame: ("
           << m_mapOriginG4.x() / CLHEP::mm << ", "
           << m_mapOriginG4.y() / CLHEP::mm << ", "
           << m_mapOriginG4.z() / CLHEP::mm << ") mm\n";
    G4cout << std::defaultfloat;
}

// ─────────────────────────────────────────────────────────────────────────────
// Geant4 field evaluation
// ─────────────────────────────────────────────────────────────────────────────
void ShipFieldMap::GetFieldValue(const G4double point[4], G4double* Bfield) const {
    // Initialise to zero. If the query point is outside the loaded grid,
    // this is what the caller sees.
    Bfield[0] = 0.0;
    Bfield[1] = 0.0;
    Bfield[2] = 0.0;

    if (!m_valid) return;

    // Transform from Geant4 world frame to map frame, then from internal units
    // to mm (the grid is stored in mm).
    const double x = (point[0] - m_mapOriginG4.x()) / CLHEP::mm;
    const double y = (point[1] - m_mapOriginG4.y()) / CLHEP::mm;
    const double z = (point[2] - m_mapOriginG4.z()) / CLHEP::mm;

    if (x < m_xMin || x > m_xMax) return;
    if (y < m_yMin || y > m_yMax) return;
    if (z < m_zMin || z > m_zMax) return;

    // Fractional grid indices
    const double fx = (m_dx > 0.0) ? (x - m_xMin) / m_dx : 0.0;
    const double fy = (m_dy > 0.0) ? (y - m_yMin) / m_dy : 0.0;
    const double fz = (m_dz > 0.0) ? (z - m_zMin) / m_dz : 0.0;

    // Lower-corner integer indices, clamped to [0, n-2].
    int ix = static_cast<int>(fx); if (ix >= m_nX - 1) ix = m_nX - 2; if (ix < 0) ix = 0;
    int iy = static_cast<int>(fy); if (iy >= m_nY - 1) iy = m_nY - 2; if (iy < 0) iy = 0;
    int iz = static_cast<int>(fz); if (iz >= m_nZ - 1) iz = m_nZ - 2; if (iz < 0) iz = 0;

    const double tx = fx - ix;
    const double ty = fy - iy;
    const double tz = fz - iz;

    // Trilinear interpolation
    auto interp = [&](const std::vector<float>& B) -> double {
        const float b000 = B[idx(ix,   iy,   iz  )];
        const float b100 = B[idx(ix+1, iy,   iz  )];
        const float b010 = B[idx(ix,   iy+1, iz  )];
        const float b110 = B[idx(ix+1, iy+1, iz  )];
        const float b001 = B[idx(ix,   iy,   iz+1)];
        const float b101 = B[idx(ix+1, iy,   iz+1)];
        const float b011 = B[idx(ix,   iy+1, iz+1)];
        const float b111 = B[idx(ix+1, iy+1, iz+1)];

        const double b00 = b000 * (1.0 - tx) + b100 * tx;
        const double b10 = b010 * (1.0 - tx) + b110 * tx;
        const double b01 = b001 * (1.0 - tx) + b101 * tx;
        const double b11 = b011 * (1.0 - tx) + b111 * tx;

        const double b0 = b00 * (1.0 - ty) + b10 * ty;
        const double b1 = b01 * (1.0 - ty) + b11 * ty;

        return b0 * (1.0 - tz) + b1 * tz;
    };

    // File values are in Tesla → multiply by CLHEP::tesla to go to Geant4 units.
    Bfield[0] = interp(m_Bx) * CLHEP::tesla;
    Bfield[1] = interp(m_By) * CLHEP::tesla;
    Bfield[2] = interp(m_Bz) * CLHEP::tesla;
}
