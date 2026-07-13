// ShipFieldMap.cpp
// Reads a 3-D tabulated magnetic field from a text file or a ROOT file
// ("Range"+"Data" TTrees) and evaluates it by trilinear interpolation.

#include "ShipFieldMap.h"

#include "G4SystemOfUnits.hh"
#include "globals.hh"

// ROOT (already a project dependency: Core, RIO, Tree)
#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <cctype>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
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
// File loader — dispatch on extension
// ─────────────────────────────────────────────────────────────────────────────
void ShipFieldMap::loadFile(const std::string& file) {
    // Lower-case the extension so ".ROOT" is handled too.
    std::string ext;
    if (auto dot = file.find_last_of('.'); dot != std::string::npos)
        ext = file.substr(dot);
    std::transform(ext.begin(), ext.end(), ext.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    if (ext == ".root")
        loadRootFile(file);
    else
        loadTextFile(file);
}

// ─────────────────────────────────────────────────────────────────────────────
// Shared: allocate flat B arrays (zero-filled) for the current grid geometry.
// Requires m_nX/m_nY/m_nZ to be set. Returns the total point count.
// ─────────────────────────────────────────────────────────────────────────────
void ShipFieldMap::allocateGrid() {
    const std::size_t expected =
        static_cast<std::size_t>(m_nX) * m_nY * m_nZ;
    m_Bx.assign(expected, 0.f);
    m_By.assign(expected, 0.f);
    m_Bz.assign(expected, 0.f);
}

// ─────────────────────────────────────────────────────────────────────────────
// Shared: bin one sample into the flat arrays by nearest grid index. This is
// robust against small floating-point drift and against the data points being
// listed in any order. Requires the grid geometry + arrays to be set up.
// ─────────────────────────────────────────────────────────────────────────────
// ─────────────────────────────────────────────────────────────────────────────
// UNITS. The SHiP ROOT field map stores POSITIONS IN CENTIMETRES (Range bounds
// and the Data x/y/z alike); the field values are in Tesla. Everything else in
// this class -- and in Geant4 -- is millimetres, so every coordinate read from
// the map is converted on the way in. Reading the map as mm shrinks the field
// region 10x (to +-500 mm instead of +-5000 mm), which silently un-bends any
// track more than half a metre off axis and destroys its momentum measurement.
static constexpr double kCmToMm = 10.0;

void ShipFieldMap::placeSample(double x, double y, double z,
                               double bx, double by, double bz) {
    const int ix = (m_dx > 0.0) ? static_cast<int>(std::lround((x - m_xMin) / m_dx)) : 0;
    const int iy = (m_dy > 0.0) ? static_cast<int>(std::lround((y - m_yMin) / m_dy)) : 0;
    const int iz = (m_dz > 0.0) ? static_cast<int>(std::lround((z - m_zMin) / m_dz)) : 0;
    if (ix < 0 || ix >= m_nX) return;
    if (iy < 0 || iy >= m_nY) return;
    if (iz < 0 || iz >= m_nZ) return;
    const std::size_t i = idx(ix, iy, iz);
    m_Bx[i] = static_cast<float>(bx);
    m_By[i] = static_cast<float>(by);
    m_Bz[i] = static_cast<float>(bz);
}

// ─────────────────────────────────────────────────────────────────────────────
// ROOT loader
//
// Reads a ROOT file with two TTrees:
//   "Range" (1 entry) : xMin,xMax,dx, yMin,yMax,dy, zMin,zMax,dz   (Float_t)
//   "Data"  (N entries): x,y,z, Bx,By,Bz                           (Float_t)
//
// Positions are taken to be in mm (magnet frame) and field values in Tesla,
// matching the text format. The grid dimensions are computed from "Range";
// if "Range" is absent, they are inferred from the unique coordinates in
// "Data" (same fallback the text loader uses).
// ─────────────────────────────────────────────────────────────────────────────
void ShipFieldMap::loadRootFile(const std::string& file) {
    std::unique_ptr<TFile> f(TFile::Open(file.c_str(), "READ"));
    if (!f || f->IsZombie()) {
        G4cerr << "[ShipFieldMap] ERROR: cannot open ROOT file '" << file << "'\n";
        return;
    }

    auto* dataTree = dynamic_cast<TTree*>(f->Get("Data"));
    if (!dataTree) {
        G4cerr << "[ShipFieldMap] ERROR: no 'Data' TTree in '" << file << "'\n";
        return;
    }

    // ── Grid geometry: prefer the "Range" tree ────────────────────────────────
    auto* rangeTree = dynamic_cast<TTree*>(f->Get("Range"));
    if (rangeTree && rangeTree->GetEntries() > 0) {
        Float_t xMin=0, xMax=0, dx=0, yMin=0, yMax=0, dy=0, zMin=0, zMax=0, dz=0;
        rangeTree->SetBranchAddress("xMin", &xMin);
        rangeTree->SetBranchAddress("xMax", &xMax);
        rangeTree->SetBranchAddress("dx",   &dx);
        rangeTree->SetBranchAddress("yMin", &yMin);
        rangeTree->SetBranchAddress("yMax", &yMax);
        rangeTree->SetBranchAddress("dy",   &dy);
        rangeTree->SetBranchAddress("zMin", &zMin);
        rangeTree->SetBranchAddress("zMax", &zMax);
        rangeTree->SetBranchAddress("dz",   &dz);
        rangeTree->GetEntry(0);

        // Map coordinates are in CENTIMETRES -> convert to mm.
        m_xMin = xMin*kCmToMm; m_xMax = xMax*kCmToMm; m_dx = dx*kCmToMm;
        m_yMin = yMin*kCmToMm; m_yMax = yMax*kCmToMm; m_dy = dy*kCmToMm;
        m_zMin = zMin*kCmToMm; m_zMax = zMax*kCmToMm; m_dz = dz*kCmToMm;

        m_nX = (m_dx > 0.0) ? static_cast<int>(std::lround((m_xMax - m_xMin) / m_dx)) + 1 : 1;
        m_nY = (m_dy > 0.0) ? static_cast<int>(std::lround((m_yMax - m_yMin) / m_dy)) + 1 : 1;
        m_nZ = (m_dz > 0.0) ? static_cast<int>(std::lround((m_zMax - m_zMin) / m_dz)) + 1 : 1;
    } else {
        // Fallback: no Range tree — infer the grid from the Data coordinates.
        G4cout << "[ShipFieldMap] No 'Range' tree; inferring grid from 'Data'.\n";
        Float_t x=0, y=0, z=0;
        dataTree->SetBranchAddress("x", &x);
        dataTree->SetBranchAddress("y", &y);
        dataTree->SetBranchAddress("z", &z);
        std::set<double> xs, ys, zs;
        const Long64_t nEnt = dataTree->GetEntries();
        for (Long64_t e = 0; e < nEnt; ++e) {
            dataTree->GetEntry(e);
            xs.insert(x*kCmToMm); ys.insert(y*kCmToMm); zs.insert(z*kCmToMm);
        }
        if (xs.empty()) {
            G4cerr << "[ShipFieldMap] ERROR: 'Data' tree is empty in '" << file << "'\n";
            return;
        }
        m_nX = static_cast<int>(xs.size());
        m_nY = static_cast<int>(ys.size());
        m_nZ = static_cast<int>(zs.size());
        m_xMin = *xs.begin();  m_xMax = *xs.rbegin();
        m_yMin = *ys.begin();  m_yMax = *ys.rbegin();
        m_zMin = *zs.begin();  m_zMax = *zs.rbegin();
        m_dx = (m_nX > 1) ? (m_xMax - m_xMin) / (m_nX - 1) : 0.0;
        m_dy = (m_nY > 1) ? (m_yMax - m_yMin) / (m_nY - 1) : 0.0;
        m_dz = (m_nZ > 1) ? (m_zMax - m_zMin) / (m_nZ - 1) : 0.0;
        dataTree->ResetBranchAddresses();
    }

    const std::size_t expected =
        static_cast<std::size_t>(m_nX) * m_nY * m_nZ;
    const Long64_t nEnt = dataTree->GetEntries();
    if (static_cast<std::size_t>(nEnt) != expected) {
        G4cerr << "[ShipFieldMap] WARNING: " << nEnt << " data entries but "
               << "inferred grid has " << expected << " points ("
               << m_nX << " x " << m_nY << " x " << m_nZ
               << "). Missing points will read as B = 0.\n";
    }

    allocateGrid();

    // ── Read the field samples ────────────────────────────────────────────────
    Float_t x=0, y=0, z=0, bx=0, by=0, bz=0;
    dataTree->SetBranchAddress("x",  &x);
    dataTree->SetBranchAddress("y",  &y);
    dataTree->SetBranchAddress("z",  &z);
    dataTree->SetBranchAddress("Bx", &bx);
    dataTree->SetBranchAddress("By", &by);
    dataTree->SetBranchAddress("Bz", &bz);

    for (Long64_t e = 0; e < nEnt; ++e) {
        dataTree->GetEntry(e);
        placeSample(x*kCmToMm, y*kCmToMm, z*kCmToMm, bx, by, bz);  // cm -> mm
    }
    dataTree->ResetBranchAddresses();

    m_valid = true;

    G4cout << std::fixed << std::setprecision(1);
    G4cout << "[ShipFieldMap] Loaded ROOT map '" << file << "'\n"
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
// Text loader
// ─────────────────────────────────────────────────────────────────────────────
void ShipFieldMap::loadTextFile(const std::string& file) {
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
    allocateGrid();
    for (const auto& r : rows)
        placeSample(r[0], r[1], r[2], r[3], r[4], r[5]);

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
