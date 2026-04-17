// dump_field.C
// =============================================================================
// Diagnostic: read a field-map file, print the Bx profile along the beam axis
// (x=0, y=0) in the lab frame, and compute straight-line integrals of Bx over
// several z-windows. Useful for sanity-checking momentum reconstruction when
// p_reco disagrees with the truth.
//
// Usage:
//   root -l -b -q 'dump_field.C("fieldmaps/V21_2000A_B_full_fieldmap_SHiP.txt")'
// =============================================================================

#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

// Same field-map reader as in measure_momentum.C — see that file for details.
class FieldMap {
public:
    bool load(const std::string& fname, double originZ_lab_mm = 31500.0) {
        m_originZ = originZ_lab_mm;
        std::ifstream in(fname);
        if (!in) { std::cerr << "[FieldMap] Cannot open '" << fname << "'\n"; return false; }
        std::vector<std::array<double, 6>> rows;
        std::string line;
        while (std::getline(in, line)) {
            auto p = line.find_first_not_of(" \t\r\n");
            if (p == std::string::npos || line[p] == '#') continue;
            std::istringstream iss(line);
            double x, y, z, bx, by, bz;
            if (!(iss >> x >> y >> z >> bx >> by >> bz)) continue;
            rows.push_back({x, y, z, bx, by, bz});
        }
        if (rows.empty()) return false;

        std::set<double> xs, ys, zs;
        for (auto& r : rows) { xs.insert(r[0]); ys.insert(r[1]); zs.insert(r[2]); }
        m_nX = (int)xs.size(); m_nY = (int)ys.size(); m_nZ = (int)zs.size();
        m_xMin = *xs.begin();  m_xMax = *xs.rbegin();
        m_yMin = *ys.begin();  m_yMax = *ys.rbegin();
        m_zMin = *zs.begin();  m_zMax = *zs.rbegin();
        m_dx = (m_nX > 1) ? (m_xMax - m_xMin) / (m_nX - 1) : 0;
        m_dy = (m_nY > 1) ? (m_yMax - m_yMin) / (m_nY - 1) : 0;
        m_dz = (m_nZ > 1) ? (m_zMax - m_zMin) / (m_nZ - 1) : 0;

        const std::size_t N = (std::size_t)m_nX * m_nY * m_nZ;
        m_Bx.assign(N, 0); m_By.assign(N, 0); m_Bz.assign(N, 0);
        for (auto& r : rows) {
            int ix = (m_dx > 0) ? (int)std::lround((r[0] - m_xMin) / m_dx) : 0;
            int iy = (m_dy > 0) ? (int)std::lround((r[1] - m_yMin) / m_dy) : 0;
            int iz = (m_dz > 0) ? (int)std::lround((r[2] - m_zMin) / m_dz) : 0;
            if (ix < 0 || ix >= m_nX || iy < 0 || iy >= m_nY || iz < 0 || iz >= m_nZ) continue;
            std::size_t i = ((std::size_t)iz * m_nY + iy) * m_nX + ix;
            m_Bx[i] = (float)r[3]; m_By[i] = (float)r[4]; m_Bz[i] = (float)r[5];
        }
        m_valid = true;
        std::cout << "[FieldMap] " << fname << ": "
                  << m_nX << "x" << m_nY << "x" << m_nZ << " grid, "
                  << "X [" << m_xMin << "," << m_xMax << "] mm, "
                  << "Y [" << m_yMin << "," << m_yMax << "] mm, "
                  << "Z [" << m_zMin << "," << m_zMax << "] mm\n"
                  << "           lab-z origin = " << m_originZ << " mm\n\n";
        return true;
    }

    bool valid() const { return m_valid; }
    double Bx(double x, double y, double z_lab) const { return interp(x, y, z_lab - m_originZ, m_Bx); }
    double By(double x, double y, double z_lab) const { return interp(x, y, z_lab - m_originZ, m_By); }
    double Bz(double x, double y, double z_lab) const { return interp(x, y, z_lab - m_originZ, m_Bz); }

private:
    double interp(double x, double y, double z, const std::vector<float>& B) const {
        if (!m_valid) return 0.0;
        if (x < m_xMin || x > m_xMax || y < m_yMin || y > m_yMax || z < m_zMin || z > m_zMax)
            return 0.0;
        const double fx = (m_dx > 0) ? (x - m_xMin) / m_dx : 0;
        const double fy = (m_dy > 0) ? (y - m_yMin) / m_dy : 0;
        const double fz = (m_dz > 0) ? (z - m_zMin) / m_dz : 0;
        int ix = (int)fx; if (ix >= m_nX - 1) ix = m_nX - 2; if (ix < 0) ix = 0;
        int iy = (int)fy; if (iy >= m_nY - 1) iy = m_nY - 2; if (iy < 0) iy = 0;
        int iz = (int)fz; if (iz >= m_nZ - 1) iz = m_nZ - 2; if (iz < 0) iz = 0;
        const double tx = fx - ix, ty = fy - iy, tz = fz - iz;
        auto at = [&](int i, int j, int k) {
            return B[((std::size_t)k * m_nY + j) * m_nX + i];
        };
        const double b00 = at(ix,iy,iz)  *(1-tx) + at(ix+1,iy,iz)  *tx;
        const double b10 = at(ix,iy+1,iz)*(1-tx) + at(ix+1,iy+1,iz)*tx;
        const double b01 = at(ix,iy,iz+1)*(1-tx) + at(ix+1,iy,iz+1)*tx;
        const double b11 = at(ix,iy+1,iz+1)*(1-tx)+at(ix+1,iy+1,iz+1)*tx;
        const double b0  = b00*(1-ty) + b10*ty;
        const double b1  = b01*(1-ty) + b11*ty;
        return b0*(1-tz) + b1*tz;
    }
    std::vector<float> m_Bx, m_By, m_Bz;
    int    m_nX{0}, m_nY{0}, m_nZ{0};
    double m_xMin{0}, m_yMin{0}, m_zMin{0}, m_xMax{0}, m_yMax{0}, m_zMax{0};
    double m_dx{0}, m_dy{0}, m_dz{0}, m_originZ{31500};
    bool   m_valid{false};
};

// Straight-line integral of Bx along (x=0, y=0, z) from z1 to z2 (in lab mm).
double intBx_straight(const FieldMap& fmap, double z1, double z2, int N = 2000) {
    const double dz = (z2 - z1) / N;
    double sum = 0;
    for (int i = 0; i <= N; ++i) {
        const double z = z1 + i * dz;
        const double w = (i == 0 || i == N) ? 0.5 : 1.0;
        sum += w * fmap.Bx(0, 0, z);
    }
    return sum * dz * 1e-3;  // T·m
}

} // anonymous namespace


void dump_field(const char*  fname         = "fieldmap.txt",
                double       originZ_lab   = 31500.0,
                double       magnetHalfZ   = 1750.0)
{
    FieldMap fmap;
    if (!fmap.load(fname, originZ_lab)) {
        std::cerr << "Could not load the field map.\n";
        return;
    }

    // ── Bx profile along the beam axis ────────────────────────────────────
    std::cout << "=== On-axis field profile: (x=0, y=0, z) ===\n";
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  z_lab [mm]     z_map [mm]     Bx [T]      By [T]      Bz [T]\n";
    for (double z_lab = originZ_lab - 4000; z_lab <= originZ_lab + 4000; z_lab += 200) {
        std::cout << "  " << std::setw(10) << z_lab
                  << "    " << std::setw(10) << z_lab - originZ_lab
                  << "    " << std::setw(10) << fmap.Bx(0, 0, z_lab)
                  << "  " << std::setw(10) << fmap.By(0, 0, z_lab)
                  << "  " << std::setw(10) << fmap.Bz(0, 0, z_lab) << "\n";
    }

    // ── Integrated Bx over several z-windows ──────────────────────────────
    std::cout << "\n=== Straight-line ∫B_x dL along (x=0, y=0) ===\n";
    std::cout << "(These are in Tesla·metre. Compare against |p·θ/0.3|.)\n\n";

    const double zEntry = originZ_lab - magnetHalfZ;
    const double zExit  = originZ_lab + magnetHalfZ;
    std::cout << "  Simulation MagnetVolume window (halfZ = " << magnetHalfZ << " mm):\n"
              << "    z in [" << zEntry << ", " << zExit << "] mm : "
              << intBx_straight(fmap, zEntry, zExit) << " T·m\n\n";

    std::cout << "  Progressively wider windows (to see fringe contributions):\n";
    for (double half : {1000.0, 1500.0, 1750.0, 2000.0, 2500.0, 3000.0, 4000.0, 5000.0}) {
        std::cout << "    halfZ = " << std::setw(5) << half
                  << " mm,  z in [" << std::setw(6) << originZ_lab - half
                  << ", " << std::setw(6) << originZ_lab + half << "] : "
                  << intBx_straight(fmap, originZ_lab - half, originZ_lab + half)
                  << " T·m\n";
    }

    // ── Peak |Bx| on axis (fine scan) ─────────────────────────────────────
    double maxAbsBx = 0, maxAt = 0;
    for (double z = originZ_lab - 5000; z <= originZ_lab + 5000; z += 10) {
        const double b = std::abs(fmap.Bx(0, 0, z));
        if (b > maxAbsBx) { maxAbsBx = b; maxAt = z; }
    }
    std::cout << "\n=== Peak on-axis |Bx| = " << maxAbsBx
              << " T at lab z = " << maxAt << " mm"
              << " (z_map = " << maxAt - originZ_lab << " mm) ===\n";

    std::cout << std::defaultfloat;
}
