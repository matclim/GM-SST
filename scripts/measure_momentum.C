// measure_momentum.C  (v2)
// =============================================================================
// Reconstruct muon momentum from the 4-station straw-tracker output.
//
// Method (same as v1, with per-event field integration):
//   1. Group hits by station and split into upstream (0,1) / downstream (2,3).
//   2. Linear fit y(z) and x(z) on each half → slopes and intercepts.
//   3. Kink angle  θ = atan(m_y_down) − atan(m_y_up)   (bending is y-z plane).
//   4. Integrated field ∫B·dL:
//        - If a field-map file is provided: integrate B_x along a Hermite-
//          interpolated trajectory running from (x_entry, y_entry, m_up)
//          at z = 29 750 mm to (x_exit, y_exit, m_dn) at z = 33 250 mm.
//          The Hermite cubic is the unique one matching both positions AND
//          both slopes at the endpoints, giving a smooth curved path through
//          the magnet that collapses to a straight line when there's no bend.
//        - Otherwise fall back to a user-supplied constant (default 0.525 T·m,
//          which corresponds to the uniform −0.15 T × 3.5 m fallback).
//   5. Momentum  p [GeV/c] ≈ 0.3 · ∫B·dL [T·m] / θ [rad].
//
// Usage:
//   root -l -b -q 'measure_momentum.C("hits.root", 0.1, "fieldmap.txt", 0.525, 10.0)'
//
// Arguments:
//   infile              Input ROOT file from run_StrawTracker
//   sigma_xy_mm         Per-hit Gaussian smearing on (x, y) in mm
//   fieldMapFile        Optional field-map file. Empty = use fallback constant.
//   fallbackBdL_Tm      ∫B·dL used when no field map is given (T·m)
//   truthP_GeV          Known beam momentum (for resolution plots). 0 disables.
//   outfile             Output ROOT file with diagnostic histograms.
//   magnetZ_lab_mm      Lab-frame z of the magnet centre (31500 default).
//   magnetHalfZ_mm      Half-length of the MagnetVolume in z (1750 default).
//   field_ripple_rel    Per-event multiplicative rescaling of ∫B·dL, drawn
//                       from N(0, σ). Models run-to-run reproducibility /
//                       ripple (e.g. 1e-3 for the SHiP magnet). Default 0.
//   field_map_noise_T   One-off Gaussian noise added to every stored grid
//                       value at load time (Tesla). Models the measurement
//                       precision of the field map (e.g. 1e-4 for ±0.1 mT).
//                       Default 0.
//
// Coordinate conventions:
//   - Positions in the tree are in Geant4's *world* frame, which is centred
//     at lab z = 31000 mm. The macro converts to lab frame (adds 31000) so
//     everything — including the field-map z origin at lab z = 31500 — is
//     expressed consistently.
//   - Magnet volume spans lab z ∈ [29 750, 33 250] mm (halfZ = 1750 mm).
//     Edit `z_entry_lab` / `z_exit_lab` below if your MagneticFieldRegion
//     uses different dimensions.
// =============================================================================

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TRandom3.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

// =============================================================================
// Lightweight 3-D field-map reader (standalone, no Geant4 dependency).
// Mirrors the ShipFieldMap class used in the simulation.
// =============================================================================
class FieldMap {
public:
    bool load(const std::string& fname, double originZ_lab_mm = 31500.0) {
        m_originZ = originZ_lab_mm;

        std::ifstream in(fname);
        if (!in) {
            std::cerr << "[FieldMap] Cannot open '" << fname << "'\n";
            return false;
        }

        std::vector<std::array<double, 6>> rows;
        rows.reserve(1u << 16);
        std::string line;
        while (std::getline(in, line)) {
            auto p = line.find_first_not_of(" \t\r\n");
            if (p == std::string::npos || line[p] == '#') continue;
            std::istringstream iss(line);
            double x, y, z, bx, by, bz;
            if (!(iss >> x >> y >> z >> bx >> by >> bz)) continue;
            rows.push_back({x, y, z, bx, by, bz});
        }
        if (rows.empty()) {
            std::cerr << "[FieldMap] No data rows in '" << fname << "'\n";
            return false;
        }

        // Infer grid extent from the set of unique coordinates.
        std::set<double> xs, ys, zs;
        for (const auto& r : rows) { xs.insert(r[0]); ys.insert(r[1]); zs.insert(r[2]); }
        m_nX = static_cast<int>(xs.size());
        m_nY = static_cast<int>(ys.size());
        m_nZ = static_cast<int>(zs.size());
        m_xMin = *xs.begin();   m_xMax = *xs.rbegin();
        m_yMin = *ys.begin();   m_yMax = *ys.rbegin();
        m_zMin = *zs.begin();   m_zMax = *zs.rbegin();
        m_dx = (m_nX > 1) ? (m_xMax - m_xMin) / (m_nX - 1) : 0.0;
        m_dy = (m_nY > 1) ? (m_yMax - m_yMin) / (m_nY - 1) : 0.0;
        m_dz = (m_nZ > 1) ? (m_zMax - m_zMin) / (m_nZ - 1) : 0.0;

        const std::size_t N = static_cast<std::size_t>(m_nX) * m_nY * m_nZ;
        m_Bx.assign(N, 0.0f);
        m_By.assign(N, 0.0f);
        m_Bz.assign(N, 0.0f);
        for (const auto& r : rows) {
            const int ix = (m_dx > 0) ? static_cast<int>(std::lround((r[0] - m_xMin) / m_dx)) : 0;
            const int iy = (m_dy > 0) ? static_cast<int>(std::lround((r[1] - m_yMin) / m_dy)) : 0;
            const int iz = (m_dz > 0) ? static_cast<int>(std::lround((r[2] - m_zMin) / m_dz)) : 0;
            if (ix < 0 || ix >= m_nX) continue;
            if (iy < 0 || iy >= m_nY) continue;
            if (iz < 0 || iz >= m_nZ) continue;
            const std::size_t i = idx(ix, iy, iz);
            m_Bx[i] = static_cast<float>(r[3]);
            m_By[i] = static_cast<float>(r[4]);
            m_Bz[i] = static_cast<float>(r[5]);
        }

        m_valid = true;
        std::cout << "[FieldMap] Loaded '" << fname << "':\n"
                  << "           grid  " << m_nX << " x " << m_nY << " x " << m_nZ << "\n"
                  << "           X in [" << m_xMin << ", " << m_xMax << "] mm  step " << m_dx << "\n"
                  << "           Y in [" << m_yMin << ", " << m_yMax << "] mm  step " << m_dy << "\n"
                  << "           Z in [" << m_zMin << ", " << m_zMax << "] mm  step " << m_dz << "\n"
                  << "           origin at lab z = " << m_originZ << " mm\n";
        return true;
    }

    bool valid() const { return m_valid; }

    // Add independent Gaussian noise of σ = sigma_T [T] to every grid value.
    // Simulates the measurement precision of the field map itself (e.g.
    // ±0.1 mT for the SHiP V21 map). Applied once per run — all events see
    // the same perturbed map, which matches the physical situation: one
    // measurement campaign produces one realisation of measurement errors.
    void addGridNoise(double sigma_T, unsigned seed = 12345) {
        if (sigma_T <= 0.0 || !m_valid) return;
        TRandom3 rng(seed);
        for (float& v : m_Bx) v += static_cast<float>(rng.Gaus(0.0, sigma_T));
        for (float& v : m_By) v += static_cast<float>(rng.Gaus(0.0, sigma_T));
        for (float& v : m_Bz) v += static_cast<float>(rng.Gaus(0.0, sigma_T));
        std::cout << "[FieldMap] Added Gaussian noise sigma = " << sigma_T
                  << " T to each of the " << m_Bx.size()
                  << " grid points (seed=" << seed << ").\n";
    }

    // Query components at lab-frame position (x, y, z) [mm]. Returns Tesla.
    // Returns 0 if the point is outside the grid.
    double Bx(double x_lab, double y_lab, double z_lab) const {
        return interp(x_lab, y_lab, z_lab - m_originZ, m_Bx);
    }
    double By(double x_lab, double y_lab, double z_lab) const {
        return interp(x_lab, y_lab, z_lab - m_originZ, m_By);
    }
    double Bz(double x_lab, double y_lab, double z_lab) const {
        return interp(x_lab, y_lab, z_lab - m_originZ, m_Bz);
    }

private:
    std::size_t idx(int ix, int iy, int iz) const {
        return (static_cast<std::size_t>(iz) * m_nY + iy) * m_nX + ix;
    }

    double interp(double x, double y, double z, const std::vector<float>& B) const {
        if (!m_valid) return 0.0;
        if (x < m_xMin || x > m_xMax) return 0.0;
        if (y < m_yMin || y > m_yMax) return 0.0;
        if (z < m_zMin || z > m_zMax) return 0.0;

        const double fx = (m_dx > 0) ? (x - m_xMin) / m_dx : 0.0;
        const double fy = (m_dy > 0) ? (y - m_yMin) / m_dy : 0.0;
        const double fz = (m_dz > 0) ? (z - m_zMin) / m_dz : 0.0;

        int ix = static_cast<int>(fx); if (ix >= m_nX - 1) ix = m_nX - 2; if (ix < 0) ix = 0;
        int iy = static_cast<int>(fy); if (iy >= m_nY - 1) iy = m_nY - 2; if (iy < 0) iy = 0;
        int iz = static_cast<int>(fz); if (iz >= m_nZ - 1) iz = m_nZ - 2; if (iz < 0) iz = 0;

        const double tx = fx - ix;
        const double ty = fy - iy;
        const double tz = fz - iz;

        const float b000 = B[idx(ix,   iy,   iz  )];
        const float b100 = B[idx(ix+1, iy,   iz  )];
        const float b010 = B[idx(ix,   iy+1, iz  )];
        const float b110 = B[idx(ix+1, iy+1, iz  )];
        const float b001 = B[idx(ix,   iy,   iz+1)];
        const float b101 = B[idx(ix+1, iy,   iz+1)];
        const float b011 = B[idx(ix,   iy+1, iz+1)];
        const float b111 = B[idx(ix+1, iy+1, iz+1)];

        const double b00 = b000 * (1 - tx) + b100 * tx;
        const double b10 = b010 * (1 - tx) + b110 * tx;
        const double b01 = b001 * (1 - tx) + b101 * tx;
        const double b11 = b011 * (1 - tx) + b111 * tx;
        const double b0  = b00  * (1 - ty) + b10  * ty;
        const double b1  = b01  * (1 - ty) + b11  * ty;
        return b0 * (1 - tz) + b1 * tz;
    }

    std::vector<float> m_Bx, m_By, m_Bz;
    int    m_nX{0}, m_nY{0}, m_nZ{0};
    double m_xMin{0}, m_yMin{0}, m_zMin{0};
    double m_xMax{0}, m_yMax{0}, m_zMax{0};
    double m_dx{0},   m_dy{0},   m_dz{0};
    double m_originZ{31500.0};
    bool   m_valid{false};
};

// Closed-form linear fit y = a + b·x. Returns (intercept, slope).
std::pair<double, double> linfit(const std::vector<double>& x,
                                 const std::vector<double>& y) {
    const std::size_t n = x.size();
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (std::size_t i = 0; i < n; ++i) {
        sx  += x[i]; sy  += y[i];
        sxx += x[i] * x[i]; sxy += x[i] * y[i];
    }
    const double denom = n * sxx - sx * sx;
    if (std::abs(denom) < 1e-30) return {0.0, 0.0};
    const double b = (n * sxy - sx * sy) / denom;
    const double a = (sy - b * sx) / n;
    return {a, b};
}

// Cubic Hermite basis — matches positions y0, y1 and slopes m0, m1 at t=0, 1.
// dz is the total span in the natural coordinate (z1 - z0).
double hermite(double t, double y0, double y1,
               double m0, double m1, double dz) {
    const double h00 = (2 * t - 3) * t * t + 1;
    const double h10 = ((t - 2) * t + 1) * t;
    const double h01 = (3 - 2 * t) * t * t;
    const double h11 = (t - 1) * t * t;
    return h00 * y0 + h10 * dz * m0 + h01 * y1 + h11 * dz * m1;
}

// ∫B_x dL along a Hermite-interpolated trajectory from z_from to z_to (lab mm).
// Returns the integral in Tesla·metres.
double integrate_BxdL(const FieldMap& fmap,
                      double z_from, double z_to,
                      double y0, double y1, double my0, double my1,
                      double x0, double x1, double mx0, double mx1,
                      int nSteps = 500) {
    if (!fmap.valid())       return 0.0;
    const double L = z_to - z_from;
    if (L <= 0.0)            return 0.0;

    const double dz = L / nSteps;
    double sum_Bx = 0.0;
    for (int i = 0; i <= nSteps; ++i) {
        const double t  = static_cast<double>(i) / nSteps;
        const double zz = z_from + i * dz;
        const double xx = hermite(t, x0, x1, mx0, mx1, L);
        const double yy = hermite(t, y0, y1, my0, my1, L);
        const double bx = fmap.Bx(xx, yy, zz);
        const double w  = (i == 0 || i == nSteps) ? 0.5 : 1.0;  // trapezoidal
        sum_Bx += w * bx;
    }
    // sum_Bx * dz has units of Tesla·mm → convert to Tesla·m.
    return sum_Bx * dz * 1e-3;
}

} // anonymous namespace


// =============================================================================
// Main entry point
// =============================================================================
void measure_momentum(
    const char* infile             = "StrawTracker_hits_t0.root",
    double      sigma_xy_mm        = 0.0,
    const char* fieldMapFile       = "",
    double      fallbackBdL_Tm     = 0.525,
    double      truthP_GeV         = 10.0,
    const char* outfile            = "momentum.root",
    double      magnetZ_lab_mm     = 31500.0,
    double      magnetHalfZ_mm     = 1750.0,
    double      field_ripple_rel   = 0.0,
    double      field_map_noise_T  = 0.0
) {
    // ── Geometry — must match the simulation ─────────────────────────────
    // magnetZ_lab_mm and magnetHalfZ_mm are arguments so you can match a
    // modified MagneticFieldRegion without editing this file. Defaults
    // correspond to the stock build (halfZ = 1750 mm, centred at lab
    // z = 31500 mm).
    constexpr double worldZOrigin_mm = 31000.0;
    const double z_entry_lab = magnetZ_lab_mm - magnetHalfZ_mm;
    const double z_exit_lab  = magnetZ_lab_mm + magnetHalfZ_mm;

    // ── Optional field map ───────────────────────────────────────────────
    FieldMap fmap;
    const bool useMap = (std::string(fieldMapFile).size() > 0)
                        && fmap.load(fieldMapFile, magnetZ_lab_mm);
    if (useMap && field_map_noise_T > 0.0) {
        fmap.addGridNoise(field_map_noise_T);
    }

    // ── Open input tree ──────────────────────────────────────────────────
    TFile* fin = TFile::Open(infile, "READ");
    if (!fin || fin->IsZombie()) {
        std::cerr << "[measure_momentum] Cannot open " << infile << "\n";
        return;
    }
    TTree* tree = (TTree*)fin->Get("Events");
    if (!tree) {
        std::cerr << "[measure_momentum] No 'Events' tree in " << infile << "\n";
        fin->Close(); return;
    }

    std::vector<int>*    trackID   = nullptr;
    std::vector<int>*    stationID = nullptr;
    std::vector<double>* hx = nullptr;
    std::vector<double>* hy = nullptr;
    std::vector<double>* hz = nullptr;
    tree->SetBranchAddress("trackID",   &trackID);
    tree->SetBranchAddress("stationID", &stationID);
    tree->SetBranchAddress("x",         &hx);
    tree->SetBranchAddress("y",         &hy);
    tree->SetBranchAddress("z",         &hz);

    // ── Output histograms ────────────────────────────────────────────────
    TFile* fout = TFile::Open(outfile, "RECREATE");

    const double pMax = (truthP_GeV > 0.0) ? 3.0 * truthP_GeV : 50.0;

    TH1F* h_p        = new TH1F("h_p",
        "Reconstructed momentum;p_{reco} [GeV/c];Events", 200, 0.0, pMax);
    TH1F* h_theta_y  = new TH1F("h_theta_y",
        "Bending kink;#theta_{y} [mrad];Events", 200, -50.0, 50.0);
    TH1F* h_theta_x  = new TH1F("h_theta_x",
        "Non-bending kink (resolution check);#theta_{x} [mrad];Events",
        200, -10.0, 10.0);

    TH1F* h_BdL = nullptr;
    if (useMap) {
        h_BdL = new TH1F("h_BdL",
            "Per-event #int B_{x} dL;#int B_{x} dL [T#upointm];Events",
            200, -2.0, 2.0);
    }

    TH1F* h_res = nullptr;
    TH2F* h_p_vs_truth = nullptr;
    if (truthP_GeV > 0.0) {
        h_res = new TH1F("h_res",
            "Momentum resolution;(p_{reco} - p_{true}) / p_{true};Events",
            200, -1.0, 1.0);
        h_p_vs_truth = new TH2F("h_p_vs_truth",
            ";p_{true} [GeV/c];p_{reco} [GeV/c]",
            50, 0.0, 1.2 * truthP_GeV,
            200, 0.0, pMax);
    }

    // ── Event loop ───────────────────────────────────────────────────────
    TRandom3 rng(42);
    const Long64_t nEv = tree->GetEntries();

    std::cout << "[measure_momentum] Reading " << nEv << " events from " << infile << "\n"
              << "[measure_momentum] sigma_xy          = " << sigma_xy_mm << " mm\n"
              << "[measure_momentum] field ripple (rel) = " << field_ripple_rel << "\n"
              << "[measure_momentum] field map noise    = " << field_map_noise_T << " T\n";
    if (useMap) {
        std::cout << "[measure_momentum] Field map      = '" << fieldMapFile << "'\n"
                  << "[measure_momentum] Integration    : lab z in ["
                  << z_entry_lab << ", " << z_exit_lab << "] mm\n";
    } else {
        std::cout << "[measure_momentum] Fallback #intB.dL = "
                  << fallbackBdL_Tm << " T.m\n";
    }
    std::cout << "[measure_momentum] Truth p        = " << truthP_GeV << " GeV/c\n";

    long   nAcc = 0, nRejHit = 0, nRejTheta = 0;
    long   totalSecondaryHits = 0;
    double sumAbsBdL = 0.0;

    for (Long64_t iEv = 0; iEv < nEv; ++iEv) {
        tree->GetEntry(iEv);

        std::vector<double> zs_up,   xs_up,   ys_up;
        std::vector<double> zs_down, xs_down, ys_down;

        long nSecondaryHits = 0;  // diagnostic, aggregated across the run
        for (std::size_t k = 0; k < stationID->size(); ++k) {
            // Keep only hits from the primary muon. Delta rays and other
            // secondaries (trackID > 1) get excluded from the fit — their
            // low-momentum curls can otherwise hijack an upstream or
            // downstream slope and spoof an anomalously large kink,
            // producing the long low-p tail in h_p.
            if ((*trackID)[k] != 1) { ++nSecondaryHits; continue; }

            const int s = (*stationID)[k];
            double xv = (*hx)[k];
            double yv = (*hy)[k];
            // Tree is in Geant4 world frame (centred at lab z=31000); convert.
            const double zv = (*hz)[k] + worldZOrigin_mm;

            if (sigma_xy_mm > 0.0) {
                xv += rng.Gaus(0.0, sigma_xy_mm);
                yv += rng.Gaus(0.0, sigma_xy_mm);
            }

            if (s == 0 || s == 1) {
                zs_up.push_back(zv);
                xs_up.push_back(xv);
                ys_up.push_back(yv);
            } else if (s == 2 || s == 3) {
                zs_down.push_back(zv);
                xs_down.push_back(xv);
                ys_down.push_back(yv);
            }
        }
        totalSecondaryHits += nSecondaryHits;

        if (zs_up.size() < 2 || zs_down.size() < 2) {
            ++nRejHit; continue;
        }

        auto [a_y_up, m_y_up]   = linfit(zs_up,   ys_up);
        auto [a_y_dn, m_y_dn]   = linfit(zs_down, ys_down);
        auto [a_x_up, m_x_up]   = linfit(zs_up,   xs_up);
        auto [a_x_dn, m_x_dn]   = linfit(zs_down, xs_down);

        const double theta_y = std::atan(m_y_dn) - std::atan(m_y_up);
        const double theta_x = std::atan(m_x_dn) - std::atan(m_x_up);

        h_theta_y->Fill(theta_y * 1000.0);
        h_theta_x->Fill(theta_x * 1000.0);

        if (std::abs(theta_y) < 1e-7) { ++nRejTheta; continue; }

        // ── Integrated field ────────────────────────────────────────────
        double BdL_Tm = fallbackBdL_Tm;
        if (useMap) {
            // Extrapolate the upstream line to z_entry and the downstream
            // line to z_exit. These are the endpoint positions and slopes
            // that define the Hermite trajectory through the magnet.
            const double y_entry = a_y_up + m_y_up * z_entry_lab;
            const double y_exit  = a_y_dn + m_y_dn * z_exit_lab;
            const double x_entry = a_x_up + m_x_up * z_entry_lab;
            const double x_exit  = a_x_dn + m_x_dn * z_exit_lab;

            BdL_Tm = integrate_BxdL(fmap,
                z_entry_lab, z_exit_lab,
                y_entry, y_exit, m_y_up, m_y_dn,
                x_entry, x_exit, m_x_up, m_x_dn);

            h_BdL->Fill(BdL_Tm);
            sumAbsBdL += std::abs(BdL_Tm);
        }

        // Per-event ripple / reproducibility: a multiplicative rescaling of
        // the whole integrated field drawn from N(1, field_ripple_rel). This
        // applies to both the field-map and the fallback-constant cases.
        if (field_ripple_rel > 0.0) {
            BdL_Tm *= (1.0 + rng.Gaus(0.0, field_ripple_rel));
        }

        if (std::abs(BdL_Tm) < 1e-9) continue;   // no field → can't measure

        const double p_GeV = 0.3 * std::abs(BdL_Tm) / std::abs(theta_y);

        h_p->Fill(p_GeV);
        if (h_res)        h_res->Fill((p_GeV - truthP_GeV) / truthP_GeV);
        if (h_p_vs_truth) h_p_vs_truth->Fill(truthP_GeV, p_GeV);
        ++nAcc;
    }

    // ── Summary ──────────────────────────────────────────────────────────
    std::cout << "[measure_momentum] Accepted            : " << nAcc << "\n"
              << "[measure_momentum] Rejected (< 2 hits) : " << nRejHit << "\n"
              << "[measure_momentum] Rejected (theta~0)  : " << nRejTheta << "\n"
              << "[measure_momentum] Secondary hits cut  : " << totalSecondaryHits << "\n";
    if (h_p->GetEntries() > 0) {
        std::cout << "[measure_momentum] <p_reco>            = "
                  << h_p->GetMean() << " GeV/c   RMS = "
                  << h_p->GetRMS() << "\n";
    }
    if (h_res && h_res->GetEntries() > 0) {
        std::cout << "[measure_momentum] <dp/p>              = "
                  << h_res->GetMean() << "   sigma = "
                  << h_res->GetRMS() << "\n";
    }
    if (useMap && nAcc > 0) {
        std::cout << "[measure_momentum] <|int B_x dL|>      = "
                  << sumAbsBdL / nAcc << " T.m\n";
    }

    fout->Write();
    fout->Close();
    fin->Close();

    std::cout << "[measure_momentum] Histograms written to " << outfile << "\n";
}
