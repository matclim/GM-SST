// src/run_reco.cpp — multi-track straw reco + single-vertex (Billoir) finding.
//
//   run_reco --hits ks.root --field <map> --straws straws.root
//            [--field-origin-z 500] [--meas-res 0.15] [--reco-out reco.root]
//            [--charge-sign 1] [--n -1]
//
// Per event: select the two decay pions (parentID==1 && |pdg|==211), fit each
// straw-drift track (per-track charge from the bend), then fit one common
// vertex with the Billoir fitter. Writes a per-event vertex tree with the
// fitted-vs-truth residual.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <map>
#include <optional>
#include <vector>

#include <Acts/Definitions/Algebra.hpp>
#include <Eigen/Eigenvalues>
#include <Acts/Definitions/TrackParametrization.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/EventData/SourceLink.hpp>
#include <Acts/Surfaces/LineSurface.hpp>
#include <Acts/Surfaces/PlaneSurface.hpp>
#include <Acts/Surfaces/RectangleBounds.hpp>
#include <Acts/Surfaces/Surface.hpp>

#include <TFile.h>
#include <TTree.h>

#include "ShipRecoContext.hpp"
#include "ShipEDM.hpp"
#include "StrawDrift.h"
#include "ActsShipField.hpp"
#include "ShipHitReader.hpp"
#include "ShipSeeder.hpp"
#include "ShipStrawGeometry.hpp"
#include "ShipTrackFitter.hpp"
#include "ShipVertexer.hpp"

using namespace shipreco;

static std::string opt(int argc, char** argv, const std::string& key,
                       const std::string& def) {
  for (int i = 1; i + 1 < argc; ++i) if (key == argv[i]) return argv[i + 1];
  return def;
}
static Acts::Vector3 muonPCA(const Acts::Vector3& E, const Acts::Vector3& m,
                             const Acts::Vector3& C, const Acts::Vector3& w) {
  const double b = m.dot(w), denom = 1.0 - b * b;
  if (std::fabs(denom) < 1e-9) return E;
  const Acts::Vector3 r = E - C;
  return E + ((b * w.dot(r) - m.dot(r)) / denom) * m;
}

// 3D closest-approach midpoint of two lines (P0,d0) and (P1,d1).
static Acts::Vector3 twoLineMidpoint(const Acts::Vector3& P0, const Acts::Vector3& d0,
                                     const Acts::Vector3& P1, const Acts::Vector3& d1) {
  const double a = d0.dot(d0), b = d0.dot(d1), c = d1.dot(d1);
  const Acts::Vector3 r = P0 - P1;
  const double d = d0.dot(r), e = d1.dot(r);
  const double denom = a * c - b * b;
  double t0 = 0, t1 = 0;
  if (std::fabs(denom) > 1e-9) { t0 = (b * e - c * d) / denom; t1 = (a * e - b * d) / denom; }
  return 0.5 * ((P0 + t0 * d0) + (P1 + t1 * d1));
}

// Closest distance between two lines (P0,d0),(P1,d1).
static double lineLineDOCA(const Acts::Vector3& P0, const Acts::Vector3& d0,
                           const Acts::Vector3& P1, const Acts::Vector3& d1) {
  const Acts::Vector3 n = d0.cross(d1);
  const double nn = n.norm();
  if (nn < 1e-9) return (P1 - P0).cross(d0).norm();   // parallel
  return std::fabs((P1 - P0).dot(n)) / nn;
}

// Least-squares point closest to N lines (minimizes sum of squared distances).
static Acts::Vector3 multiLineVertex(const std::vector<Acts::Vector3>& P,
                                     const std::vector<Acts::Vector3>& D,
                                     double* condOut = nullptr) {
  Acts::SquareMatrix<3> A = Acts::SquareMatrix<3>::Zero();
  Acts::Vector3 b = Acts::Vector3::Zero();
  for (std::size_t i = 0; i < P.size(); ++i) {
    const Acts::Vector3 d = D[i].normalized();
    Acts::SquareMatrix<3> M = Acts::SquareMatrix<3>::Identity() - d * d.transpose();
    A += M; b += M * P[i];
  }
  // Conditioning: ratio of smallest to largest eigenvalue. For forward,
  // near-parallel tracks the smallest eigenvalue (along the beam) is tiny --
  // that direction is essentially unconstrained, and the seed z is unreliable.
  Eigen::SelfAdjointEigenSolver<Acts::SquareMatrix<3>> es(A);
  const auto ev = es.eigenvalues();
  if (condOut)
    *condOut = (ev(2) > 0) ? ev(0) / ev(2) : 0.0;   // 0 = degenerate, 1 = ideal
  if (std::fabs(A.determinant()) < 1e-9) return P.empty() ? Acts::Vector3::Zero() : P[0];
  return A.inverse() * b;
}

// Refine the seed's z by SCANNING along the beam. The least-squares crossing is
// ill-conditioned in z for forward, near-parallel tracks (seedCond ~ 1e-5), so
// its z is unreliable and Billoir -- which returns the seed in z -- inherits the
// error. But the TRANSVERSE positions are well measured, so at each candidate z
// we extrapolate every track (straight line: p is huge, curvature over the decay
// volume is negligible) and take the RMS transverse spread of the crossings. The
// z that minimises that spread is where the tracks actually converge. A
// well-conditioned 1-D minimisation replacing a near-singular 3-D inversion.
static double scanVertexZ(const std::vector<Acts::Vector3>& P,
                          const std::vector<Acts::Vector3>& D,
                          double zLo, double zHi) {
  auto spreadAt = [&](double z) -> double {
    double sx=0, sy=0, sxx=0, syy=0; int n=0;
    for (std::size_t i = 0; i < P.size(); ++i) {
      if (std::fabs(D[i].z()) < 1e-9) continue;
      const double t = (z - P[i].z()) / D[i].z();
      const double x = P[i].x() + t * D[i].x();
      const double y = P[i].y() + t * D[i].y();
      sx+=x; sy+=y; sxx+=x*x; syy+=y*y; ++n;
    }
    if (n < 2) return 1e30;
    return std::sqrt(std::max((sxx/n - (sx/n)*(sx/n)) + (syy/n - (sy/n)*(sy/n)), 0.0));
  };
  auto minOn = [&](double lo, double hi, double step) -> double {
    double bz=lo, bs=1e30;
    for (double z=lo; z<=hi; z+=step) { double sp=spreadAt(z); if (sp<bs){bs=sp;bz=z;} }
    return bz;
  };
  const double z1 = minOn(zLo, zHi, 1000.0);          // 1 m
  const double z2 = minOn(z1-1000.0, z1+1000.0, 50.0); // 5 cm
  return       minOn(z2-50.0, z2+50.0, 2.0);           // 2 mm
}

int main(int argc, char** argv) {
  const std::string hitsFile  = opt(argc, argv, "--hits", "ks.root");
  const std::string fieldFile = opt(argc, argv, "--field", "");
  const std::string strawFile = opt(argc, argv, "--straws", "straws.root");
  const std::string recoOut   = opt(argc, argv, "--reco-out", "reco.root");
  const double originZ        = std::stod(opt(argc, argv, "--field-origin-z", "500"));
  const double measResMm      = std::stod(opt(argc, argv, "--meas-res", "0.15"));
  const double chargeSign     = std::stod(opt(argc, argv, "--charge-sign", "1"));
  const std::string pSeedsStr = opt(argc, argv, "--p-seeds", "1,2,5,10,20,50,100");
  const double chi2Max        = std::stod(opt(argc, argv, "--chi2-max", "50"));
  // Track selection:
  //   "signal"  = K_S decay pions   (parentID==1 && |pdg|==211)
  //   "llp"     = LLP daughters     (parentID==0 && charged) -- they are fired
  //               as PRIMARIES from the LLP decay vertex, so each has its own
  //               trackID and parentID 0. Works for any number of daughters.
  //   "primary" = the single gun particle (trackID==1)
  //   "all"     = every track
  const std::string select    = opt(argc, argv, "--select", "signal");
  // Print the first N failing propagations in full (0 = tally only). The tally
  // below always runs: it groups ACTS's own error messages, so a new failure
  // mode shows up as a named bucket instead of a silent drop.
  const int verboseFail       = std::stoi(opt(argc, argv, "--verbose-fail", "0"));
  // OUTPUT FRAME. Internally everything is the WORLD frame (= SHiP - 60000 mm),
  // because that is what the straw table and the hits are in. But the generator
  // file, the experiment, and every conversation about this detector use SHiP
  // coordinates -- so the output trees are written in SHIP, and only positions
  // are shifted (residuals and covariances are frame-independent).
  const double shipZOrigin    = std::stod(opt(argc, argv, "--ship-z-origin", "60000"));
  // Target position for the reconstructed-LLP impact parameter, in SHiP z (mm).
  // Default: SHiP origin (0,0,0). The parent line (vertex + summed daughter
  // momentum) is tested for closest approach to this point -- a displaced-decay
  // discriminant: does the reconstructed decay point back to the target?
  const double ipOriginZ      = std::stod(opt(argc, argv, "--ip-origin-z", "0"));
  // Effective field integral [T*m]. Re-derive whenever the field map changes:
  //   fire a straight mu-, measure its deflection d at station 3,
  //   theta = d / (z_ST3 - z_magnet),  kBL = p * theta / 0.3.
  const double kBL            = std::stod(opt(argc, argv, "--kbl", "0.49"));
  // Dump the RECO's field along the z axis and exit. Compare against the map:
  //   Data->Draw("Bx:z","abs(x)<50 && abs(y)<50")   (z there is in CM)
  const bool dumpField        = [&]{ for (int i=1;i<argc;++i)
        if (std::string(argv[i])=="--dump-field") return true; return false; }();
  // DIAGNOSTIC ONLY (not detector-level): take the drift SIGN from the truth
  // momentum direction instead of the entry->exit chord. Tests whether the
  // chord-derived left/right assignment is what breaks large-angle tracks.
  const bool truthDirMeas     = [&]{ for (int i=1;i<argc;++i)
        if (std::string(argv[i])=="--truth-dir-meas") return true; return false; }();
  const long   nMax           = std::stol(opt(argc, argv, "--n", "-1"));
  if (fieldFile.empty()) { std::cerr << "need --field <map.root>\n"; return 1; }

  // Seed momentum grid (GeV). Curvature ~ 1/p, so scan logarithmically.
  std::vector<double> pSeeds;
  { std::stringstream ss(pSeedsStr); std::string tok;
    while (std::getline(ss, tok, ',')) if (!tok.empty()) pSeeds.push_back(std::stod(tok)); }
  if (pSeeds.empty()) pSeeds = {1., 2., 5., 10., 20., 50., 100.};
  std::cout << "seed grid:"; for (double v : pSeeds) std::cout << " " << v; 
  std::cout << " GeV   chi2Max=" << chi2Max
            << "   select=" << select
            << (truthDirMeas ? "   [DIAG: truth-dir drift sign]" : "")
            << "\n"
            << "output frame: SHiP  (world z + " << shipZOrigin << " mm)\n";

  Contexts ctx;
  auto field = makeShipFieldFromRootMap(fieldFile, {0.0, 0.0, originZ});

  // ---- --dump-field: what does the RECO actually see? ----------------------
  if (dumpField) {
    Acts::MagneticFieldContext mctx{};
    auto cache = field->makeCache(mctx);
    std::cout << "\n# reco field along the z axis (x=y=0), origin z = "
              << originZ << " mm\n";
    std::cout << "#      z[mm]        Bx[T]        By[T]        Bz[T]\n";
    double integral = 0.0;      // int Bx dz  [T*mm]
    const double zA = originZ - 9000.0, zB = originZ + 9000.0, dz = 250.0;
    double bxPrev = 0.0; bool first = true;
    for (double z = zA; z <= zB; z += dz) {
      auto r = field->getField(Acts::Vector3(0.0, 0.0, z), cache);
      const Acts::Vector3 B = r.ok() ? *r : Acts::Vector3::Zero();
      const double bx = B.x() / Acts::UnitConstants::T;
      std::printf("%12.1f %12.6f %12.6f %12.6f\n", z, bx,
                  B.y()/Acts::UnitConstants::T, B.z()/Acts::UnitConstants::T);
      if (!first) integral += 0.5 * (bx + bxPrev) * dz;   // trapezoid
      bxPrev = bx; first = false;
    }
    std::cout << "\n# integral Bx.dl = " << integral / 1000.0 << " T*m"
              << "    (muon measurement implies ~0.49 T*m)\n";
    // also probe transversely, at the magnet centre
    std::cout << "\n# transverse profile at z = " << originZ << " mm\n";
    std::cout << "#      y[mm]        Bx[T]\n";
    for (double y = -5500; y <= 5500; y += 500) {
      auto r = field->getField(Acts::Vector3(0.0, y, originZ), cache);
      const Acts::Vector3 B = r.ok() ? *r : Acts::Vector3::Zero();
      std::printf("%12.1f %12.6f\n", y, B.x()/Acts::UnitConstants::T);
    }
    return 0;
  }
  ShipStrawGeometry geom; geom.loadTable(strawFile);
  std::cout << "loaded " << geom.size() << " straw surfaces\n";

  MeasStore store;
  ShipTrackFitter fitter(field, geom, store);
  ShipVertexer vertexer(field);

  TFile fout(recoOut.c_str(), "RECREATE");
  TTree vtree("Vertices", "Reconstructed vertices");
  int    b_event=0, b_nTrk=0;
  double b_vx=0,b_vy=0,b_vz=0, b_sx=0,b_sy=0,b_sz=0;      // fitted xyz + sigma
  double b_tx=0,b_ty=0,b_tz=0, b_rx=0,b_ry=0,b_rz=0;      // truth xyz + residual
  vtree.Branch("event",&b_event); vtree.Branch("nTrk",&b_nTrk);
  vtree.Branch("vx",&b_vx); vtree.Branch("vy",&b_vy); vtree.Branch("vz",&b_vz);
  vtree.Branch("sx",&b_sx); vtree.Branch("sy",&b_sy); vtree.Branch("sz",&b_sz);
  vtree.Branch("tx",&b_tx); vtree.Branch("ty",&b_ty); vtree.Branch("tz",&b_tz);
  vtree.Branch("rx",&b_rx); vtree.Branch("ry",&b_ry); vtree.Branch("rz",&b_rz);
  double b_p0=0,b_p1=0,b_q0=0,b_q1=0;   // per-track momentum (GeV) and charge sign
  vtree.Branch("p0",&b_p0); vtree.Branch("p1",&b_p1);
  vtree.Branch("q0",&b_q0); vtree.Branch("q1",&b_q1);
  // DOCA: track-to-track (consistency) and track-to-vertex (impact parameters)
  double b_docaMax=0, b_docaMean=0, b_ipMax=0, b_ipMean=0;
  int    b_nTruthAcc=0, b_nTruthHit=0, b_nFitted=0;
  vtree.Branch("docaMax",&b_docaMax);   // max pairwise track-track DOCA [mm]
  vtree.Branch("docaMean",&b_docaMean);
  vtree.Branch("ipMax",&b_ipMax);       // max track impact param wrt vertex [mm]
  vtree.Branch("ipMean",&b_ipMean);
  vtree.Branch("nTruthAcc",&b_nTruthAcc);  // truth tracks pointing INTO the detector
  vtree.Branch("nTruthHit",&b_nTruthHit);  // truth tracks that actually left hits
  vtree.Branch("nFitted",&b_nFitted);      // tracks successfully fitted
  // Vertex SEED (multiLineVertex: least-squares point closest to all track lines).
  // For forward, near-parallel tracks this is ill-conditioned along z, so a small
  // direction error can throw it by metres -- which would place the perigee
  // surfaces in the wrong spot and break Billoir's linearisation.
  double b_sdx=0, b_sdy=0, b_sdz=0;   // seed position
  double b_seedRz=0;                  // seed z residual vs truth
  double b_seedCond=0;                // conditioning of the least-squares matrix
  int    b_fitOK=0;                   // 1 = Billoir converged, 0 = failed
  vtree.Branch("sdx",&b_sdx); vtree.Branch("sdy",&b_sdy); vtree.Branch("sdz",&b_sdz);
  vtree.Branch("seedRz",&b_seedRz);
  vtree.Branch("seedCond",&b_seedCond);
  vtree.Branch("fitOK",&b_fitOK);   // rows are written for FAILURES too
  int b_vtxFail=0, b_nProp=0;
  vtree.Branch("vtxFail",&b_vtxFail);  // 0 ok, 1 <2 tracks, 2 PROPAGATION failed, 3 BILLOIR failed
  vtree.Branch("nProp",&b_nProp);      // tracks that reached the perigee surface
  // Reconstructed-LLP pointing: parent = sum of fitted daughter momenta; its
  // line from the vertex is tested against the target at (0,0,ipOriginZ) [SHiP].
  double b_ipToOrigin=0, b_ipCApZ=0, b_parentPx=0, b_parentPy=0, b_parentPz=0;
  vtree.Branch("ipToOrigin",&b_ipToOrigin);  // transverse miss to the target [mm]
  vtree.Branch("ipCApZ",&b_ipCApZ);          // SHiP z of closest approach [mm]
  vtree.Branch("parentPx",&b_parentPx);      // reconstructed parent momentum [GeV]
  vtree.Branch("parentPy",&b_parentPy);
  vtree.Branch("parentPz",&b_parentPz);
  // per-station: HOW MANY truth tracks are in acceptance at each station
  int b_nGeo[4] = {0,0,0,0}, b_nHit[4] = {0,0,0,0};
  vtree.Branch("nGeoAcc0",&b_nGeo[0]); vtree.Branch("nGeoAcc1",&b_nGeo[1]);
  vtree.Branch("nGeoAcc2",&b_nGeo[2]); vtree.Branch("nGeoAcc3",&b_nGeo[3]);
  vtree.Branch("nHitAcc0",&b_nHit[0]); vtree.Branch("nHitAcc1",&b_nHit[1]);
  vtree.Branch("nHitAcc2",&b_nHit[2]); vtree.Branch("nHitAcc3",&b_nHit[3]);

  // ---- per-track tree: fitted vs truth (angular + momentum + charge) ------
  TTree ttree("Tracks", "Fitted tracks vs truth");
  int    t_event=0, t_pdg=0, t_nMeas=0;
  double t_pFit=0,t_pTrue=0,t_pRes=0;
  double t_phiFit=0,t_phiTrue=0,t_dphi=0;
  double t_thetaFit=0,t_thetaTrue=0,t_dtheta=0,t_dangle=0;
  double t_qFit=0,t_qTrue=0,t_chi2=0;
  ttree.Branch("event",&t_event);   ttree.Branch("pdg",&t_pdg);
  ttree.Branch("pFit",&t_pFit);     ttree.Branch("pTrue",&t_pTrue);
  ttree.Branch("pRes",&t_pRes);
  ttree.Branch("phiFit",&t_phiFit); ttree.Branch("phiTrue",&t_phiTrue);
  ttree.Branch("dphi",&t_dphi);
  ttree.Branch("thetaFit",&t_thetaFit); ttree.Branch("thetaTrue",&t_thetaTrue);
  ttree.Branch("dtheta",&t_dtheta); ttree.Branch("dangle",&t_dangle);
  ttree.Branch("qFit",&t_qFit);     ttree.Branch("qTrue",&t_qTrue);
  ttree.Branch("chi2",&t_chi2);     ttree.Branch("nMeas",&t_nMeas);
  double t_bestChi2=0;
  ttree.Branch("bestChi2",&t_bestChi2);
  int t_inAcc=0;
  ttree.Branch("inAcc",&t_inAcc);   // truth: crosses ALL 4 stations (geometric)
  // per-station acceptance: geometric (straight-line truth) vs actual hits
  int t_geoAcc[4] = {0,0,0,0};
  int t_hitAcc[4] = {0,0,0,0};
  int t_nGeoAcc=0, t_nHitAcc=0;
  ttree.Branch("geoAcc0",&t_geoAcc[0]); ttree.Branch("geoAcc1",&t_geoAcc[1]);
  ttree.Branch("geoAcc2",&t_geoAcc[2]); ttree.Branch("geoAcc3",&t_geoAcc[3]);
  ttree.Branch("hitAcc0",&t_hitAcc[0]); ttree.Branch("hitAcc1",&t_hitAcc[1]);
  ttree.Branch("hitAcc2",&t_hitAcc[2]); ttree.Branch("hitAcc3",&t_hitAcc[3]);
  ttree.Branch("nGeoAcc",&t_nGeoAcc);   // 0-4 stations geometrically crossed
  ttree.Branch("nHitAcc",&t_nHitAcc);   // 0-4 stations actually hit

  // ---- per-STRAW measurement diagnostics (truth-level; not used in the fit) --
  // driftMeas : the signed drift we hand the fit (from the entry->exit chord)
  // driftTrue : exact signed DOCA of the TRUTH trajectory to the same wire
  // diff      : driftMeas - driftTrue  -> a coherent z-dependence here is a
  //             fake curvature, i.e. exactly what would bias the momentum.
  TTree htree("Hits", "Per-straw measurement diagnostics");
  int    h_event=0, h_trk=0, h_station=0, h_layer=0, h_sub=0, h_straw=0;
  double h_z=0, h_driftMeas=0, h_driftTrue=0, h_diff=0;
  double h_theta=0, h_pTrue=0, h_chordDot=0, h_pathLen=0;
  int    h_nSteps=0;
  htree.Branch("event",&h_event);     htree.Branch("trk",&h_trk);
  htree.Branch("station",&h_station); htree.Branch("layer",&h_layer);
  htree.Branch("sub",&h_sub);         htree.Branch("straw",&h_straw);
  htree.Branch("z",&h_z);
  htree.Branch("driftMeas",&h_driftMeas);
  double h_truthLoc0=0; int h_wrongSide=0;
  htree.Branch("truthLoc0",&h_truthLoc0);   // SIGNED truth drift
  htree.Branch("wrongSide",&h_wrongSide);   // 1 if the fit chose the wrong side
  htree.Branch("driftTrue",&h_driftTrue);
  htree.Branch("diff",&h_diff);
  htree.Branch("theta",&h_theta);     htree.Branch("pTrue",&h_pTrue);
  htree.Branch("chordDot",&h_chordDot); // cos(angle) between chord and truth dir
  htree.Branch("pathLen",&h_pathLen);   // |exit-entry| inside the straw
  htree.Branch("nSteps",&h_nSteps);

  auto events = readEvents(hitsFile, /*primaryOnly=*/false);
  std::cout << "read " << events.size() << " events from " << hitsFile << "\n";

  int nVtx = 0, nRejected = 0, nKept = 0; long processed = 0;
  int nVtxFail[4] = {0,0,0,0};   // 1 = <2 tracks, 2 = propagation, 3 = Billoir
  std::map<std::string,int> propErrTally;   // ACTS error message -> count
  int nPropOK = 0, nPropBad = 0, nVerbosePrinted = 0;
  // efficiency breakdown
  int nEvTruth2 = 0;      // events with >=2 truth pions pointing at the detector
  int nEvHits2  = 0;      // events with >=2 pions leaving hits
  int nEvFit2   = 0;      // events with >=2 fitted tracks
  int nEvVtx    = 0;      // events with a reconstructed vertex

  // Station z (WORLD frame) derived from the straw table -- never hardcoded,
  // so a geometry change (e.g. the SHiP re-frame) can't silently desync.
  double kStationZ[4] = {0,0,0,0};
  {
    std::array<double,4> zsum{0,0,0,0}; std::array<int,4> zn{0,0,0,0};
    for (int st = 0; st < 4; ++st)
      for (int l = 0; l < 4; ++l)
        for (int sub = 0; sub < 2; ++sub) {
          const Acts::Surface* sf = geom.surfaceFor(st, l, sub, 150);
          if (sf) { zsum[st] += sf->center(ctx.gctx).z(); ++zn[st]; }
        }
    for (int st = 0; st < 4; ++st) kStationZ[st] = (zn[st] > 0) ? zsum[st]/zn[st] : 0.0;
    std::cout << "station z: world";
    for (double z : kStationZ) std::cout << " " << z;
    std::cout << "  |  SHiP";
    for (double z : kStationZ) std::cout << " " << z + shipZOrigin;
    std::cout << " mm\n";
  }
  const double kHalfX = 2000.0, kHalfY = 3000.0;

  // GEOMETRIC acceptance, PER STATION: straight-line swim of the truth momentum
  // from the truth vertex. Exact upstream of the magnet (field-free); for the
  // downstream stations (2,3) it ignores bending, so compare with hitAcc below:
  // the DIFFERENCE (geoAcc - hitAcc) is the bending/interaction loss.
  auto geoAccPerStation = [&](const Acts::Vector3& vtx, const Acts::Vector3& p,
                              std::array<int,4>& acc) {
    acc = {0,0,0,0};
    if (p.z() <= 0) return;
    const Acts::Vector3 d = p.normalized();
    for (int k = 0; k < 4; ++k) {
      const double t = (kStationZ[k] - vtx.z()) / d.z();
      if (t < 0) continue;
      const Acts::Vector3 at = vtx + t * d;
      acc[k] = (std::fabs(at.x()) <= kHalfX && std::fabs(at.y()) <= kHalfY) ? 1 : 0;
    }
  };
  // convenience: all four stations
  auto truthInAcceptance = [&](const Acts::Vector3& vtx, const Acts::Vector3& p) {
    std::array<int,4> a; geoAccPerStation(vtx, p, a);
    return a[0] && a[1] && a[2] && a[3];
  };
  for (const auto& ev : events) {
    if (nMax >= 0 && processed >= nMax) break;
    const int eventId = static_cast<int>(processed); ++processed;

    // group selected hits by trackID
    auto keepHit = [&](const RawHit& h) {
      if (select == "primary") return h.trackID == 1;
      if (select == "all")     return true;
      if (select == "llp") {
        // LLP daughters are primaries (parentID==0). Keep the charged, tracked
        // species; reject the interaction secondaries (parentID>0).
        const int a = std::abs(h.pdg);
        return h.parentID == 0 &&
               (a == 211 || a == 321 || a == 13 || a == 11 || a == 2212);
      }
      return h.parentID == 1 && std::abs(h.pdg) == 211;   // "signal"
    };
    std::map<int, std::vector<const RawHit*>> byTrack;
    for (const auto& h : ev.hits)
      if (keepHit(h)) byTrack[h.trackID].push_back(&h);

    // ---- TRUTH acceptance: how many decay pions POINT AT the detector? -----
    // (one entry per distinct truth track; uses the truth vertex + truth momentum)
    std::map<int, Acts::Vector3> truthMom;   // trackID -> truth p (MeV)
    for (const auto& h : ev.hits)
      if (keepHit(h) && !truthMom.count(h.trackID))
        truthMom[h.trackID] = Acts::Vector3(h.vpx, h.vpy, h.vpz);
    // Truth vertex: for "signal" it's the K_S decay point (mean over daughters);
    // for a primary/gun particle it's that track's own production vertex.
    Acts::Vector3 tvtx(ev.truthVtxX, ev.truthVtxY, ev.truthVtxZ);
    if (select == "primary" || select == "llp" || !ev.hasTruthVtx) {
      for (const auto& h : ev.hits)
        if (keepHit(h)) { tvtx = Acts::Vector3(h.vtxX, h.vtxY, h.vtxZ); break; }
    }
    int nTruthAcc = 0;
    std::array<int,4> nGeoStn{0,0,0,0}, nHitStn{0,0,0,0};
    for (const auto& [tid, tp] : truthMom) {
      std::array<int,4> ga; geoAccPerStation(tvtx, tp, ga);
      for (int k = 0; k < 4; ++k) nGeoStn[k] += ga[k];
      if (ga[0] && ga[1] && ga[2] && ga[3]) ++nTruthAcc;
      // stations actually hit by this truth track
      if (byTrack.count(tid))
        for (int k = 0; k < 4; ++k)
          for (const auto* h : byTrack.at(tid))
            if (h->stationID == k) { ++nHitStn[k]; break; }
    }
    const int nTruthHit = static_cast<int>(byTrack.size());
    for (int k = 0; k < 4; ++k) { b_nGeo[k] = nGeoStn[k]; b_nHit[k] = nHitStn[k]; }
    if (nTruthAcc >= 2) ++nEvTruth2;
    if (nTruthHit >= 2) ++nEvHits2;

    if (byTrack.empty()) continue;   // vertexing needs >=2, enforced later

    std::vector<Acts::BoundTrackParameters> fitted;   // must outlive vertex fit

    for (auto& [tid, ths] : byTrack) {
      // build a per-track RawEvent (seeder + DOCA reuse the same structures)
      RawEvent te; for (auto* h : ths) te.hits.push_back(*h);

      store.clear();
      std::vector<Acts::SourceLink> sls;
      std::vector<const Acts::Surface*> seq;
      std::map<uint64_t, std::vector<const RawHit*>> byStraw;
      for (auto* h : ths)
        byStraw[ShipStrawGeometry::geoId(h->stationID,h->layerID,h->subLayerID,h->strawID).value()]
            .push_back(h);
      for (auto& [key, hs] : byStraw) {
        Acts::GeometryIdentifier gid(key);
        const Acts::Surface* surf = geom.surfaceForGeoId(gid);
        if (!surf) continue;
        std::sort(hs.begin(), hs.end(),
                  [](const RawHit* a, const RawHit* b){ return a->z < b->z; });
        Acts::Vector3 E(hs.front()->xe,hs.front()->ye,hs.front()->ze);
        Acts::Vector3 X(hs.back()->xx, hs.back()->yx, hs.back()->zx);
        Acts::Vector3 m = X - E; if (m.norm()==0) continue; m.normalize();

        // DIAGNOSTIC: override the direction used to SIGN the drift distance
        // with the truth momentum direction. The chord (entry->exit over a
        // 20 mm straw) is a poor estimate of the track direction, and the sign
        // convention (measX = e_z x p_hat) rotates with it -- suspected to
        // break the left/right assignment for obliquely-crossing tracks.
        Acts::Vector3 mSign = m;
        if (truthDirMeas) {
          Acts::Vector3 tp(hs.front()->vpx, hs.front()->vpy, hs.front()->vpz);
          if (tp.norm() > 0) mSign = tp.normalized();
        }

        const auto* line = dynamic_cast<const Acts::LineSurface*>(surf);
        if (!line) continue;
        const Acts::Vector3 C = surf->center(ctx.gctx);
        const Acts::Vector3 w = line->lineDirection(ctx.gctx);
        auto loc = surf->globalToLocal(ctx.gctx, muonPCA(E,m,C,w), mSign);
        if (!loc.ok()) continue;
        // THE MEASUREMENT. A real straw gives a drift TIME; we invert it with
        // the calibrated r-t relation to an UNSIGNED radius. The side (left or
        // right of the wire) is NOT ours to supply -- the fit resolves it (see
        // ShipCalibrator). Feeding the chord-derived signed drift, as we used
        // to, hands the fit a truth-derived answer to the hardest part of the
        // problem.
        const double driftTruth = (*loc)[Acts::eBoundLoc0];   // SIGNED, truth
        // The sim writes one driftTime per STEP; the steps of a straw are
        // grouped in `hs`, and a straw has one wire and one TDC -- so take the
        // EARLIEST arrival across them, which is exactly what the TDC records.
        double tDrift = -1.0;
        for (const RawHit* rh : hs)
          if (rh->driftTime >= 0.0 && (tDrift < 0.0 || rh->driftTime < tDrift))
            tDrift = rh->driftTime;
        if (tDrift < 0.0) continue;        // straw did not fire (no clusters)

        const double radius = strawdrift::radiusFromTime(tDrift);

        const double driftMeas = radius;   // unsigned, what the reco sees
        store.push_back(Meas{gid, radius, measResMm, driftTruth});

        // ---- DIAGNOSTIC: the same drift, from the TRUTH trajectory ----------
        {
          Acts::Vector3 tp(hs.front()->vpx, hs.front()->vpy, hs.front()->vpz);
          Acts::Vector3 tv(hs.front()->vtxX, hs.front()->vtxY, hs.front()->vtxZ);
          if (tp.norm() > 0) {
            const Acts::Vector3 td = tp.normalized();
            // straight-line truth swim from the production vertex
            const Acts::Vector3 pcaT = muonPCA(tv, td, C, w);
            auto locT = surf->globalToLocal(ctx.gctx, pcaT, td);
            h_driftTrue = locT.ok() ? (*locT)[Acts::eBoundLoc0]
                                    : std::numeric_limits<double>::quiet_NaN();
            h_theta   = std::acos(std::clamp(td.z(), -1.0, 1.0));
            h_pTrue   = tp.norm() / 1000.0;
            h_chordDot= m.dot(td);
          }
          h_event = eventId; h_trk = tid;
          h_station = hs.front()->stationID; h_layer = hs.front()->layerID;
          h_sub = hs.front()->subLayerID;    h_straw = hs.front()->strawID;
          h_z = C.z();
          h_driftMeas  = driftMeas;
          h_truthLoc0  = driftTruth;
          h_diff = driftMeas - h_driftTrue;
          h_pathLen = (X - E).norm();
          h_nSteps = static_cast<int>(hs.size());
          htree.Fill();
        }
        sls.emplace_back(ShipSourceLink{gid, store.size()-1});
        seq.push_back(surf);
      }
      if (seq.size() < 6) continue;
      std::sort(seq.begin(), seq.end(), [&](const Acts::Surface* a, const Acts::Surface* b){
        return a->center(ctx.gctx).z() < b->center(ctx.gctx).z(); });

      const double zRef = seq.front()->center(ctx.gctx).z() - 50.0;
      auto refSurf = Acts::Surface::makeShared<Acts::PlaneSurface>(
          Acts::Transform3(Acts::Translation3(0.0,0.0,zRef)),
          std::make_shared<Acts::RectangleBounds>(5000.0,5000.0));

      // ---- ITERATIVE SEEDING: scan (p_seed x charge), keep best chi2/ndf ---
      // The sagitta momentum/sign estimate is noise-limited, so don't trust it.
      // Instead try a grid of seed momenta with BOTH charges, let the full
      // Kalman fit (all straw hits, real field) arbitrate, and REJECT the track
      // if no hypothesis achieves an acceptable chi2/ndf. No kinematic clamp.
      TrackContainer tracks{std::make_shared<Acts::VectorTrackContainer>(),
                            std::make_shared<Acts::VectorMultiTrajectory>()};
      double bestChi2 = std::numeric_limits<double>::max();
      bool   haveFit  = false;

      for (double pSeedGeV : pSeeds) {
        for (double qTry : {+1.0, -1.0}) {
          auto seedQ = seedFromUpstream(ctx, te, *refSurf,
                                        pSeedGeV * 1000.0 /*MeV*/, qTry,
                                        0, 1, 2, 3, kBL);
          if (!seedQ) continue;
          // override the seeder's own (noisy) momentum with the grid value
          Acts::BoundVector pv0 = seedQ->parameters();
          const double qopG = qTry / (pSeedGeV * Acts::UnitConstants::GeV);
          pv0[Acts::eBoundQOverP] = qopG;
          // Rebuild the q/p variance for THIS grid momentum. Reusing the
          // seeder's (computed from its unreliable sagitta estimate) pins the
          // fit to the seed and it cannot move off it.
          auto covG = seedQ->covariance();
          if (covG) (*covG)(Acts::eBoundQOverP, Acts::eBoundQOverP) =
                        std::pow(qopG * 1.0, 2);   // loose prior: let the fit find p
          Acts::BoundTrackParameters seedG(refSurf->getSharedPtr(), pv0,
                                           covG,
                                           (select == "primary")
                                             ? Acts::ParticleHypothesis::muon()
                                             : Acts::ParticleHypothesis::pion());

          TrackContainer tq{std::make_shared<Acts::VectorTrackContainer>(),
                            std::make_shared<Acts::VectorMultiTrajectory>()};
          if (!fitter.fit(ctx, sls, seedG, seq, tq)) continue;
          for (const auto& t : tq) {
            const int ndf = std::max(1, static_cast<int>(t.nMeasurements()) - 5);
            const double c2 = t.chi2() / ndf;
            if (std::isfinite(c2) && c2 < bestChi2) {
              bestChi2 = c2;
              tracks   = std::move(tq);
              haveFit  = true;
            }
            break;
          }
        }
      }
      // REJECT the track if no hypothesis fits acceptably.
      if (!haveFit || bestChi2 > chi2Max) { ++nRejected; continue; }
      for (const auto& trk : tracks) {
        auto par = trk.createParametersAtReference();
        if (!par.parameters().allFinite()) continue;
        if (par.covariance() && !par.covariance()->allFinite()) continue;
        const double pGeV = par.absoluteMomentum() / Acts::UnitConstants::GeV;
        if (!std::isfinite(pGeV)) continue;

        // ---- fitted vs truth for THIS track ----------------------------
        const auto pv = par.parameters();
        const double phiF = pv[Acts::eBoundPhi], thF = pv[Acts::eBoundTheta];
        Acts::Vector3 dF(std::sin(thF)*std::cos(phiF),
                         std::sin(thF)*std::sin(phiF), std::cos(thF));
        Acts::Vector3 pT(ths.front()->vpx, ths.front()->vpy, ths.front()->vpz);
        const double pTrueGeV = pT.norm() / 1000.0;   // MeV -> GeV
        Acts::Vector3 dT = (pT.norm() > 0) ? pT.normalized() : Acts::Vector3::UnitZ();
        const double phiT = std::atan2(dT.y(), dT.x());
        const double thT  = std::acos(std::clamp(dT.z(), -1.0, 1.0));

        t_event = eventId; t_pdg = ths.front()->pdg;
        t_pFit = pGeV;  t_pTrue = pTrueGeV;
        t_pRes = (pTrueGeV > 0) ? (pGeV - pTrueGeV) / pTrueGeV : 0.0;
        t_phiFit = phiF; t_phiTrue = phiT;
        t_dphi = std::remainder(phiF - phiT, 2.0 * M_PI);
        t_thetaFit = thF; t_thetaTrue = thT; t_dtheta = thF - thT;
        t_dangle = std::acos(std::clamp(dF.dot(dT), -1.0, 1.0));   // 3D opening angle
        t_qFit = (pv[Acts::eBoundQOverP] > 0) ? +1 : -1;
        // truth charge from PDG (pi+ = +211 -> +1;  mu- = 13 -> -1)
        {
          const int pd = ths.front()->pdg;
          if (std::abs(pd) == 13)       t_qTrue = (pd > 0) ? -1 : +1;   // mu-/mu+
          else                          t_qTrue = (pd > 0) ? +1 : -1;   // pi+/pi-
        }
        t_chi2 = trk.chi2(); t_nMeas = static_cast<int>(trk.nMeasurements());
        t_bestChi2 = bestChi2;
        {
          std::array<int,4> ga; geoAccPerStation(tvtx, pT, ga);
          std::array<int,4> ha{0,0,0,0};
          for (const auto* h : ths) if (h->stationID >= 0 && h->stationID < 4) ha[h->stationID] = 1;
          t_nGeoAcc = 0; t_nHitAcc = 0;
          for (int k = 0; k < 4; ++k) {
            t_geoAcc[k] = ga[k]; t_hitAcc[k] = ha[k];
            t_nGeoAcc += ga[k]; t_nHitAcc += ha[k];
          }
          t_inAcc = (t_nGeoAcc == 4) ? 1 : 0;
        }
        ttree.Fill(); ++nKept;

        if (pGeV < 0.05) continue;   // only reject degenerate fits
        fitted.push_back(std::move(par));
      }
    }

    if (fitted.size() < 2) continue;
    ++nEvFit2;

    // ---- N-track vertex seed: least-squares point closest to all tracks ----
    auto dirOf = [](const Acts::BoundTrackParameters& p) {
      const double th = p.parameters()[Acts::eBoundTheta];
      const double ph = p.parameters()[Acts::eBoundPhi];
      return Acts::Vector3(std::sin(th)*std::cos(ph),
                           std::sin(th)*std::sin(ph), std::cos(th));
    };
    std::vector<Acts::Vector3> Pv, Dv;
    for (const auto& f : fitted) { Pv.push_back(f.position(ctx.gctx)); Dv.push_back(dirOf(f)); }
    double seedCond = 0.0;
    Acts::Vector3 seedPos = multiLineVertex(Pv, Dv, &seedCond);
    // The LS z is unreliable here; keep its transverse crossing, refine z by the
    // scan. Bracket by the tracker (kStationZ in WORLD frame) and 60 m upstream.
    {
      const double zHi = kStationZ[0] - 200.0;
      const double zLo = zHi - 60000.0;
      seedPos.z() = scanVertexZ(Pv, Dv, zLo, zHi);
    }

    // ---- DOCA 1: pairwise track-to-track (vertex consistency) --------------
    double docaMax = 0.0, docaSum = 0.0; int nPair = 0;
    for (std::size_t i = 0; i < Pv.size(); ++i)
      for (std::size_t j = i + 1; j < Pv.size(); ++j) {
        const double dd = lineLineDOCA(Pv[i], Dv[i], Pv[j], Dv[j]);
        docaMax = std::max(docaMax, dd); docaSum += dd; ++nPair;
      }
    const double docaMean = (nPair > 0) ? docaSum / nPair : 0.0;

    auto vres = vertexer.fit(ctx, fitted, seedPos);

    // ---- propagation diagnostics: tally ACTS's own error messages -----------
    if (vres) {
      for (const auto& d : vres->propDiag) {
        if (d.ok) { ++nPropOK; continue; }
        ++nPropBad;
        ++propErrTally[d.error];
        if (nVerbosePrinted < verboseFail) {
          std::cout << "  [prop fail] event " << eventId
                    << (d.backward ? "  backward" : "  forward")
                    << "  dist=" << d.distToVtx << " mm"
                    << "  :: " << d.error << "\n";
          ++nVerbosePrinted;
        }
      }
    }

    b_sdx = seedPos.x(); b_sdy = seedPos.y(); b_sdz = seedPos.z();
    b_seedCond = seedCond;
    const bool ok = vres && vres->fail == shipreco::VertexFail::None;
    b_fitOK   = ok ? 1 : 0;
    b_vtxFail = vres ? static_cast<int>(vres->fail) : 1;
    b_nProp   = vres ? vres->nPropagated : 0;
    if (!ok) ++nVtxFail[b_vtxFail < 4 ? b_vtxFail : 0];

    if (!ok) {
      // failed fit: fill truth + seed, leave the fitted vertex at the seed
      // (no parent line without a fitted vertex, so the IP is undefined)
      b_ipToOrigin = -1.0; b_ipCApZ = 0.0;
      b_parentPx = b_parentPy = b_parentPz = 0.0;
      b_tx = tvtx.x(); b_ty = tvtx.y(); b_tz = tvtx.z();
      b_seedRz = b_sdz - b_tz;               // residual: frame-free
      b_vx = b_sdx; b_vy = b_sdy; b_vz = b_sdz;
      b_rx = 0; b_ry = 0; b_rz = 0;
      b_vz  += shipZOrigin;                  // -> SHiP frame
      b_tz  += shipZOrigin;
      b_sdz += shipZOrigin;
      b_nTrk = static_cast<int>(fitted.size());
      b_docaMax = docaMax; b_docaMean = docaMean;
      b_ipMax = 0; b_ipMean = 0;
      b_nTruthAcc = nTruthAcc; b_nTruthHit = nTruthHit;
      b_nFitted = static_cast<int>(fitted.size());
      b_p0 = fitted[0].absoluteMomentum()/Acts::UnitConstants::GeV;
      b_p1 = (fitted.size()>1) ? fitted[1].absoluteMomentum()/Acts::UnitConstants::GeV : 0.0;
      b_q0 = fitted[0].parameters()[Acts::eBoundQOverP]>0?+1:-1;
      b_q1 = (fitted.size()>1) ? (fitted[1].parameters()[Acts::eBoundQOverP]>0?+1:-1) : 0;
      vtree.Fill();
      continue;
    }

    // ---- DOCA 2: each track's impact parameter wrt the FITTED vertex -------
    double ipMax = 0.0, ipSum = 0.0;
    for (std::size_t i = 0; i < Pv.size(); ++i) {
      const Acts::Vector3 r = vres->position - Pv[i];
      const double ip = r.cross(Dv[i].normalized()).norm();   // perpendicular distance
      ipMax = std::max(ipMax, ip); ipSum += ip;
    }
    const double ipMean = Pv.empty() ? 0.0 : ipSum / Pv.size();

    b_event=eventId; b_nTrk=vres->nTracks;
    b_vx=vres->position.x(); b_vy=vres->position.y(); b_vz=vres->position.z();
    b_sx=std::sqrt(std::max(0.0,vres->covariance(0,0)));
    b_sy=std::sqrt(std::max(0.0,vres->covariance(1,1)));
    b_sz=std::sqrt(std::max(0.0,vres->covariance(2,2)));
    // Use the truth vertex computed for THIS selection mode (tvtx), not the
    // K_S-specific mean over parentID==1 daughters -- which is meaningless
    // for an LLP sample and silently placed the reference in the tracker.
    b_tx=tvtx.x(); b_ty=tvtx.y(); b_tz=tvtx.z();
    b_rx=b_vx-b_tx; b_ry=b_vy-b_ty; b_rz=b_vz-b_tz;   // residuals: frame-free
    b_seedRz = b_sdz - b_tz;
    // -> SHiP frame for every POSITION written out
    b_vz  += shipZOrigin;
    b_tz  += shipZOrigin;
    b_sdz += shipZOrigin;
    // ---- reconstructed-LLP impact parameter to the target -----------------
    // Parent momentum = vector sum of the fitted daughters (physical: available
    // in data). The parent trajectory is the line through the fitted vertex V
    // along p_hat; the target O is at (0,0,ipOriginZ) in SHiP = world z minus
    // shipZOrigin. Work in the WORLD frame, where V and the momenta live.
    Acts::Vector3 pSum = Acts::Vector3::Zero();
    for (const auto& f : fitted)
      pSum += (f.absoluteMomentum()/Acts::UnitConstants::GeV) * dirOf(f);
    b_parentPx = pSum.x(); b_parentPy = pSum.y(); b_parentPz = pSum.z();

    {
      const Acts::Vector3 O(0.0, 0.0, ipOriginZ - shipZOrigin);   // target, world
      const Acts::Vector3 V = vres->position;                      // vertex, world
      if (pSum.norm() > 0) {
        const Acts::Vector3 ph = pSum.normalized();
        const Acts::Vector3 d  = O - V;
        b_ipToOrigin = d.cross(ph).norm();                 // transverse miss [mm]
        const double tCA = d.dot(ph);                      // param to closest approach
        b_ipCApZ = (V.z() + tCA*ph.z()) + shipZOrigin;     // -> SHiP z
      } else {
        b_ipToOrigin = -1.0; b_ipCApZ = 0.0;
      }
    }

    b_docaMax=docaMax; b_docaMean=docaMean; b_ipMax=ipMax; b_ipMean=ipMean;
    b_nTruthAcc=nTruthAcc; b_nTruthHit=nTruthHit;
    b_nFitted=static_cast<int>(fitted.size());
    b_p0=fitted[0].absoluteMomentum()/Acts::UnitConstants::GeV;
    b_p1=fitted[1].absoluteMomentum()/Acts::UnitConstants::GeV;
    b_q0=fitted[0].parameters()[Acts::eBoundQOverP]>0?+1:-1;
    b_q1=fitted[1].parameters()[Acts::eBoundQOverP]>0?+1:-1;
    vtree.Fill(); ++nVtx; ++nEvVtx;
  }

  fout.cd(); vtree.Write(); ttree.Write(); htree.Write(); fout.Close();
  std::cout << "\n=== efficiency breakdown (of " << processed << " events) ===\n"
            << "  >=2 truth pions in acceptance : " << nEvTruth2 << "\n"
            << "  >=2 pions leaving hits        : " << nEvHits2  << "\n"
            << "  >=2 tracks fitted             : " << nEvFit2   << "\n"
            << "  vertex reconstructed          : " << nEvVtx    << "\n"
            << "  vertex fit failures           : "
            << " <2trk=" << nVtxFail[1]
            << "  propagation=" << nVtxFail[2]
            << "  Billoir=" << nVtxFail[3] << "\n"
            << "  tracks kept/rejected          : " << nKept << " / " << nRejected
            << "  (chi2/ndf > " << chi2Max << ")\n"
            << "  perigee propagation           : " << nPropOK << " ok, "
            << nPropBad << " failed\n"
            << "wrote " << recoOut << "\n";
  if (!propErrTally.empty()) {
    std::cout << "\n  propagation failure modes (ACTS):\n";
    for (const auto& [msg, n] : propErrTally)
      std::cout << "    " << std::setw(6) << n << "  " << msg << "\n";
    std::cout << "  (re-run with --verbose-fail 10 for per-track detail)\n";
  }
  return 0;
}
