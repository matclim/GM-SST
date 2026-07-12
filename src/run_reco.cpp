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
#include <iostream>
#include <limits>
#include <sstream>
#include <map>
#include <optional>
#include <vector>

#include <Acts/Definitions/Algebra.hpp>
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
                                     const std::vector<Acts::Vector3>& D) {
  Acts::SquareMatrix<3> A = Acts::SquareMatrix<3>::Zero();
  Acts::Vector3 b = Acts::Vector3::Zero();
  for (std::size_t i = 0; i < P.size(); ++i) {
    const Acts::Vector3 d = D[i].normalized();
    Acts::SquareMatrix<3> M = Acts::SquareMatrix<3>::Identity() - d * d.transpose();
    A += M; b += M * P[i];
  }
  if (std::fabs(A.determinant()) < 1e-9) return P.empty() ? Acts::Vector3::Zero() : P[0];
  return A.inverse() * b;
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
  const long   nMax           = std::stol(opt(argc, argv, "--n", "-1"));
  if (fieldFile.empty()) { std::cerr << "need --field <map.root>\n"; return 1; }

  // Seed momentum grid (GeV). Curvature ~ 1/p, so scan logarithmically.
  std::vector<double> pSeeds;
  { std::stringstream ss(pSeedsStr); std::string tok;
    while (std::getline(ss, tok, ',')) if (!tok.empty()) pSeeds.push_back(std::stod(tok)); }
  if (pSeeds.empty()) pSeeds = {1., 2., 5., 10., 20., 50., 100.};
  std::cout << "seed grid:"; for (double v : pSeeds) std::cout << " " << v; 
  std::cout << " GeV   chi2Max=" << chi2Max << "\n";

  Contexts ctx;
  auto field = makeShipFieldFromRootMap(fieldFile, {0.0, 0.0, originZ});
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

  auto events = readEvents(hitsFile, /*primaryOnly=*/false);
  std::cout << "read " << events.size() << " events from " << hitsFile << "\n";

  int nVtx = 0, nRejected = 0, nKept = 0; long processed = 0;
  // efficiency breakdown
  int nEvTruth2 = 0;      // events with >=2 truth pions pointing at the detector
  int nEvHits2  = 0;      // events with >=2 pions leaving hits
  int nEvFit2   = 0;      // events with >=2 fitted tracks
  int nEvVtx    = 0;      // events with a reconstructed vertex

  // Station z (world frame) and half-apertures.
  const double kStationZ[4] = {-4500.0, -2000.0, 3000.0, 4500.0};
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

    // group SIGNAL pion hits by trackID
    std::map<int, std::vector<const RawHit*>> byTrack;
    for (const auto& h : ev.hits)
      if (h.parentID == 1 && std::abs(h.pdg) == 211) byTrack[h.trackID].push_back(&h);

    // ---- TRUTH acceptance: how many decay pions POINT AT the detector? -----
    // (one entry per distinct truth track; uses the truth vertex + truth momentum)
    std::map<int, Acts::Vector3> truthMom;   // trackID -> truth p (MeV)
    for (const auto& h : ev.hits)
      if (h.parentID == 1 && std::abs(h.pdg) == 211 && !truthMom.count(h.trackID))
        truthMom[h.trackID] = Acts::Vector3(h.vpx, h.vpy, h.vpz);
    const Acts::Vector3 tvtx(ev.truthVtxX, ev.truthVtxY, ev.truthVtxZ);
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

    if (byTrack.size() < 2) continue;

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
        const auto* line = dynamic_cast<const Acts::LineSurface*>(surf);
        if (!line) continue;
        const Acts::Vector3 C = surf->center(ctx.gctx);
        const Acts::Vector3 w = line->lineDirection(ctx.gctx);
        auto loc = surf->globalToLocal(ctx.gctx, muonPCA(E,m,C,w), m);
        if (!loc.ok()) continue;
        store.push_back(Meas{gid,(*loc)[Acts::eBoundLoc0],measResMm});
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
                                        pSeedGeV * 1000.0 /*MeV*/, qTry);
          if (!seedQ) continue;
          // override the seeder's own (noisy) momentum with the grid value
          Acts::BoundVector pv0 = seedQ->parameters();
          pv0[Acts::eBoundQOverP] = qTry / (pSeedGeV * Acts::UnitConstants::GeV);
          Acts::BoundTrackParameters seedG(refSurf->getSharedPtr(), pv0,
                                           seedQ->covariance(),
                                           Acts::ParticleHypothesis::pion());

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
        t_qTrue = (ths.front()->pdg > 0) ? +1 : -1;   // pi+ = 211
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
    Acts::Vector3 seedPos = multiLineVertex(Pv, Dv);

    // ---- DOCA 1: pairwise track-to-track (vertex consistency) --------------
    double docaMax = 0.0, docaSum = 0.0; int nPair = 0;
    for (std::size_t i = 0; i < Pv.size(); ++i)
      for (std::size_t j = i + 1; j < Pv.size(); ++j) {
        const double dd = lineLineDOCA(Pv[i], Dv[i], Pv[j], Dv[j]);
        docaMax = std::max(docaMax, dd); docaSum += dd; ++nPair;
      }
    const double docaMean = (nPair > 0) ? docaSum / nPair : 0.0;

    auto vres = vertexer.fit(ctx, fitted, seedPos);
    if (!vres) continue;

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
    b_tx=ev.truthVtxX; b_ty=ev.truthVtxY; b_tz=ev.truthVtxZ;
    b_rx=b_vx-b_tx; b_ry=b_vy-b_ty; b_rz=b_vz-b_tz;
    b_docaMax=docaMax; b_docaMean=docaMean; b_ipMax=ipMax; b_ipMean=ipMean;
    b_nTruthAcc=nTruthAcc; b_nTruthHit=nTruthHit;
    b_nFitted=static_cast<int>(fitted.size());
    b_p0=fitted[0].absoluteMomentum()/Acts::UnitConstants::GeV;
    b_p1=fitted[1].absoluteMomentum()/Acts::UnitConstants::GeV;
    b_q0=fitted[0].parameters()[Acts::eBoundQOverP]>0?+1:-1;
    b_q1=fitted[1].parameters()[Acts::eBoundQOverP]>0?+1:-1;
    vtree.Fill(); ++nVtx; ++nEvVtx;
  }

  fout.cd(); vtree.Write(); ttree.Write(); fout.Close();
  std::cout << "\n=== efficiency breakdown (of " << processed << " events) ===\n"
            << "  >=2 truth pions in acceptance : " << nEvTruth2 << "\n"
            << "  >=2 pions leaving hits        : " << nEvHits2  << "\n"
            << "  >=2 tracks fitted             : " << nEvFit2   << "\n"
            << "  vertex reconstructed          : " << nEvVtx    << "\n"
            << "  tracks kept/rejected          : " << nKept << " / " << nRejected
            << "  (chi2/ndf > " << chi2Max << ")\n"
            << "wrote " << recoOut << "\n";
  return 0;
}
