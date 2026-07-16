// tests/evdisplay.C — event display in the bend plane (z,y) from REAL hit data.
//
//   root -l -b -q 'tests/evdisplay.C("build/pit_0.0_t-1.root","build/pit_0.06_t-1.root")'
//
// Selects, from each file, an event whose PRIMARY track is fully inside the
// tracker acceptance -- i.e. it left hits in all FOUR stations. That is the
// point of the display: the tilted track is a perfectly good, fully-accepted,
// 4-station, ~32-hit track, and its momentum is STILL unmeasurable, because it
// has drifted outside the +-500 mm field map by the time it reaches the magnet.
//
//   Tracker acceptance  (hits in all 4 stations)      -> the tilted track PASSES
//   Field-map coverage  (|y| < 500 mm at the magnet)  -> the tilted track FAILS
//
// The gap between those two is the whole finding.

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TBox.h>
#include <TLine.h>
#include <TLatex.h>
#include <TH2F.h>
#include <TStyle.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdio>

namespace {
// world-frame geometry (mm) = SHiP frame - kWorldZOrigin(60000)
constexpr double kStationZ[4] = {24320., 26820., 31820., 33320.};
constexpr double kMagnetZ     = 29220.;
constexpr double kMapHalfZ    = 750.;    // field map half-length in z
constexpr double kMapHalfY    = 500.;    // field map half-height in y
constexpr double kAperture    = 3000.;   // straw station half-height in y
constexpr double kGunZ        = 19320.;

struct Ev {
  std::vector<double> z, y;
  int    nStations = 0;      // how many of the 4 stations were hit
  double yMag      = 0;      // y at the magnet plane (interpolated across the gap)
  bool   inMap     = false;  // |yMag| < 500 -> the field actually acts on it
  bool   ok        = false;
};

// First event whose PRIMARY track hits ALL FOUR stations.
Ev findAcceptedEvent(const char* file, int maxScan = 200) {
  Ev out;
  TFile* f = TFile::Open(file);
  if (!f || f->IsZombie()) { printf("  cannot open %s\n", file); return out; }
  TTree* t = (TTree*)f->Get("Events");
  if (!t) { printf("  no Events tree in %s\n", file); return out; }

  std::vector<int>    *trkID = nullptr, *stID = nullptr;
  std::vector<double> *hy = nullptr, *hz = nullptr;
  t->SetBranchAddress("trackID",   &trkID);
  t->SetBranchAddress("stationID", &stID);
  t->SetBranchAddress("y", &hy);
  t->SetBranchAddress("z", &hz);

  const Long64_t n = std::min<Long64_t>(t->GetEntries(), maxScan);
  for (Long64_t i = 0; i < n; ++i) {
    t->GetEntry(i);
    bool hit[4] = {false, false, false, false};
    std::vector<double> zz, yy;
    for (size_t h = 0; h < trkID->size(); ++h) {
      if ((*trkID)[h] != 1) continue;                 // the gun particle only
      const int s = (*stID)[h];
      if (s >= 0 && s < 4) hit[s] = true;
      zz.push_back((*hz)[h]); yy.push_back((*hy)[h]);
    }
    const int nSt = hit[0] + hit[1] + hit[2] + hit[3];
    if (nSt < 4) continue;                            // NOT in acceptance -> skip

    // y at the magnet: there are no hits inside it, so interpolate ST1 -> ST2
    double y1 = 0, y2 = 0; bool h1 = false, h2 = false;
    for (size_t k = 0; k < zz.size(); ++k) {
      if (std::fabs(zz[k] - kStationZ[1]) < 250) { y1 = yy[k]; h1 = true; }
      if (std::fabs(zz[k] - kStationZ[2]) < 250) { y2 = yy[k]; h2 = true; }
    }
    if (!(h1 && h2)) continue;
    const double fr = (kMagnetZ - kStationZ[1]) / (kStationZ[2] - kStationZ[1]);

    out.z = zz; out.y = yy; out.nStations = nSt;
    out.yMag  = y1 + fr * (y2 - y1);
    out.inMap = std::fabs(out.yMag) <= kMapHalfY;
    out.ok    = true;

    printf("  %s\n", file);
    printf("     event %lld : %zu hits, %d/4 stations  ->  IN TRACKER ACCEPTANCE\n",
           (long long)i, zz.size(), nSt);
    printf("     y at magnet = %+.0f mm  ->  %s\n\n", out.yMag,
           out.inMap ? "inside field map   (bends: p measurable)"
                     : "OUTSIDE field map  (no bend: p UNMEASURABLE)");
    return out;
  }
  printf("  %s : no event with all 4 stations in the first %lld\n",
         file, (long long)n);
  return out;
}

void drawScene(const char* title, int tag) {
  auto* fr = new TH2F(Form("fr_%d", tag),
                      Form("%s;z  (mm, world frame);y  (mm)", title),
                      100, 18500, 34800, 100, -3400, 3400);
  fr->SetStats(0);
  fr->Draw();

  // ---- STATION ACCEPTANCE: the straws span the full +-3000 mm in y ----------
  for (int i = 0; i < 4; ++i) {
    auto* b = new TBox(kStationZ[i] - 110, -kAperture, kStationZ[i] + 110, kAperture);
    b->SetFillColorAlpha(kGray + 1, 0.22);
    b->SetLineColor(kGray + 2); b->SetLineWidth(1);
    b->Draw("l same");
    auto* l = new TLatex(kStationZ[i], -kAperture - 240, Form("ST%d", i));
    l->SetTextAlign(22); l->SetTextSize(0.028); l->SetTextColor(kGray + 3);
    l->Draw();
  }
  // the acceptance edges, drawn right across so the scale is unmistakable
  for (double a : {+kAperture, -kAperture}) {
    auto* l = new TLine(18500, a, 34800, a);
    l->SetLineColor(kGray + 2); l->SetLineStyle(1); l->SetLineWidth(2);
    l->Draw();
  }
  auto* apl = new TLatex(18800, kAperture + 130,
                         "straw station acceptance  #pm3000 mm");
  apl->SetTextSize(0.026); apl->SetTextColor(kGray + 3); apl->Draw();

  // the field map -- the whole story
  auto* m = new TBox(kMagnetZ - kMapHalfZ, -kMapHalfY,
                     kMagnetZ + kMapHalfZ,  kMapHalfY);
  m->SetFillColorAlpha(kBlue - 9, 0.35);
  m->SetLineColor(kBlue + 1); m->SetLineWidth(2);
  m->Draw("l same");
  auto* ml = new TLatex(kMagnetZ + 950, kMapHalfY + 380, "field map  #pm500 mm");
  ml->SetTextAlign(12); ml->SetTextSize(0.027); ml->SetTextColor(kBlue + 1);
  ml->Draw();
  auto* arr = new TLine(kMagnetZ + 900, kMapHalfY + 330, kMagnetZ + 300, kMapHalfY + 40);
  arr->SetLineColor(kBlue + 1); arr->Draw();
  // the mismatch, stated
  auto* gap = new TLatex(18800, -kAperture - 700,
      "#color[600]{field map covers only 1/6 of the station aperture}");
  gap->SetTextSize(0.026); gap->Draw();

  for (double s : {+kMapHalfY, -kMapHalfY}) {
    auto* l = new TLine(18500, s, 34800, s);
    l->SetLineColor(kBlue + 1); l->SetLineStyle(2);
    l->Draw();
  }

  auto* ax = new TLine(18500, 0, 34800, 0);
  ax->SetLineColor(kGray + 1); ax->SetLineStyle(3); ax->Draw();
  auto* g = new TGraph(1); g->SetPoint(0, kGunZ, 0);
  g->SetMarkerStyle(29); g->SetMarkerSize(1.7); g->Draw("P same");
}
}  // namespace

void evdisplay(const char* fileStraight = "build/pit_0.0_t-1.root",
               const char* fileTilted   = "build/pit_0.06_t-1.root") {
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

  printf("\n=== selecting events with the primary in FULL tracker acceptance ===\n\n");

  auto* c = new TCanvas("ev", "SHiP straw tracker -- bend plane", 1560, 660);
  c->Divide(2, 1);

  struct Case { const char* file; const char* title; int col; };
  Case cases[2] = {
    {fileStraight, "on-axis  #theta = 0,  5 GeV  #minus  inside map, bends",       kTeal + 2},
    {fileTilted,   "tilted  #theta = 0.06,  5 GeV  #minus  misses map, no bend",   kOrange + 8}
  };

  for (int p = 0; p < 2; ++p) {
    c->cd(p + 1);
    gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
    drawScene(cases[p].title, p);

    Ev e = findAcceptedEvent(cases[p].file);
    if (!e.ok) continue;

    auto* g = new TGraph((int)e.z.size(), e.z.data(), e.y.data());
    g->SetMarkerStyle(20); g->SetMarkerSize(0.6);
    g->SetMarkerColor(cases[p].col);
    g->SetLineColor(cases[p].col); g->SetLineWidth(2);
    g->Draw("PL same");

    auto* mk = new TGraph(1); mk->SetPoint(0, kMagnetZ, e.yMag);
    mk->SetMarkerStyle(24); mk->SetMarkerSize(2.2);
    mk->SetMarkerColor(cases[p].col); mk->SetLineWidth(2);
    mk->Draw("P same");

    auto* tx = new TLatex(kMagnetZ + 350, e.yMag,
                          Form("y = %+.0f mm  %s", e.yMag,
                               e.inMap ? "(in map)" : "(OUTSIDE)"));
    tx->SetTextSize(0.029); tx->SetTextColor(cases[p].col); tx->Draw();

    auto* acc = new TLatex(18800, -2400,
                           Form("#bf{%d/4 stations hit} #minus fully in tracker acceptance",
                                e.nStations));
    acc->SetTextSize(0.027); acc->SetTextColor(kGray + 3); acc->Draw();
  }

  c->SaveAs("evdisplay.png");
  printf("wrote evdisplay.png\n\n"
         "Both tracks are fully inside the TRACKER acceptance (4/4 stations).\n"
         "Only the on-axis one is inside the FIELD MAP -- that gap is the finding.\n");
}
