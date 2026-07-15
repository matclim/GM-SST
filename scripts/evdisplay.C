// tests/evdisplay.C — event display for the SHiP straw tracker.
//
//   root -l 'tests/evdisplay.C("dpvtx.root","straws.root","<fieldmap>.root",0)'
//
// Two side views of one event, drawn from the OUTPUT files:
//   * z-y (the bend plane -- the field bends in x, so the deflection shows here)
//   * z-x
// In each panel:
//   * the field map as a coloured |B| heatmap (from the map's Data tree)
//   * the four straw stations (from the straw table)
//   * the real hits of the event (straw address -> wire position via the table)
//   * the fitted tracks, drawn as rays from the fitted vertex along (theta,phi)
//   * the fitted vertex (star) and the truth vertex (open circle)
//   * the reconstructed parent line to the target (pointing)
//
// Coordinates: the reco writes vertices in the SHiP frame; the straw table and
// the field map are in the WORLD frame (= SHiP - kZOrigin). Everything is drawn
// in the WORLD frame, so the vertex is shifted back on the way in.

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TH2F.h>
#include <TBox.h>
#include <TLine.h>
#include <TMarker.h>
#include <TLatex.h>
#include <TStyle.h>
#include <vector>
#include <map>
#include <cmath>
#include <cstdio>

namespace {

constexpr double kZOrigin = 60000.0;
constexpr double kStationZ[4] = {84320., 86820., 91820., 93320.};  // SHiP
constexpr double kMagnetZ     = 89220.;  // SHiP
constexpr double kMapHalfZ    = 7500.;
constexpr double kMapHalfXY   = 5000.;
constexpr double kAperture    = 3000.;

struct View { int axis; const char* vlabel; const char* title; };
struct Wire { double cx, cy, cz; };

std::map<long,Wire> g_wires;
long wireKey(int st,int la,int su,int sw){ return ((((long)st*8 + la)*2 + su)*400 + sw); }

void loadStraws(const char* file) {
  TFile* f = TFile::Open(file);
  if (!f || f->IsZombie()) { printf("cannot open straw table %s\n", file); return; }
  TTree* t = (TTree*)f->Get("Straws");
  int st,la,su,sw; double cx,cy,cz;
  t->SetBranchAddress("station",&st); t->SetBranchAddress("layer",&la);
  t->SetBranchAddress("subLayer",&su); t->SetBranchAddress("straw",&sw);
  t->SetBranchAddress("cx",&cx); t->SetBranchAddress("cy",&cy); t->SetBranchAddress("cz",&cz);
  for (Long64_t i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    g_wires[wireKey(st,la,su,sw)] = {cx,cy,cz + kZOrigin}; }  // -> SHiP z
  printf("loaded %zu wires\n", g_wires.size());
  f->Close();
}

TH2F* fieldHeatmap(const char* /*mapFile*/, const View& /*v*/, int /*tag*/) {
  // Heatmap disabled for now -- the field REGION is still drawn as a box in
  // drawScene, so the bending is visible in context without the colour map.
  return nullptr;
}

void drawScene(const View& v, TH2F* heat, int tag) {
  auto* frame = new TH2F(Form("fr%d",tag),
      Form("%s;z  (mm, SHiP);%s  (mm)", v.title, v.vlabel),
      100, 30000, 96000, 100, -3400, 3400);   // SHiP z: decay volume + tracker
  frame->SetStats(0);
  frame->Draw();
  if (heat) {
    heat->SetContour(40);
    heat->SetMinimum(1e-3);         // floor
    heat->SetMaximum(0.5);          // cap: on-axis ~0.15 T stays visible against
                                    // the ~2 T off-axis peaks (which just saturate)
    heat->Draw("COLZ SAME");
    gPad->Update();
    frame->Draw("AXIS SAME");       // axes back on top of the colour
  }
  auto* mbox = new TBox(kMagnetZ-kMapHalfZ, -kMapHalfXY, kMagnetZ+kMapHalfZ, kMapHalfXY);
  mbox->SetFillStyle(0); mbox->SetLineColor(kAzure+1); mbox->SetLineStyle(2); mbox->Draw("l");
  for (int i=0;i<4;++i){
    auto* b=new TBox(kStationZ[i]-90,-kAperture,kStationZ[i]+90,kAperture);
    b->SetFillStyle(0); b->SetLineColor(kGray+2); b->Draw("l");
    auto* l=new TLatex(kStationZ[i],-kAperture-230,Form("ST%d",i));
    l->SetTextAlign(22); l->SetTextSize(0.028); l->SetTextColor(kGray+3); l->Draw();
  }
  auto* ax=new TLine(30000,0,96000,0); ax->SetLineColor(kGray+1); ax->SetLineStyle(3); ax->Draw();
}

}  // namespace

void evdisplay(const char* recoFile = "dpvtx.root",
               const char* strawFile = "straws.root",
               const char* mapFile   = "",
               int event = 0) {
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kBird);
  gStyle->SetNumberContours(40);

  loadStraws(strawFile);

  TFile* rf = TFile::Open(recoFile);
  if (!rf || rf->IsZombie()) { printf("cannot open reco file %s\n", recoFile); return; }
  TTree* vt = (TTree*)rf->Get("Vertices");
  TTree* tt = (TTree*)rf->Get("Tracks");
  TTree* ht = (TTree*)rf->Get("Hits");

  int    v_event, v_nTrk, v_fitOK;
  double v_vx,v_vy,v_vz, v_tx,v_ty,v_tz, v_ppx,v_ppy,v_ppz;
  vt->SetBranchAddress("event",&v_event); vt->SetBranchAddress("nTrk",&v_nTrk);
  vt->SetBranchAddress("fitOK",&v_fitOK);
  vt->SetBranchAddress("vx",&v_vx); vt->SetBranchAddress("vy",&v_vy); vt->SetBranchAddress("vz",&v_vz);
  vt->SetBranchAddress("tx",&v_tx); vt->SetBranchAddress("ty",&v_ty); vt->SetBranchAddress("tz",&v_tz);
  const bool hasParent = vt->GetBranch("parentPx");
  if (hasParent) {
    vt->SetBranchAddress("parentPx",&v_ppx); vt->SetBranchAddress("parentPy",&v_ppy);
    vt->SetBranchAddress("parentPz",&v_ppz);
  }

  bool found=false; double Vx=0,Vy=0,Vz=0,Tx=0,Ty=0,Tz=0,Ppx=0,Ppy=0,Ppz=0; int nTrk=0;
  for (Long64_t i=0;i<vt->GetEntries();++i){ vt->GetEntry(i);
    if (v_event!=event) continue;
    Vx=v_vx; Vy=v_vy; Vz=v_vz;            // tree is already SHiP
    Tx=v_tx; Ty=v_ty; Tz=v_tz;
    if (hasParent){ Ppx=v_ppx; Ppy=v_ppy; Ppz=v_ppz; }
    nTrk=v_nTrk; found=(v_fitOK==1); break; }
  if (!found) printf("event %d: no successful vertex (drawing what exists)\n", event);

  int    t_event; double t_theta,t_phi,t_p;
  // qFit has been stored as either int or double across reco versions -- detect.
  int    t_qi = 0; double t_qd = 0.0; int t_q = 0;
  TBranch* qb = tt->GetBranch("qFit");
  const bool qIsDouble = qb && TString(qb->GetLeaf("qFit")->GetTypeName())=="Double_t";
  tt->SetBranchAddress("event",&t_event);
  tt->SetBranchAddress("thetaFit",&t_theta); tt->SetBranchAddress("phiFit",&t_phi);
  tt->SetBranchAddress("pFit",&t_p);
  if (qIsDouble) tt->SetBranchAddress("qFit",&t_qd);
  else           tt->SetBranchAddress("qFit",&t_qi);
  struct Trk{ double dx,dy,dz,p; int q; };
  std::vector<Trk> trks;
  for (Long64_t i=0;i<tt->GetEntries();++i){ tt->GetEntry(i);
    if (t_event!=event) continue;
    t_q = qIsDouble ? (t_qd < 0 ? -1 : 1) : t_qi;
    trks.push_back({ std::sin(t_theta)*std::cos(t_phi),
                     std::sin(t_theta)*std::sin(t_phi),
                     std::cos(t_theta), t_p, t_q }); }

  int h_event,h_st,h_la,h_su,h_sw;
  ht->SetBranchAddress("event",&h_event); ht->SetBranchAddress("station",&h_st);
  ht->SetBranchAddress("layer",&h_la); ht->SetBranchAddress("sub",&h_su);
  ht->SetBranchAddress("straw",&h_sw);
  std::vector<Wire> hits;
  for (Long64_t i=0;i<ht->GetEntries();++i){ ht->GetEntry(i);
    if (h_event!=event) continue;
    auto it=g_wires.find(wireKey(h_st,h_la,h_su,h_sw));
    if (it!=g_wires.end()) hits.push_back(it->second); }

  printf("event %d:  nTrk=%d  tracks=%zu  hits=%zu  vtx(SHiP)=(%.0f,%.0f,%.0f)\n",
         event, nTrk, trks.size(), hits.size(), Vx, Vy, Vz);

  View views[2] = { {1,"y","z-y  (bend plane)"}, {0,"x","z-x"} };
  auto* c = new TCanvas("ev","SHiP straw tracker event display",1500,900);
  c->Divide(1,2);

  for (int p=0;p<2;++p){
    c->cd(p+1); gPad->SetRightMargin(0.12); gPad->SetLeftMargin(0.08);
    const View& v = views[p];
    TH2F* heat = fieldHeatmap(mapFile, v, p);
    drawScene(v, heat, p);

    auto vert=[&](double x,double y){ return (v.axis==1)?y:x; };

    if (!hits.empty()){
      auto* g=new TGraph((int)hits.size());
      for (size_t i=0;i<hits.size();++i) g->SetPoint(i, hits[i].cz, vert(hits[i].cx,hits[i].cy));
      g->SetMarkerStyle(20); g->SetMarkerSize(0.5); g->SetMarkerColor(kBlack); g->Draw("P");
    }

    // representative on-axis field for the bending-radius estimate
    const double kBrepr = 0.15;   // T (on-axis; |B| peaks ~2 T off-axis)
    int trkIdx = 0;
    for (auto& tr : trks){
      const double L = 62000.0;   // reach from the decay volume to the tracker
      auto* line=new TLine(Vz, vert(Vx,Vy), Vz+tr.dz*L, vert(Vx,Vy)+vert(tr.dx,tr.dy)*L);
      line->SetLineColor(tr.q>0?kRed+1:kAzure+2); line->SetLineWidth(2); line->Draw();

      if (p==0) {   // annotate momentum + bending radius, top panel only
        const double R = tr.p / (0.3 * kBrepr);         // metres
        // stagger labels vertically so they don't overlap near the axis
        const double zlab = 44000.0 + trkIdx*9000.0;
        const double ylab = 700.0 + trkIdx*450.0;
        auto* lab = new TLatex(zlab, ylab,
          Form("#color[%d]{q%+d  p=%.0f GeV  R#approx%.0f m}",
               (tr.q>0?kRed+1:kAzure+2), tr.q, tr.p, R));
        lab->SetTextSize(0.033); lab->Draw();
        // thin dotted leader from the label down to the track at that z
        const double ytrk = vert(Vx,Vy) + vert(tr.dx,tr.dy)*(zlab-Vz);
        auto* lead = new TLine(zlab+3000.0, ylab, zlab+3000.0, ytrk);
        lead->SetLineColor(tr.q>0?kRed+1:kAzure+2);
        lead->SetLineStyle(3); lead->Draw();
      }
      ++trkIdx;
    }

    if (Ppx||Ppy||Ppz){
      const double pn=std::sqrt(Ppx*Ppx+Ppy*Ppy+Ppz*Ppz);
      const double ux=Ppx/pn, uy=Ppy/pn, uz=Ppz/pn;
      const double L=Vz-30000.0;   // back toward the target, to the frame edge
      auto* line=new TLine(Vz, vert(Vx,Vy), Vz-uz*L, vert(Vx,Vy)-vert(ux,uy)*L);
      line->SetLineColor(kGreen+2); line->SetLineStyle(2); line->SetLineWidth(1); line->Draw();
    }

    auto* mv=new TMarker(Vz, vert(Vx,Vy), 29); mv->SetMarkerColor(kOrange+7);
    mv->SetMarkerSize(2.4); mv->Draw();
    auto* mt=new TMarker(Tz, vert(Tx,Ty), 24); mt->SetMarkerColor(kGreen+3);
    mt->SetMarkerSize(1.8); mt->SetMarkerStyle(24); mt->Draw();

    if (p==0){
      auto* leg=new TLatex(31000,3050,
        Form("event %d   n_{trk}=%d   #star fitted vtx   #circ truth vtx", event,nTrk));
      leg->SetTextSize(0.030); leg->Draw();
    }
  }

  c->SaveAs(Form("evdisplay_ev%d.png",event));
  printf("wrote evdisplay_ev%d.png\n", event);
}
