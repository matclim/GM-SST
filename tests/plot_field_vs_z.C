// plot_field_vs_z.C
// Plot Bx vs z and By vs z along the central axis (x = 0, y = 0) of a
// ROOT field map (the "Data" TTree with branches x,y,z,Bx,By,Bz).
//
// Run (interactive):   root -l 'plot_field_vs_z.C("yourmap.root")'
// Run (batch, no GUI): root -l -b -q 'plot_field_vs_z.C("yourmap.root")'
// Compile with ACLiC:  root -l 'plot_field_vs_z.C+("yourmap.root")'

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <vector>

void plot_field_vs_z(const char* fname = "field.root") {
    TFile* f = TFile::Open(fname, "READ");
    if (!f || f->IsZombie()) { std::cerr << "Cannot open " << fname << "\n"; return; }

    TTree* data = (TTree*)f->Get("Data");
    if (!data) { std::cerr << "No 'Data' tree in " << fname << "\n"; return; }

    // Tolerance for "x = 0" / "y = 0": a quarter of the grid step, taken from
    // the Range tree if present, else 1 mm. (Grid nodes are usually exact.)
    double tol = 1.0;
    if (TTree* rng = (TTree*)f->Get("Range")) {
        Float_t dx = 0, dy = 0;
        rng->SetBranchAddress("dx", &dx);
        rng->SetBranchAddress("dy", &dy);
        rng->GetEntry(0);
        if (dx > 0 && dy > 0) tol = 0.25 * std::min(dx, dy);
    }

    Float_t x, y, z, Bx, By, Bz;
    data->SetBranchAddress("x",  &x);
    data->SetBranchAddress("y",  &y);
    data->SetBranchAddress("z",  &z);
    data->SetBranchAddress("Bx", &Bx);
    data->SetBranchAddress("By", &By);
    data->SetBranchAddress("Bz", &Bz);

    // Collect the on-axis line (x = 0, y = 0).  ← change these two cuts to slice
    // at a different transverse position, e.g. std::fabs(y - 500) < tol.
    std::vector<std::array<double, 3>> pts;  // {z, Bx, By}
    const Long64_t n = data->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        data->GetEntry(i);
        if (std::fabs(x) < tol && std::fabs(y) < tol)
            pts.push_back({(double)z, (double)Bx, (double)By});
    }
    if (pts.empty()) { std::cerr << "No points found with x = y = 0.\n"; return; }

    // Sort by z so the connecting line is monotonic.
    std::sort(pts.begin(), pts.end(),
              [](const auto& a, const auto& b) { return a[0] < b[0]; });

    TGraph* gBx = new TGraph((int)pts.size());
    TGraph* gBy = new TGraph((int)pts.size());
    for (size_t i = 0; i < pts.size(); ++i) {
        gBx->SetPoint((int)i, pts[i][0], pts[i][1]);
        gBy->SetPoint((int)i, pts[i][0], pts[i][2]);
    }

    gBx->SetTitle("B_{x} vs z at x = y = 0;z [mm];B_{x} [T]");
    gBy->SetTitle("B_{y} vs z at x = y = 0;z [mm];B_{y} [T]");
    gBx->SetLineColor(kBlue + 1); gBx->SetLineWidth(2);
    gBy->SetLineColor(kRed + 1);  gBy->SetLineWidth(2);
    gBx->SetMarkerStyle(20); gBx->SetMarkerSize(0.4); gBx->SetMarkerColor(kBlue + 1);
    gBy->SetMarkerStyle(20); gBy->SetMarkerSize(0.4); gBy->SetMarkerColor(kRed + 1);

    TCanvas* c = new TCanvas("cField", "Field along z", 1100, 450);
    c->Divide(2, 1);
    c->cd(1); gBx->Draw("ALP");
    c->cd(2); gBy->Draw("ALP");
    c->SaveAs("field_vs_z.png");

    std::cout << "Plotted " << pts.size() << " on-axis points; "
              << "wrote field_vs_z.png\n";
}
