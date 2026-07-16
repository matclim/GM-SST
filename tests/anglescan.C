// tests/anglescan.C — compare reco quality vs incident angle (fixed momentum).
//   root -l -b -q 'tests/anglescan.C("build")'
double stat(TTree* t, const char* expr, const char* cut, bool rms) {
  static int n = 0;
  TString h = Form("h_%d", n++);
  t->Draw(Form("%s>>%s", expr, h.Data()), cut, "goff");
  TH1* hp = (TH1*)gDirectory->Get(h);
  if (!hp || hp->GetEntries() == 0) return 0.0;
  double v = rms ? hp->GetRMS() : hp->GetMean();
  return v;
}
// core sigma from a Gaussian fit (robust to the stuck-at-seed tail)
double coreSigma(TTree* t, const char* expr, double lo, double hi) {
  static int n = 0;
  TString h = Form("g_%d", n++);
  t->Draw(Form("%s>>%s(200,%f,%f)", expr, h.Data(), 4*lo, 4*hi), "", "goff");
  TH1* hp = (TH1*)gDirectory->Get(h);
  if (!hp || hp->GetEntries() < 10) return 0.0;
  hp->Fit("gaus", "Q0", "", lo, hi);
  TF1* f = hp->GetFunction("gaus");
  return f ? f->GetParameter(2) : 0.0;
}
void anglescan(const char* dir = "build") {
  const char* thetas[] = {"0.0", "0.02", "0.04", "0.06"};
  printf("%-8s %6s  %9s %9s %9s  %9s  %8s\n",
         "theta", "N", "meanPRes", "rmsPRes", "corePRes", "<bestChi2>", "stuck%");
  for (auto th : thetas) {
    TString path = Form("%s/pitreco_%s.root", dir, th);
    TFile f(path);
    if (f.IsZombie()) { printf("%-8s  (missing)\n", th); continue; }
    TTree* t = (TTree*)f.Get("Tracks");
    if (!t || t->GetEntries() == 0) { printf("%-8s  (no tracks)\n", th); continue; }
    Long64_t n     = t->GetEntries();
    Long64_t stuck = t->GetEntries("abs(pRes)>0.4");    // parked / badly wrong
    printf("%-8s %6lld  %+9.3f %9.3f %9.3f  %9.2f  %7.1f%%\n",
           th, n,
           stat(t, "pRes", "", false),
           stat(t, "pRes", "", true),
           coreSigma(t, "pRes", -0.3, 0.3),
           stat(t, "bestChi2", "", false),
           100.0 * stuck / n);
  }
}
