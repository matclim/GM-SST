// tests/summary.C — tabulate reco quality for every reco_*.root in a dir.
//   root -l -b -q 'summary.C(".")'
double stat(TTree* t, const char* expr, bool rms) {
  t->Draw(Form("%s>>h_%p", expr, (void*)expr), "", "goff");
  TH1* h = (TH1*)gDirectory->Get(Form("h_%p", (void*)expr));
  double v = h ? (rms ? h->GetRMS() : h->GetMean()) : 0.0;
  if (h) delete h;
  return v;
}
void summary(const char* dir = ".") {
  TSystemDirectory sdir(dir, dir);
  TList* files = sdir.GetListOfFiles();
  if (!files) { printf("no files in %s\n", dir); return; }
  std::vector<TString> names;
  TIter next(files); TSystemFile* f;
  while ((f = (TSystemFile*)next())) {
    TString n = f->GetName();
    if (!f->IsDirectory() && n.BeginsWith("reco_") && n.EndsWith(".root"))
      names.push_back(n);
  }
  std::sort(names.begin(), names.end());
  printf("%-22s %5s  %8s %8s  %10s  %10s\n",
         "file", "N", "meanRes", "rmsRes", "<sigmaP>", "<chi2>/ndf");
  for (auto& n : names) {
    TString path = TString(dir) + "/" + n;
    TFile fin(path); TTree* t = (TTree*)fin.Get("Tracks");
    if (!t) { printf("%-22s  no tree\n", n.Data()); continue; }
    double nd = stat(t, "nMeas", false) - 5;
    printf("%-22s %5lld  %+8.3f %8.3f  %10.3f  %10.2f\n",
           n.Data(), (long long)t->GetEntries(),
           stat(t, "pRes", false), stat(t, "pRes", true),
           stat(t, "sigmaP", false), nd > 0 ? stat(t, "chi2", false) / nd : 0.0);
  }
}
