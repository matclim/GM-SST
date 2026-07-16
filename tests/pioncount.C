void pioncount(const char* fn="build/ks20_t-1.root") {
  TFile f(fn);
  TTree* T = (TTree*)f.Get("Events");
  std::vector<int> *trackID=nullptr,*parentID=nullptr,*pdg=nullptr;
  T->SetBranchAddress("trackID",&trackID);
  T->SetBranchAddress("parentID",&parentID);
  T->SetBranchAddress("pdg",&pdg);
  int n2=0,n1=0,n0=0;
  for (Long64_t i=0;i<T->GetEntries();++i) {
    T->GetEntry(i);
    std::set<int> ids;
    for (size_t j=0;j<trackID->size();++j)
      if ((*parentID)[j]==1 && abs((*pdg)[j])==211) ids.insert((*trackID)[j]);
    printf("event %lld: %zu pion tracks\n", i, ids.size());
    if (ids.size()>=2) ++n2; else if (ids.size()==1) ++n1; else ++n0;
  }
  printf("\nSUMMARY: 2+ pions: %d   1 pion: %d   0 pions: %d  (of %lld)\n",
         n2, n1, n0, T->GetEntries());
}
