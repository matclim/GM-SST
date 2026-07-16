// check_bending.C  (v2 - primary-track only, line-fit slopes)
// Verify spectrometer bending from a StrawTracker run.
//
//   root -l 'check_bending.C("field_muon.root", "yourmap.root", 10000)'
//     arg1: hits file (Events tree)
//     arg2: field map (optional) -> analytic cross-check
//     arg3: muon kinetic energy in MeV (must match the run; default 10000)
//
// Uses ONLY the primary muon (trackID==1) and fits a straight line to its hits
// upstream (stations 0,1) and downstream (stations 2,3). The bend is the change
// in slope. A dipole along x deflects a +z muon in y: dy' != 0, dx' ~ 0.
// Reports the error on the mean so you can judge significance.

#include "TFile.h"
#include "TTree.h"
#include <cmath>
#include <cstdio>
#include <vector>

namespace {
double pGeV(double keMeV){ const double m=105.6583745; double E=keMeV+m;
                          return std::sqrt(E*E-m*m)/1000.0; }
// least-squares slope of v vs z; returns false if z has no spread
bool slope(const std::vector<double>& z,const std::vector<double>& v,double& s){
    size_t n=z.size(); if(n<2) return false;
    double zb=0,vb=0; for(size_t i=0;i<n;++i){zb+=z[i];vb+=v[i];} zb/=n; vb/=n;
    double szz=0,szv=0; for(size_t i=0;i<n;++i){szz+=(z[i]-zb)*(z[i]-zb); szv+=(z[i]-zb)*(v[i]-vb);}
    if(szz<=0) return false; s=szv/szz; return true;
}
}

void check_bending(const char* hitsFile="field_muon.root",
                   const char* mapFile=nullptr, double keMeV=10000.0){
    TFile* fh=TFile::Open(hitsFile,"READ");
    if(!fh||fh->IsZombie()){printf("Cannot open %s\n",hitsFile);return;}
    TTree* T=(TTree*)fh->Get("Events");
    if(!T){printf("No 'Events' tree\n");return;}

    std::vector<int>* trackID=nullptr; std::vector<int>* stationID=nullptr;
    std::vector<double>* x=nullptr; std::vector<double>* y=nullptr; std::vector<double>* z=nullptr;
    T->SetBranchAddress("trackID",&trackID);
    T->SetBranchAddress("stationID",&stationID);
    T->SetBranchAddress("x",&x); T->SetBranchAddress("y",&y); T->SetBranchAddress("z",&z);

    double sY=0,sY2=0,sX=0,sX2=0; long nUsed=0;
    const Long64_t nEv=T->GetEntries();
    for(Long64_t e=0;e<nEv;++e){
        T->GetEntry(e);
        std::vector<double> zu,yu,xu, zd,yd,xd;
        for(size_t h=0;h<stationID->size();++h){
            if((*trackID)[h]!=1) continue;           // primary muon only
            int s=(*stationID)[h]; if(s<0||s>3) continue;
            if(s<2){ zu.push_back((*z)[h]); yu.push_back((*y)[h]); xu.push_back((*x)[h]); }
            else   { zd.push_back((*z)[h]); yd.push_back((*y)[h]); xd.push_back((*x)[h]); }
        }
        double syu,syd,sxu,sxd;
        if(!slope(zu,yu,syu)||!slope(zd,yd,syd)) continue;
        if(!slope(zu,xu,sxu)||!slope(zd,xd,sxd)) continue;
        double dY=syd-syu, dX=sxd-sxu;
        sY+=dY; sY2+=dY*dY; sX+=dX; sX2+=dX*dX; ++nUsed;
    }
    if(!nUsed){printf("No usable events (primary hits in up- and down-stream).\n");return;}
    double mY=sY/nUsed, rY=std::sqrt(std::max(0.0,sY2/nUsed-mY*mY)), eY=rY/std::sqrt((double)nUsed);
    double mX=sX/nUsed, rX=std::sqrt(std::max(0.0,sX2/nUsed-mX*mX)), eX=rX/std::sqrt((double)nUsed);

    printf("\n===== BENDING CHECK (primary track): %s =====\n",hitsFile);
    printf("events used         : %ld / %lld\n",nUsed,(long long)nEv);
    printf("bend dy'            : %+.3e +/- %.1e rad  (rms %.1e)  significance %.1f sigma\n",
           mY,eY,rY, eY>0?std::fabs(mY)/eY:0);
    printf("bend dx'            : %+.3e +/- %.1e rad  (rms %.1e)  significance %.1f sigma\n",
           mX,eX,rX, eX>0?std::fabs(mX)/eX:0);

    if(mapFile){
        TFile* fm=TFile::Open(mapFile,"READ");
        if(fm&&!fm->IsZombie()){
            TTree* R=(TTree*)fm->Get("Range"); TTree* D=(TTree*)fm->Get("Data");
            if(R&&D){
                Float_t dx,dy,dz,zmin,zmax; R->SetBranchAddress("dx",&dx);R->SetBranchAddress("dy",&dy);
                R->SetBranchAddress("dz",&dz);R->SetBranchAddress("zMin",&zmin);R->SetBranchAddress("zMax",&zmax);
                R->GetEntry(0);
                Float_t xx,yy,zz,bx,by,bz;
                D->SetBranchAddress("x",&xx);D->SetBranchAddress("y",&yy);D->SetBranchAddress("z",&zz);
                D->SetBranchAddress("Bx",&bx);D->SetBranchAddress("By",&by);D->SetBranchAddress("Bz",&bz);
                double sBx=0; long n=0; const double tx=0.5*dx,ty=0.5*dy;
                for(Long64_t i=0;i<D->GetEntries();++i){D->GetEntry(i);
                    if(std::fabs(xx)<tx&&std::fabs(yy)<ty){sBx+=bx;++n;}}
                double IBx=sBx*dz*1e-3, p=pGeV(keMeV), expY=-0.299792458*IBx/p;
                printf("\nmap Int Bx dl       : %+.4g T*m   (map z half-extent = %.0f mm)\n",
                       IBx, 0.5*(zmax-zmin));
                printf("field box half-z    : 1750 mm   <- field only applied inside this box\n");
                printf("momentum            : %.4g GeV/c\n",p);
                printf("predicted dy'       : %+.3e rad\n",expY);
                printf("measured/predicted  : %.2f  (measured error is %.0f%% of prediction)\n",
                       expY!=0?mY/expY:0, expY!=0?100*eY/std::fabs(expY):0);
                if(0.5*(zmax-zmin) > 1750.0)
                    printf("NOTE: map is longer than the field box -> the muon only sees the\n"
                           "      central +/-1750 mm, so measured < full-map prediction is expected.\n");
            }
        }
    }
    printf("=========================================\n\n");
}
