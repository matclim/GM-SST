// check_orientation.C
// Diagnostics to verify a SHiP ROOT field map (Range + Data trees) is
// correctly oriented and in the units your reader expects.
//
// Run:  root -l 'check_orientation.C("yourmap.root", 0.15)'
//         2nd arg = expected on-axis peak |B| in Tesla (your dipole strength).
//
// It reports, and flags, the seven checks discussed:
//   1 beam axis (longest dimension)      2 length units (cm vs mm)
//   3 field units (T vs kG)              4 dominant on-axis component
//   5 fringe profile + z-symmetry        6 sign of dominant component
//   7 Maxwell closure  div B ~ 0  in the field region
//
// Nothing here assumes a particular unit; it prints both interpretations and
// lets the physics (magnet size, dipole strength, div B) decide.

#include "TFile.h"
#include "TTree.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {
struct Grid {
    int    N[3];
    double lo[3], hi[3], d[3];
    long   idx(int i, int j, int k) const {
        return (long(k) * N[1] + j) * N[0] + i;   // z-major, matches ShipFieldMap
    }
};
const char* AX[3] = {"x", "y", "z"};
}

void check_orientation(const char* fname = "field.root",
                       double expPeakTesla = 0.15) {
    TFile* f = TFile::Open(fname, "READ");
    if (!f || f->IsZombie()) { printf("Cannot open %s\n", fname); return; }

    // ---- Grid geometry from the Range tree --------------------------------
    Grid g{};
    if (TTree* R = (TTree*)f->Get("Range")) {
        Float_t xMin,xMax,dx,yMin,yMax,dy,zMin,zMax,dz;
        R->SetBranchAddress("xMin",&xMin); R->SetBranchAddress("xMax",&xMax); R->SetBranchAddress("dx",&dx);
        R->SetBranchAddress("yMin",&yMin); R->SetBranchAddress("yMax",&yMax); R->SetBranchAddress("dy",&dy);
        R->SetBranchAddress("zMin",&zMin); R->SetBranchAddress("zMax",&zMax); R->SetBranchAddress("dz",&dz);
        R->GetEntry(0);
        g.lo[0]=xMin; g.hi[0]=xMax; g.d[0]=dx;
        g.lo[1]=yMin; g.hi[1]=yMax; g.d[1]=dy;
        g.lo[2]=zMin; g.hi[2]=zMax; g.d[2]=dz;
        for (int a=0;a<3;++a)
            g.N[a] = (g.d[a]>0) ? int(std::lround((g.hi[a]-g.lo[a])/g.d[a]))+1 : 1;
    } else { printf("No 'Range' tree\n"); return; }

    TTree* D = (TTree*)f->Get("Data");
    if (!D) { printf("No 'Data' tree\n"); return; }

    const long Ntot = long(g.N[0])*g.N[1]*g.N[2];
    printf("\n===== FIELD MAP ORIENTATION CHECK: %s =====\n", fname);
    printf("grid            : %d x %d x %d = %ld points  (Data entries: %lld)\n",
           g.N[0],g.N[1],g.N[2],Ntot,(long long)D->GetEntries());
    for (int a=0;a<3;++a)
        printf("%s range        : [%.3g, %.3g], step %.3g  ->  %.4g wide "
               "(= %.4g cm  OR  %.4g mm)\n",
               AX[a], g.lo[a], g.hi[a], g.d[a], g.hi[a]-g.lo[a],
               g.hi[a]-g.lo[a], g.hi[a]-g.lo[a]);

    // ---- CHECK 1: beam axis = longest dimension ---------------------------
    int lng = 0; for (int a=1;a<3;++a) if (g.N[a] > g.N[lng]) lng = a;
    printf("\n[1] longest axis  : %s  %s\n", AX[lng],
           lng==2 ? "(OK - z is the beam axis)"
                  : "(CHECK - expected z to be longest; axes may be swapped)");

    // ---- Load field into dense arrays (bin by nearest index) --------------
    std::vector<float> Bx(Ntot,0), By(Ntot,0), Bz(Ntot,0);
    std::vector<char>  hit(Ntot,0);
    Float_t x,y,z,bx,by,bz;
    D->SetBranchAddress("x",&x); D->SetBranchAddress("y",&y); D->SetBranchAddress("z",&z);
    D->SetBranchAddress("Bx",&bx); D->SetBranchAddress("By",&by); D->SetBranchAddress("Bz",&bz);
    double gmax[3]={0,0,0}, gmaxB=0; float bmaxAt[3]={0,0,0};
    const long nE = D->GetEntries();
    for (long e=0;e<nE;++e) {
        D->GetEntry(e);
        int i = int(std::lround((x-g.lo[0])/g.d[0]));
        int j = int(std::lround((y-g.lo[1])/g.d[1]));
        int k = int(std::lround((z-g.lo[2])/g.d[2]));
        if (i<0||i>=g.N[0]||j<0||j>=g.N[1]||k<0||k>=g.N[2]) continue;
        long id=g.idx(i,j,k); Bx[id]=bx; By[id]=by; Bz[id]=bz; hit[id]=1;
        gmax[0]=std::max(gmax[0],std::fabs((double)bx));
        gmax[1]=std::max(gmax[1],std::fabs((double)by));
        gmax[2]=std::max(gmax[2],std::fabs((double)bz));
        double mag=std::sqrt(double(bx)*bx+double(by)*by+double(bz)*bz);
        if (mag>gmaxB){gmaxB=mag; bmaxAt[0]=x; bmaxAt[1]=y; bmaxAt[2]=z;}
    }
    long miss=0; for (char h:hit) if(!h) ++miss;
    if (miss) printf("    note: %ld grid points had no data (read as 0)\n", miss);

    // ---- CHECK 3: field units (peak |B|) ----------------------------------
    printf("\n[3] peak |B|      : %.4g   at (x,y,z)=(%.3g, %.3g, %.3g)\n",
           gmaxB, bmaxAt[0],bmaxAt[1],bmaxAt[2]);
    printf("    per-component max: |Bx|=%.4g  |By|=%.4g  |Bz|=%.4g\n",
           gmax[0],gmax[1],gmax[2]);
    printf("    if file is Tesla -> peak %.3g T ; if kilogauss -> peak %.3g T\n",
           gmaxB, gmaxB*0.1);
    { double rT=gmaxB/expPeakTesla, rG=gmaxB*0.1/expPeakTesla;
      printf("    vs expected %.3g T:  Tesla gives x%.2g , kG gives x%.2g  -> %s\n",
             expPeakTesla, rT, rG,
             (std::fabs(rT-1)<std::fabs(rG-1)) ? "file looks like TESLA"
                                               : "file looks like KILOGAUSS"); }

    // ---- On-axis line x=y=0 ------------------------------------------------
    int i0 = int(std::lround((0-g.lo[0])/g.d[0]));
    int j0 = int(std::lround((0-g.lo[1])/g.d[1]));
    i0=std::min(std::max(i0,0),g.N[0]-1); j0=std::min(std::max(j0,0),g.N[1]-1);
    printf("\n    on-axis line taken at x=%.3g, y=%.3g (nearest nodes to 0)\n",
           g.lo[0]+i0*g.d[0], g.lo[1]+j0*g.d[1]);

    double axMax[3]={0,0,0}; int axPeakK=0; double axPeakVal=0;
    for (int k=0;k<g.N[2];++k){
        double c[3]={Bx[g.idx(i0,j0,k)],By[g.idx(i0,j0,k)],Bz[g.idx(i0,j0,k)]};
        for (int a=0;a<3;++a) axMax[a]=std::max(axMax[a],std::fabs(c[a]));
    }
    int dom=0; for (int a=1;a<3;++a) if (axMax[a]>axMax[dom]) dom=a;
    // peak location & sign of dominant component
    for (int k=0;k<g.N[2];++k){
        double v = (dom==0?Bx:dom==1?By:Bz)[g.idx(i0,j0,k)];
        if (std::fabs(v)>std::fabs(axPeakVal)){axPeakVal=v; axPeakK=k;}
    }

    // ---- CHECK 4: dominant component --------------------------------------
    printf("\n[4] on-axis peaks : |Bx|=%.4g  |By|=%.4g  |Bz|=%.4g\n",
           axMax[0],axMax[1],axMax[2]);
    double sub = std::max({axMax[0],axMax[1],axMax[2]});
    sub = ( (dom==0?std::max(axMax[1],axMax[2]):
             dom==1?std::max(axMax[0],axMax[2]):
                    std::max(axMax[0],axMax[1])) ) / (axMax[dom]>0?axMax[dom]:1);
    printf("    dominant = B%s  (others are %.1f%% of it)  %s\n",
           AX[dom], 100*sub,
           dom==0 ? "(OK - dipole along x, matches the -0.15 T convention)"
                  : "(CHECK - expected Bx to dominate; field axes may be rotated)");

    // ---- CHECK 6: sign -----------------------------------------------------
    printf("\n[6] dominant sign : B%s peaks at %+.4g  (z=%.3g)\n",
           AX[dom], axPeakVal, g.lo[2]+axPeakK*g.d[2]);

    // ---- CHECK 5: fringe profile + z-symmetry -----------------------------
    // even-symmetry of dominant component about the map centre in z
    bool zSym = std::fabs(g.lo[2]+g.hi[2]) < 0.51*g.d[2];   // zMin ~ -zMax ?
    double endMax=0, asym=0;
    const std::vector<float>& Bd = (dom==0?Bx:dom==1?By:Bz);
    for (int k=0;k<g.N[2];++k){
        double v=Bd[g.idx(i0,j0,k)];
        if (k<2 || k>g.N[2]-3) endMax=std::max(endMax,std::fabs(v));
        if (zSym){ double vm=Bd[g.idx(i0,j0,g.N[2]-1-k)];
                   asym=std::max(asym, std::fabs(v-vm)); }
    }
    printf("\n[5] fringe        : |B%s| at the z-ends = %.4g  (%.1f%% of peak)  %s\n",
           AX[dom], endMax, 100*endMax/(axMax[dom]>0?axMax[dom]:1),
           endMax < 0.2*axMax[dom] ? "(OK - falls off at ends)"
                                   : "(CHECK - not falling off; origin or extent may be wrong)");
    if (zSym)
        printf("    z-symmetry    : max |B%s(z) - B%s(-z)| = %.4g  (%.1f%% of peak)  %s\n",
               AX[dom],AX[dom], asym, 100*asym/(axMax[dom]>0?axMax[dom]:1),
               asym < 0.1*axMax[dom] ? "(OK - even about centre)"
                                     : "(CHECK - asymmetric; possible index/origin issue)");
    else
        printf("    z-symmetry    : skipped (z range not centred on 0)\n");

    // ---- CHECK 7: Maxwell closure  div B ~ 0  in the strong-field region --
    // Sampled where |B| > 30%% of peak (most likely source-free gap).
    double sumDiv2=0, sumScale2=0; long nDiv=0;
    for (int k=1;k<g.N[2]-1;++k)
      for (int j=1;j<g.N[1]-1;++j)
        for (int i=1;i<g.N[0]-1;++i){
            long id=g.idx(i,j,k);
            double mag=std::sqrt(double(Bx[id])*Bx[id]+double(By[id])*By[id]+double(Bz[id])*Bz[id]);
            if (mag < 0.3*gmaxB) continue;
            double dBxdx=(Bx[g.idx(i+1,j,k)]-Bx[g.idx(i-1,j,k)])/(2*g.d[0]);
            double dBydy=(By[g.idx(i,j+1,k)]-By[g.idx(i,j-1,k)])/(2*g.d[1]);
            double dBzdz=(Bz[g.idx(i,j,k+1)]-Bz[g.idx(i,j,k-1)])/(2*g.d[2]);
            double div=dBxdx+dBydy+dBzdz;
            // scale = typical magnitude of the individual gradient terms
            double scale=std::fabs(dBxdx)+std::fabs(dBydy)+std::fabs(dBzdz);
            sumDiv2+=div*div; sumScale2+=scale*scale; ++nDiv;
        }
    if (nDiv){
        double rms = std::sqrt(sumDiv2/nDiv);
        double scl = std::sqrt(sumScale2/nDiv);
        double rel = (scl>0)? rms/scl : 0;
        printf("\n[7] div B         : RMS(div B)/RMS(|dB|) = %.3g  over %ld points  %s\n",
               rel, nDiv,
               rel < 0.05 ? "(OK - divergence ~ 0, field is consistent)"
                          : "(CHECK - large; a component may be sign-flipped or an axis swapped)");
    } else printf("\n[7] div B         : no strong-field interior points to sample\n");

    printf("\n===== end of report =====\n"
           "Reminder: [2] length units are settled by comparing the ranges above\n"
           "to the real magnet size; [3] settles the field units.\n\n");
}
