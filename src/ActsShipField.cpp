// src/ActsShipField.cpp
// Interpolated ACTS field from a SHiP ROOT map (Range+Data), wrapped so lookups
// outside the map return B = 0. The wrapper bounds the query against the TRUE
// data extent (xMin+origin .. xMax+origin), inset by one full grid cell -- this
// is immune to fieldMapXYZ's phantom boundary bin (which makes getMax() report
// one cell beyond the data) and to boundary floating-point overshoot. The lost
// outer cell carries only negligible fringe field.
#include "ActsShipField.hpp"

#include <Acts/Definitions/Algebra.hpp>
#include <Acts/Definitions/Units.hpp>
#include <Acts/MagneticField/BFieldMapUtils.hpp>
#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>
#include <Acts/MagneticField/MagneticFieldProvider.hpp>

#include <TFile.h>
#include <TTree.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <vector>

namespace shipreco {

namespace {

template <typename Map>
class ZeroOutsideField final : public Acts::MagneticFieldProvider {
 public:
  ZeroOutsideField(Map map, std::array<double, 3> lo, std::array<double, 3> hi)
      : m_map(std::move(map)), m_lo(lo), m_hi(hi) {}

  Acts::MagneticFieldProvider::Cache makeCache(
      const Acts::MagneticFieldContext& mctx) const override {
    return m_map.makeCache(mctx);
  }

  Acts::Result<Acts::Vector3> getField(
      const Acts::Vector3& position,
      Acts::MagneticFieldProvider::Cache& cache) const override {
    for (int i = 0; i < 3; ++i) {
      if (position[i] < m_lo[i] || position[i] > m_hi[i]) {
        return Acts::Result<Acts::Vector3>::success(Acts::Vector3::Zero());
      }
    }
    return m_map.getField(position, cache);
  }

 private:
  Map m_map;
  std::array<double, 3> m_lo{};
  std::array<double, 3> m_hi{};
};

}  // namespace

std::shared_ptr<const Acts::MagneticFieldProvider>
makeShipFieldFromRootMap(const std::string& rootFile,
                         std::array<double, 3> originMm, double bScaleToTesla) {
  std::unique_ptr<TFile> f(TFile::Open(rootFile.c_str(), "READ"));
  if (!f || f->IsZombie())
    throw std::runtime_error("ActsShipField: cannot open '" + rootFile + "'");
  auto* R = dynamic_cast<TTree*>(f->Get("Range"));
  auto* D = dynamic_cast<TTree*>(f->Get("Data"));
  if (!R || !D)
    throw std::runtime_error("ActsShipField: '" + rootFile + "' missing Range/Data");

  Float_t xMin,xMax,dx,yMin,yMax,dy,zMin,zMax,dz;
  R->SetBranchAddress("xMin",&xMin); R->SetBranchAddress("xMax",&xMax); R->SetBranchAddress("dx",&dx);
  R->SetBranchAddress("yMin",&yMin); R->SetBranchAddress("yMax",&yMax); R->SetBranchAddress("dy",&dy);
  R->SetBranchAddress("zMin",&zMin); R->SetBranchAddress("zMax",&zMax); R->SetBranchAddress("dz",&dz);
  R->GetEntry(0);
  auto nOf=[](double lo,double hi,double d){return d>0?(int)std::lround((hi-lo)/d)+1:1;};
  const int nx=nOf(xMin,xMax,dx), ny=nOf(yMin,yMax,dy), nz=nOf(zMin,zMax,dz);
  const long N=(long)nx*ny*nz;
  auto gidx=[ny,nz](int i,int j,int k){return ((long)i*ny+j)*nz+k;};

  std::vector<Acts::Vector3> bField(N, Acts::Vector3::Zero());
  Float_t x,y,z,bx,by,bz;
  D->SetBranchAddress("x",&x); D->SetBranchAddress("y",&y); D->SetBranchAddress("z",&z);
  D->SetBranchAddress("Bx",&bx); D->SetBranchAddress("By",&by); D->SetBranchAddress("Bz",&bz);
  const long nE=D->GetEntries();
  for(long e=0;e<nE;++e){ D->GetEntry(e);
    int i=(dx>0)?(int)std::lround((x-xMin)/dx):0;
    int j=(dy>0)?(int)std::lround((y-yMin)/dy):0;
    int k=(dz>0)?(int)std::lround((z-zMin)/dz):0;
    if(i<0||i>=nx||j<0||j>=ny||k<0||k>=nz) continue;
    bField[gidx(i,j,k)] = Acts::Vector3(bx,by,bz)*bScaleToTesla; }

  std::vector<double> xPos(nx),yPos(ny),zPos(nz);
  for(int i=0;i<nx;++i) xPos[i]=xMin+i*dx+originMm[0];
  for(int j=0;j<ny;++j) yPos[j]=yMin+j*dy+originMm[1];
  for(int k=0;k<nz;++k) zPos[k]=zMin+k*dz+originMm[2];
  auto l2g=[](std::array<std::size_t,3> b,std::array<std::size_t,3> n){
    return b[0]*(n[1]*n[2])+b[1]*n[2]+b[2]; };

  auto map = Acts::fieldMapXYZ(l2g,xPos,yPos,zPos,bField,
                               Acts::UnitConstants::mm, Acts::UnitConstants::T, false);

  // True data extent (mm, ACTS frame), inset by ~one full cell.
  std::array<double,3> lo{ xPos.front() + 1.001*std::abs(dx),
                           yPos.front() + 1.001*std::abs(dy),
                           zPos.front() + 1.001*std::abs(dz) };
  std::array<double,3> hi{ xPos.back()  - 1.001*std::abs(dx),
                           yPos.back()  - 1.001*std::abs(dy),
                           zPos.back()  - 1.001*std::abs(dz) };

  return std::make_shared<ZeroOutsideField<decltype(map)>>(std::move(map), lo, hi);
}

}  // namespace shipreco
