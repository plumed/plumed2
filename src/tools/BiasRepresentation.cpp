/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "BiasRepresentation.h"
#include "core/Value.h"
#include "Communicator.h"
#include <iostream>
#include "KernelFunctions.h"
#include "File.h"
#include "Grid.h"
#include "Tools.h"
#include <algorithm>
#include <cmath>

namespace PLMD {

namespace {
/// factor converting the heights written in a hills file into the bias:
/// with well-tempered metadynamics (biasf>1) the heights are scaled by biasf/(biasf-1),
/// otherwise they are the bias itself
double biasFactorToBias(double biasf) {
  return biasf>1.0 ? (biasf-1.)/biasf : 1.0;
}
}

/// the constructor here
BiasRepresentation::BiasRepresentation(const std::vector<Value*> & tmpvalues, Communicator &cc ):hasgrid(false),rescaledToBias(false),mycomm(cc) {
  lowI_=0.0;
  uppI_=0.0;
  doInt_=false;
  ndim=tmpvalues.size();
  for(int i=0; i<ndim; i++) {
    values.push_back(tmpvalues[i]);
    names.push_back(values[i]->getName());
  }
  resetSupport();
}

/// overload the constructor: add the sigma  at constructor time
BiasRepresentation::BiasRepresentation(const std::vector<Value*> & tmpvalues, Communicator &cc, const std::vector<double> & sigma ):
  hasgrid(false), rescaledToBias(false), histosigma(sigma),mycomm(cc) {
  lowI_=0.0;
  uppI_=0.0;
  doInt_=false;
  ndim=tmpvalues.size();
  for(int i=0; i<ndim; i++) {
    values.push_back(tmpvalues[i]);
    names.push_back(values[i]->getName());
  }
  resetSupport();
}

/// overload the constructor: add the grid at constructor time
BiasRepresentation::BiasRepresentation(const std::vector<Value*> & tmpvalues, Communicator &cc, const std::vector<std::string> & gmin, const std::vector<std::string> & gmax,
                                       const std::vector<unsigned> & nbin, bool doInt, double lowI, double uppI):
  hasgrid(false), rescaledToBias(false), mycomm(cc) {
  ndim=tmpvalues.size();
  for(int i=0; i<ndim; i++) {
    values.push_back(tmpvalues[i]);
    names.push_back(values[i]->getName());
  }
  resetSupport();
  doInt_=doInt;
  lowI_=lowI;
  uppI_=uppI;
  // initialize the grid
  addGrid(gmin,gmax,nbin);
}

/// overload the constructor with some external sigmas: needed for histogram
BiasRepresentation::BiasRepresentation(const std::vector<Value*> & tmpvalues, Communicator &cc, const std::vector<std::string> & gmin, const std::vector<std::string> & gmax,
                                       const std::vector<unsigned> & nbin, const std::vector<double> & sigma):
  hasgrid(false), rescaledToBias(false),histosigma(sigma),mycomm(cc) {
  lowI_=0.0;
  uppI_=0.0;
  doInt_=false;
  ndim=tmpvalues.size();
  for(int  i=0; i<ndim; i++) {
    values.push_back(tmpvalues[i]);
    names.push_back(values[i]->getName());
  }
  resetSupport();
  // initialize the grid
  addGrid(gmin,gmax,nbin);
}

void BiasRepresentation::addGrid(const std::vector<std::string> & gmin, const std::vector<std::string> & gmax, const std::vector<unsigned> & nbin ) {
  plumed_massert(nkernels==0,"you can set the grid before loading the hills");
  plumed_massert(hasgrid==false,"to build the grid you should not having the grid in this bias representation");
  std::string ss;
  ss="file.free";
  std::vector<Value*> vv;
  for(unsigned i=0; i<values.size(); i++) {
    vv.push_back(values[i]);
  }
  BiasGrid_=Tools::make_unique<Grid>(ss,vv,gmin,gmax,nbin,false,true);
  hasgrid=true;
  // geometry of the grid, used in addSeparableGaussian
  gridnbin=BiasGrid_->getNbin();
  gridpbc=BiasGrid_->getIsPeriodic();
  griddx=BiasGrid_->getDx();
  gridmin=BiasGrid_->getPoint(std::vector<unsigned>(ndim,0));
  gridstride.assign(ndim,1);
  for(int d=1; d<ndim; d++) {
    gridstride[d]=gridstride[d-1]*gridnbin[d-1];
  }
  sepIndex.resize(ndim);
  sepT.resize(ndim);
  sepT2.resize(ndim);
  sepExp.resize(ndim);
}

bool BiasRepresentation::hasSigmaInInput() {
  if(histosigma.size()==0) {
    return false;
  } else {
    return true;
  }
}

void BiasRepresentation::setRescaledToBias(bool rescaled) {
  plumed_massert(nkernels==0,"you can set the rescaling function only before loading hills");
  rescaledToBias=rescaled;
}

const bool & BiasRepresentation::isRescaledToBias() {
  return rescaledToBias;
}

unsigned BiasRepresentation::getNumberOfDimensions() {
  return values.size();
}

std::vector<std::string> BiasRepresentation::getNames() {
  return names;
}

const std::string & BiasRepresentation::getName(unsigned i) {
  return names[i];
}

const std::vector<Value*>& BiasRepresentation::getPtrToValues() {
  return values;
}

Value*  BiasRepresentation::getPtrToValue(unsigned i) {
  return values[i];
}

std::unique_ptr<KernelFunctions> BiasRepresentation::readFromPoint(IFile *ifile) {
  std::vector<double> cc( names.size() );
  for(unsigned i=0; i<names.size(); ++i) {
    ifile->scanField(names[i],cc[i]);
  }
  double h=1.0;
  return Tools::make_unique<KernelFunctions>(cc,histosigma,"stretched-gaussian","DIAGONAL",h);
}

void BiasRepresentation::pushKernel( IFile *ifile ) {
  double bf;
  auto kk=readKernel(ifile,bf);
  addKernel(std::move(kk),bf);
}

std::unique_ptr<KernelFunctions> BiasRepresentation::readKernel( IFile *ifile, double & bf ) {
  std::unique_ptr<KernelFunctions> kk;
  // here below the reading of the kernel is completely hidden
  if(histosigma.size()==0) {
    ifile->allowIgnoredFields();
    kk=KernelFunctions::read(ifile,true,names);
  } else {
    // when doing histogram assume gaussian with a given diagonal sigma
    // and neglect all the rest
    kk=readFromPoint(ifile);
  }
  // the bias factor is not something about the kernels but
  // must be stored to keep the  bias/free energy duality
  if(ifile->FieldExist("biasf")) {
    std::string dummy;
    ifile->scanField("biasf",dummy);
    Tools::convert(dummy,bf);
  } else {
    bf=1.0;
  }
  // the domain does not pertain to the kernel but to the values here defined
  std::string mins,maxs,minv,maxv,mini,maxi;
  mins="min_";
  maxs="max_";
  for(int i=0 ; i<ndim; i++) {
    if(values[i]->isPeriodic()) {
      ifile->scanField(mins+names[i],minv);
      ifile->scanField(maxs+names[i],maxv);
      // verify that the domain is correct
      values[i]->getDomain(mini,maxi);
      plumed_massert(mini==minv,"the input periodicity in hills and in value definition does not match"  );
      plumed_massert(maxi==maxv,"the input periodicity in hills and in value definition does not match"  );
    }
  }
  return kk;
}

void BiasRepresentation::addKernel( std::unique_ptr<KernelFunctions> kk, double bf ) {
  nkernels++;
  if(!hasgrid) {
    // without a grid only the support of the kernels is needed (see getMinMaxBin)
    std::vector<double> ss;
    if(histosigma.size()!=0) {
      ss=histosigma;
    } else {
      ss=kk->getContinuousSupport();
    }
    const std::vector<double> cc=kk->getCenter();
    const int ndiv=10; // adjustable parameter: division per support
    for(int j=0; j<ndim; j++) {
      supmin[j]=std::min(supmin[j],cc[j]-ss[j]);
      supmax[j]=std::max(supmax[j],cc[j]+ss[j]);
      minbinsize[j]=std::min(minbinsize[j],ss[j]/double(ndiv));
    }
    return;
  }
  // with MPI each process adds a subset of the kernels to its own grid,
  // the grids are summed in getGridPtr
  if(mycomm.Get_size()>1 && (nkernels-1)%mycomm.Get_size()!=mycomm.Get_rank()) {
    return;
  }
  const double f=(rescaledToBias ? biasFactorToBias(bf) : 1.0);
  double h,a,b,r2cut;
  std::vector<double> w;
  if(!doInt_ && kk->getSeparableGaussian(h,a,b,r2cut,w)) {
    addSeparableGaussian(*kk,f*h,a,b,r2cut,w);
  } else {
    addGenericKernel(*kk,f);
  }
}

void BiasRepresentation::addGenericKernel( const KernelFunctions & kk, double f ) {
  std::vector<unsigned> nneighb;
  if(doInt_&&(kk.getCenter()[0]+kk.getContinuousSupport()[0] > uppI_ || kk.getCenter()[0]-kk.getContinuousSupport()[0] < lowI_ )) {
    nneighb=BiasGrid_->getNbin();
  } else {
    nneighb=kk.getSupport(BiasGrid_->getDx());
  }
  std::vector<Grid::index_t> neighbors=BiasGrid_->getNeighbors(kk.getCenter(),nneighb);
  std::vector<double> der(ndim);
  std::vector<double> xx(ndim);
  for(unsigned i=0; i<neighbors.size(); ++i) {
    Grid::index_t ineigh=neighbors[i];
    BiasGrid_->getPoint(ineigh,xx);
    for(int j=0; j<ndim; ++j) {
      values[j]->set(xx[j]);
    }
    double bias;
    if(doInt_) {
      bias=kk.evaluate(values,der,true,doInt_,lowI_,uppI_);
    } else {
      bias=kk.evaluate(values,der,true);
    }
    bias*=f;
    for(int j=0; j<ndim; ++j) {
      der[j]*=f;
    }
    BiasGrid_->addValueAndDerivatives(ineigh,bias,der);
  }
}

void BiasRepresentation::addSeparableGaussian( const KernelFunctions & kk, double h, double a, double b, double r2cut, const std::vector<double> & w ) {
  // The kernel is evaluated on the same grid points as in addGenericKernel, but the
  // Gaussian factors are computed once per dimension and then multiplied, instead of
  // computing an exponential at each point of the grid.
  const std::vector<double> c=kk.getCenter();
  const std::vector<unsigned> nneigh=kk.getSupport(griddx);
  for(int d=0; d<ndim; d++) {
    auto & idx=sepIndex[d];
    auto & t=sepT[d];
    auto & q=sepT2[d];
    auto & e=sepExp[d];
    idx.clear();
    t.clear();
    q.clear();
    e.clear();
    const int nb=gridnbin[d];
    const int ic=static_cast<int>(std::floor((c[d]-gridmin[d])/griddx[d]));
    int span=2*nneigh[d]+1;
    // along a periodic direction visit each point of the grid only once
    if(gridpbc[d] && span>nb) {
      span=nb;
    }
    for(int k=0; k<span; k++) {
      int i=ic-static_cast<int>(nneigh[d])+k;
      if(!gridpbc[d] && (i<0 || i>=nb)) {
        continue;
      }
      if(gridpbc[d]) {
        i=((i%nb)+nb)%nb;
      }
      const double x=gridmin[d]+static_cast<double>(i)*griddx[d];
      const double tt=-values[d]->difference(x,c[d])/w[d];
      idx.push_back(i);
      t.push_back(tt);
      q.push_back(tt*tt);
      e.push_back(std::exp(-0.5*tt*tt));
    }
    if(idx.empty()) {
      return;
    }
  }
  // prefix (over the dimensions from the last to the second) of r2, of the product of the exponentials and of the index
  std::vector<unsigned> pos(ndim,0);
  std::vector<double> pr2(ndim+1,0.0), pe(ndim+1,1.0);
  std::vector<Grid::index_t> pidx(ndim+1,0);
  auto setLevel=[&](int l) {
    pr2[l]=pr2[l+1]+sepT2[l][pos[l]];
    pe[l]=pe[l+1]*sepExp[l][pos[l]];
    pidx[l]=pidx[l+1]+static_cast<Grid::index_t>(sepIndex[l][pos[l]])*gridstride[l];
  };
  for(int l=ndim-1; l>=1; --l) {
    setLevel(l);
  }
  std::vector<double> der(ndim), dd(ndim);
  const auto & idx0=sepIndex[0];
  const auto & t0=sepT[0];
  const auto & q0=sepT2[0];
  const auto & e0=sepExp[0];
  const double ha=h*a;
  const double hb=h*b;
  while(true) {
    const double r2u=pr2[1];
    const double eu=pe[1];
    const Grid::index_t iu=pidx[1];
    for(int l=1; l<ndim; ++l) {
      der[l]=sepT[l][pos[l]]/w[l];
    }
    for(unsigned k=0; k<idx0.size(); ++k) {
      if(r2u+q0[k]>=r2cut) {
        continue;
      }
      const double ee=eu*e0[k];
      const double kderiv=-ha*ee;
      dd[0]=t0[k]/w[0]*kderiv;
      for(int l=1; l<ndim; ++l) {
        dd[l]=der[l]*kderiv;
      }
      BiasGrid_->addValueAndDerivatives(iu+idx0[k],ha*ee+hb,dd.data(),ndim);
    }
    // next combination of the indices of the other dimensions
    int l=1;
    while(l<ndim) {
      pos[l]++;
      if(pos[l]<sepIndex[l].size()) {
        break;
      }
      pos[l]=0;
      l++;
    }
    if(l>=ndim) {
      break;
    }
    for(int m=l; m>=1; --m) {
      setLevel(m);
    }
  }
}

int BiasRepresentation::getNumberOfKernels() {
  return nkernels;
}

Grid* BiasRepresentation::getGridPtr() {
  plumed_massert(hasgrid,"if you want the grid pointer then you should have defined a grid before");
  if(mycomm.Get_size()==1) {
    return BiasGrid_.get();
  }
  // sum the grids of all the processes
  const auto n=BiasGrid_->getSize();
  std::vector<double> val(n), der(n*ndim);
  for(Grid::index_t i=0; i<n; i++) {
    val[i]=BiasGrid_->getValueAndDerivatives(i,&der[i*ndim],ndim);
  }
  mycomm.Sum(val);
  mycomm.Sum(der);
  reducedGrid_=Tools::make_unique<Grid>(*BiasGrid_);
  std::vector<double> dd(ndim);
  for(Grid::index_t i=0; i<n; i++) {
    for(int j=0; j<ndim; j++) {
      dd[j]=der[i*ndim+j];
    }
    reducedGrid_->setValueAndDerivatives(i,val[i],dd);
  }
  return reducedGrid_.get();
}

void BiasRepresentation::getMinMaxBin(std::vector<double> &vmin, std::vector<double> &vmax, std::vector<unsigned> &vbin) {
  vmin=supmin;
  vmax=supmax;
  vbin.assign(ndim,0);
  for(int j=0; j<ndim; j++) {
    // reset to periodicity
    if(values[j]->isPeriodic()) {
      double minv,maxv;
      values[j]->getDomain(minv,maxv);
      if(minv>vmin[j]) {
        vmin[j]=minv;
      }
      if(maxv<vmax[j]) {
        vmax[j]=maxv;
      }
    }
    vbin[j]=static_cast<unsigned>(std::ceil((vmax[j]-vmin[j])/minbinsize[j]) );
  }
}

void BiasRepresentation::clear() {
  nkernels=0;
  resetSupport();
  // clear the grid
  if(hasgrid) {
    BiasGrid_->clear();
  }
}

void BiasRepresentation::resetSupport() {
  supmin.assign(ndim,10.e20);
  supmax.assign(ndim,-10.e20);
  minbinsize.assign(ndim,10.e20);
}


}
