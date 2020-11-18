/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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

namespace PLMD {

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
}

/// overload the constructor: add the sigma  at constructor time
BiasRepresentation::BiasRepresentation(const std::vector<Value*> & tmpvalues, Communicator &cc, const std::vector<double> & sigma ):
  hasgrid(false), rescaledToBias(false), histosigma(sigma),mycomm(cc)
{
  lowI_=0.0;
  uppI_=0.0;
  doInt_=false;
  ndim=tmpvalues.size();
  for(int i=0; i<ndim; i++) {
    values.push_back(tmpvalues[i]);
    names.push_back(values[i]->getName());
  }
}

/// overload the constructor: add the grid at constructor time
BiasRepresentation::BiasRepresentation(const std::vector<Value*> & tmpvalues, Communicator &cc, const std::vector<std::string> & gmin, const std::vector<std::string> & gmax,
                                       const std::vector<unsigned> & nbin, bool doInt, double lowI, double uppI):
  hasgrid(false), rescaledToBias(false), mycomm(cc)
{
  ndim=tmpvalues.size();
  for(int i=0; i<ndim; i++) {
    values.push_back(tmpvalues[i]);
    names.push_back(values[i]->getName());
  }
  doInt_=doInt;
  lowI_=lowI;
  uppI_=uppI;
  // initialize the grid
  addGrid(gmin,gmax,nbin);
}

/// overload the constructor with some external sigmas: needed for histogram
BiasRepresentation::BiasRepresentation(const std::vector<Value*> & tmpvalues, Communicator &cc, const std::vector<std::string> & gmin, const std::vector<std::string> & gmax,
                                       const std::vector<unsigned> & nbin, const std::vector<double> & sigma):
  hasgrid(false), rescaledToBias(false),histosigma(sigma),mycomm(cc)
{
  lowI_=0.0;
  uppI_=0.0;
  doInt_=false;
  ndim=tmpvalues.size();
  for(int  i=0; i<ndim; i++) {
    values.push_back(tmpvalues[i]);
    names.push_back(values[i]->getName());
  }
  // initialize the grid
  addGrid(gmin,gmax,nbin);
}

void BiasRepresentation::addGrid(const std::vector<std::string> & gmin, const std::vector<std::string> & gmax, const std::vector<unsigned> & nbin ) {
  plumed_massert(hills.size()==0,"you can set the grid before loading the hills");
  plumed_massert(hasgrid==false,"to build the grid you should not having the grid in this bias representation");
  std::string ss; ss="file.free";
  std::vector<Value*> vv; for(unsigned i=0; i<values.size(); i++) vv.push_back(values[i]);
  BiasGrid_=Tools::make_unique<Grid>(ss,vv,gmin,gmax,nbin,false,true);
  hasgrid=true;
}

bool BiasRepresentation::hasSigmaInInput() {
  if(histosigma.size()==0) {return false;} else {return true;}
}

void BiasRepresentation::setRescaledToBias(bool rescaled) {
  plumed_massert(hills.size()==0,"you can set the rescaling function only before loading hills");
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
  return Tools::make_unique<KernelFunctions>(cc,histosigma,"gaussian","DIAGONAL",h);
}

void BiasRepresentation::pushKernel( IFile *ifile ) {
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
  std::string dummy; double dummyd;
  if(ifile->FieldExist("biasf")) {
    ifile->scanField("biasf",dummy);
    Tools::convert(dummy,dummyd);
  } else {dummyd=1.0;}
  biasf.push_back(dummyd);
  // the domain does not pertain to the kernel but to the values here defined
  std::string mins,maxs,minv,maxv,mini,maxi; mins="min_"; maxs="max_";
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
  // if grid is defined then it should be added on the grid
  //cerr<<"now with "<<hills.size()<<endl;
  if(hasgrid) {
    std::vector<unsigned> nneighb;
    if(doInt_&&(kk->getCenter()[0]+kk->getContinuousSupport()[0] > uppI_ || kk->getCenter()[0]-kk->getContinuousSupport()[0] < lowI_ )) {
      nneighb=BiasGrid_->getNbin();
    } else nneighb=kk->getSupport(BiasGrid_->getDx());
    std::vector<Grid::index_t> neighbors=BiasGrid_->getNeighbors(kk->getCenter(),nneighb);
    std::vector<double> der(ndim);
    std::vector<double> xx(ndim);
    if(mycomm.Get_size()==1) {
      for(unsigned i=0; i<neighbors.size(); ++i) {
        Grid::index_t ineigh=neighbors[i];
        for(int j=0; j<ndim; ++j) {der[j]=0.0;}
        BiasGrid_->getPoint(ineigh,xx);
        // assign xx to a new vector of values
        for(int j=0; j<ndim; ++j) {values[j]->set(xx[j]);}
        double bias;
        if(doInt_) bias=kk->evaluate(values,der,true,doInt_,lowI_,uppI_);
        else bias=kk->evaluate(values,der,true);
        if(rescaledToBias) {
          double f=(biasf.back()-1.)/(biasf.back());
          bias*=f;
          for(int j=0; j<ndim; ++j) {der[j]*=f;}
        }
        BiasGrid_->addValueAndDerivatives(ineigh,bias,der);
      }
    } else {
      unsigned stride=mycomm.Get_size();
      unsigned rank=mycomm.Get_rank();
      std::vector<double> allder(ndim*neighbors.size(),0.0);
      std::vector<double> allbias(neighbors.size(),0.0);
      std::vector<double> tmpder(ndim);
      for(unsigned i=rank; i<neighbors.size(); i+=stride) {
        Grid::index_t ineigh=neighbors[i];
        BiasGrid_->getPoint(ineigh,xx);
        for(int j=0; j<ndim; ++j) {values[j]->set(xx[j]);}
        if(doInt_) allbias[i]=kk->evaluate(values,der,true,doInt_,lowI_,uppI_);
        else allbias[i]=kk->evaluate(values,der,true);
        if(rescaledToBias) {
          double f=(biasf.back()-1.)/(biasf.back());
          allbias[i]*=f;
          for(int j=0; j<ndim; ++j) {tmpder[j]*=f;}
        }
        // this solution with the temporary vector is rather bad, probably better to take
        // a pointer of double as it was in old gaussian
        for(int j=0; j<ndim; ++j) { allder[ndim*i+j]=tmpder[j]; tmpder[j]=0.;}
      }
      mycomm.Sum(allbias);
      mycomm.Sum(allder);
      for(unsigned i=0; i<neighbors.size(); ++i) {
        Grid::index_t ineigh=neighbors[i];
        for(int j=0; j<ndim; ++j) {der[j]=allder[ndim*i+j];}
        BiasGrid_->addValueAndDerivatives(ineigh,allbias[i],der);
      }
    }
  }
  hills.emplace_back(std::move(kk));
}

int BiasRepresentation::getNumberOfKernels() {
  return hills.size();
}

Grid* BiasRepresentation::getGridPtr() {
  plumed_massert(hasgrid,"if you want the grid pointer then you should have defined a grid before");
  return BiasGrid_.get();
}

void BiasRepresentation::getMinMaxBin(std::vector<double> &vmin, std::vector<double> &vmax, std::vector<unsigned> &vbin) {
  std::vector<double> ss,cc,binsize;
  vmin.clear(); vmin.resize(ndim,10.e20);
  vmax.clear(); vmax.resize(ndim,-10.e20);
  vbin.clear(); vbin.resize(ndim);
  binsize.clear(); binsize.resize(ndim,10.e20);
  int ndiv=10; // adjustable parameter: division per support
  for(unsigned i=0; i<hills.size(); i++) {
    if(histosigma.size()!=0) {
      ss=histosigma;
    } else {
      ss=hills[i]->getContinuousSupport();
    }
    cc=hills[i]->getCenter();
    for(int j=0; j<ndim; j++) {
      double dmin=cc[j]-ss[j];
      double dmax=cc[j]+ss[j];
      double ddiv=ss[j]/double(ndiv);
      if(dmin<vmin[j])vmin[j]=dmin;
      if(dmax>vmax[j])vmax[j]=dmax;
      if(ddiv<binsize[j])binsize[j]=ddiv;
    }
  }
  for(int j=0; j<ndim; j++) {
    // reset to periodicity
    if(values[j]->isPeriodic()) {
      double minv,maxv;
      values[j]->getDomain(minv,maxv);
      if(minv>vmin[j])vmin[j]=minv;
      if(maxv<vmax[j])vmax[j]=maxv;
    }
    vbin[j]=static_cast<unsigned>(std::ceil((vmax[j]-vmin[j])/binsize[j]) );
  }
}

void BiasRepresentation::clear() {
  hills.clear();
  // clear the grid
  if(hasgrid) {
    BiasGrid_->clear();
  }
}


}
