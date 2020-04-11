/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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

#include "Grid.h"
#include "Tools.h"
#include "core/Value.h"
#include "File.h"
#include "Exception.h"
#include "KernelFunctions.h"
#include "RootFindingBase.h"
#include "Communicator.h"

#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cfloat>
#include <array>

using namespace std;
namespace PLMD {

constexpr size_t GridBase::maxdim;

GridBase::GridBase(const std::string& funcl, const std::vector<Value*> & args, const vector<std::string> & gmin,
                   const vector<std::string> & gmax, const vector<unsigned> & nbin, bool dospline, bool usederiv) {
// various checks
  plumed_assert(args.size()<=maxdim) << "grid dim cannot exceed "<<maxdim;
  plumed_massert(args.size()==gmin.size(),"grid min dimensions in input do not match number of arguments");
  plumed_massert(args.size()==nbin.size(),"number of bins on input do not match number of arguments");
  plumed_massert(args.size()==gmax.size(),"grid max dimensions in input do not match number of arguments");
  unsigned dim=gmax.size();
  std::vector<std::string> names;
  std::vector<bool> isperiodic;
  std::vector<string> pmin,pmax;
  names.resize( dim );
  isperiodic.resize( dim );
  pmin.resize( dim );
  pmax.resize( dim );
  for(unsigned int i=0; i<dim; ++i) {
    names[i]=args[i]->getName();
    if( args[i]->isPeriodic() ) {
      isperiodic[i]=true;
      args[i]->getDomain( pmin[i], pmax[i] );
    } else {
      isperiodic[i]=false;
      pmin[i]="0.";
      pmax[i]="0.";
    }
  }
// this is a value-independent initializator
  Init(funcl,names,gmin,gmax,nbin,dospline,usederiv,isperiodic,pmin,pmax);
}

GridBase::GridBase(const std::string& funcl, const std::vector<string> &names, const std::vector<std::string> & gmin,
                   const vector<std::string> & gmax, const std::vector<unsigned> & nbin, bool dospline, bool usederiv, const std::vector<bool> &isperiodic, const std::vector<std::string> &pmin, const std::vector<std::string> &pmax ) {
// this calls the initializator
  Init(funcl,names,gmin,gmax,nbin,dospline,usederiv,isperiodic,pmin,pmax);
}

void GridBase::Init(const std::string& funcl, const std::vector<std::string> &names, const vector<std::string> & gmin,
                    const std::vector<std::string> & gmax, const std::vector<unsigned> & nbin, bool dospline, bool usederiv,
                    const std::vector<bool> &isperiodic, const std::vector<std::string> &pmin, const std::vector<std::string> &pmax ) {
  fmt_="%14.9f";
// various checks
  plumed_assert(names.size()<=maxdim) << "grid size cannot exceed "<<maxdim;
  plumed_massert(names.size()==gmin.size(),"grid dimensions in input do not match number of arguments");
  plumed_massert(names.size()==nbin.size(),"grid dimensions in input do not match number of arguments");
  plumed_massert(names.size()==gmax.size(),"grid dimensions in input do not match number of arguments");
  dimension_=gmax.size();
  str_min_=gmin; str_max_=gmax;
  argnames.resize( dimension_ );
  min_.resize( dimension_ );
  max_.resize( dimension_ );
  pbc_.resize( dimension_ );
  for(unsigned int i=0; i<dimension_; ++i) {
    argnames[i]=names[i];
    if( isperiodic[i] ) {
      pbc_[i]=true;
      str_min_[i]=pmin[i];
      str_max_[i]=pmax[i];
    } else {
      pbc_[i]=false;
    }
    Tools::convert(str_min_[i],min_[i]);
    Tools::convert(str_max_[i],max_[i]);
    funcname=funcl;
    plumed_massert(max_[i]>min_[i],"maximum in grid must be larger than minimum");
    plumed_massert(nbin[i]>0,"number of grid points must be greater than zero");
  }
  nbin_=nbin;
  dospline_=dospline;
  usederiv_=usederiv;
  if(dospline_) plumed_assert(dospline_==usederiv_);
  maxsize_=1;
  for(unsigned int i=0; i<dimension_; ++i) {
    dx_.push_back( (max_[i]-min_[i])/static_cast<double>( nbin_[i] ) );
    if( !pbc_[i] ) { max_[i] += dx_[i]; nbin_[i] += 1; }
    maxsize_*=nbin_[i];
  }
}

vector<std::string> GridBase::getMin() const {
  return str_min_;
}

vector<std::string> GridBase::getMax() const {
  return str_max_;
}

vector<double> GridBase::getDx() const {
  return dx_;
}

double GridBase::getDx(index_t j) const {
  return dx_[j];
}

double GridBase::getBinVolume() const {
  double vol=1.;
  for(unsigned i=0; i<dx_.size(); ++i) vol*=dx_[i];
  return vol;
}

vector<bool> GridBase::getIsPeriodic() const {
  return pbc_;
}

vector<unsigned> GridBase::getNbin() const {
  return nbin_;
}

vector<string> GridBase::getArgNames() const {
  return argnames;
}


unsigned GridBase::getDimension() const {
  return dimension_;
}

// we are flattening arrays using a column-major order
GridBase::index_t GridBase::getIndex(const vector<unsigned> & indices) const {
  plumed_dbg_assert(indices.size()==dimension_);
  for(unsigned int i=0; i<dimension_; i++)
    if(indices[i]>=nbin_[i]) {
      std::string is;
      Tools::convert(i,is);
      std::string msg="ERROR: the system is looking for a value outside the grid along the " + is + " ("+getArgNames()[i]+")";
      plumed_merror(msg+" index!");
    }
  index_t index=indices[dimension_-1];
  for(unsigned int i=dimension_-1; i>0; --i) {
    index=index*nbin_[i-1]+indices[i-1];
  }
  return index;
}

GridBase::index_t GridBase::getIndex(const vector<double> & x) const {
  plumed_dbg_assert(x.size()==dimension_);
  return getIndex(getIndices(x));
}

// we are flattening arrays using a column-major order
vector<unsigned> GridBase::getIndices(index_t index) const {
  vector<unsigned> indices(dimension_);
  index_t kk=index;
  indices[0]=(index%nbin_[0]);
  for(unsigned int i=1; i<dimension_-1; ++i) {
    kk=(kk-indices[i-1])/nbin_[i-1];
    indices[i]=(kk%nbin_[i]);
  }
  if(dimension_>=2) {
    indices[dimension_-1]=((kk-indices[dimension_-2])/nbin_[dimension_-2]);
  }
  return indices;
}

void GridBase::getIndices(index_t index, std::vector<unsigned>& indices) const {
  if (indices.size()!=dimension_) indices.resize(dimension_);
  index_t kk=index;
  indices[0]=(index%nbin_[0]);
  for(unsigned int i=1; i<dimension_-1; ++i) {
    kk=(kk-indices[i-1])/nbin_[i-1];
    indices[i]=(kk%nbin_[i]);
  }
  if(dimension_>=2) {
    indices[dimension_-1]=((kk-indices[dimension_-2])/nbin_[dimension_-2]);
  }
}

vector<unsigned> GridBase::getIndices(const vector<double> & x) const {
  plumed_dbg_assert(x.size()==dimension_);
  vector<unsigned> indices(dimension_);
  for(unsigned int i=0; i<dimension_; ++i) {
    indices[i] = unsigned(floor((x[i]-min_[i])/dx_[i]));
  }
  return indices;
}

void GridBase::getIndices(const vector<double> & x, std::vector<unsigned>& indices) const {
  plumed_dbg_assert(x.size()==dimension_);
  if (indices.size()!=dimension_) indices.resize(dimension_);
  for(unsigned int i=0; i<dimension_; ++i) {
    indices[i] = unsigned(floor((x[i]-min_[i])/dx_[i]));
  }
}

vector<double> GridBase::getPoint(const vector<unsigned> & indices) const {
  plumed_dbg_assert(indices.size()==dimension_);
  vector<double> x(dimension_);
  for(unsigned int i=0; i<dimension_; ++i) {
    x[i]=min_[i]+(double)(indices[i])*dx_[i];
  }
  return x;
}

vector<double> GridBase::getPoint(index_t index) const {
  plumed_dbg_assert(index<maxsize_);
  return getPoint(getIndices(index));
}

vector<double> GridBase::getPoint(const vector<double> & x) const {
  plumed_dbg_assert(x.size()==dimension_);
  return getPoint(getIndices(x));
}

void GridBase::getPoint(index_t index,std::vector<double> & point) const {
  plumed_dbg_assert(index<maxsize_);
  getPoint(getIndices(index),point);
}

void GridBase::getPoint(const std::vector<unsigned> & indices,std::vector<double> & point) const {
  plumed_dbg_assert(indices.size()==dimension_);
  plumed_dbg_assert(point.size()==dimension_);
  for(unsigned int i=0; i<dimension_; ++i) {
    point[i]=min_[i]+(double)(indices[i])*dx_[i];
  }
}

void GridBase::getPoint(const std::vector<double> & x,std::vector<double> & point) const {
  plumed_dbg_assert(x.size()==dimension_);
  getPoint(getIndices(x),point);
}

vector<GridBase::index_t> GridBase::getNeighbors
(const vector<unsigned> &indices,const vector<unsigned> &nneigh)const {
  plumed_dbg_assert(indices.size()==dimension_ && nneigh.size()==dimension_);

  vector<index_t> neighbors;
  vector<unsigned> small_bin(dimension_);

  unsigned small_nbin=1;
  for(unsigned j=0; j<dimension_; ++j) {
    small_bin[j]=(2*nneigh[j]+1);
    small_nbin*=small_bin[j];
  }

  vector<unsigned> small_indices(dimension_);
  vector<unsigned> tmp_indices;
  for(unsigned index=0; index<small_nbin; ++index) {
    tmp_indices.resize(dimension_);
    unsigned kk=index;
    small_indices[0]=(index%small_bin[0]);
    for(unsigned i=1; i<dimension_-1; ++i) {
      kk=(kk-small_indices[i-1])/small_bin[i-1];
      small_indices[i]=(kk%small_bin[i]);
    }
    if(dimension_>=2) {
      small_indices[dimension_-1]=((kk-small_indices[dimension_-2])/small_bin[dimension_-2]);
    }
    unsigned ll=0;
    for(unsigned i=0; i<dimension_; ++i) {
      int i0=small_indices[i]-nneigh[i]+indices[i];
      if(!pbc_[i] && i0<0)         continue;
      if(!pbc_[i] && i0>=static_cast<int>(nbin_[i])) continue;
      if( pbc_[i] && i0<0)         i0=nbin_[i]-(-i0)%nbin_[i];
      if( pbc_[i] && i0>=static_cast<int>(nbin_[i])) i0%=nbin_[i];
      tmp_indices[ll]=static_cast<unsigned>(i0);
      ll++;
    }
    tmp_indices.resize(ll);
    if(tmp_indices.size()==dimension_) {neighbors.push_back(getIndex(tmp_indices));}
  }
  return neighbors;
}

vector<GridBase::index_t> GridBase::getNeighbors
(const vector<double> & x,const vector<unsigned> & nneigh)const {
  plumed_dbg_assert(x.size()==dimension_ && nneigh.size()==dimension_);
  return getNeighbors(getIndices(x),nneigh);
}

vector<GridBase::index_t> GridBase::getNeighbors
(index_t index,const vector<unsigned> & nneigh)const {
  plumed_dbg_assert(index<maxsize_ && nneigh.size()==dimension_);
  return getNeighbors(getIndices(index),nneigh);
}

void GridBase::getSplineNeighbors(const vector<unsigned> & indices, vector<GridBase::index_t>& neighbors, unsigned& nneighbors)const {
  plumed_dbg_assert(indices.size()==dimension_);
  unsigned nneigh=unsigned(pow(2.0,int(dimension_)));
  if (neighbors.size()!=nneigh) neighbors.resize(nneigh);

  vector<unsigned> nindices(dimension_);
  unsigned inind; nneighbors = 0;
  for(unsigned int i=0; i<nneigh; ++i) {
    unsigned tmp=i; inind=0;
    for(unsigned int j=0; j<dimension_; ++j) {
      unsigned i0=tmp%2+indices[j];
      tmp/=2;
      if(!pbc_[j] && i0==nbin_[j]) continue;
      if( pbc_[j] && i0==nbin_[j]) i0=0;
      nindices[inind++]=i0;
    }
    if(inind==dimension_) neighbors[nneighbors++]=getIndex(nindices);
  }
}

vector<GridBase::index_t> GridBase::getNearestNeighbors(const index_t index) const {
  vector<index_t> nearest_neighs = vector<index_t>();
  for (unsigned i = 0; i < dimension_; i++) {
    vector<unsigned> neighsneeded = vector<unsigned>(dimension_, 0);
    neighsneeded[i] = 1;
    vector<index_t> singledim_nearest_neighs = getNeighbors(index, neighsneeded);
    for (unsigned j = 0; j < singledim_nearest_neighs.size(); j++) {
      index_t neigh = singledim_nearest_neighs[j];
      if (neigh != index) {
        nearest_neighs.push_back(neigh);
      }
    }
  }
  return nearest_neighs;
}

vector<GridBase::index_t> GridBase::getNearestNeighbors(const vector<unsigned> &indices) const {
  plumed_dbg_assert(indices.size() == dimension_);
  return getNearestNeighbors(getIndex(indices));
}


void GridBase::addKernel( const KernelFunctions& kernel ) {
  plumed_dbg_assert( kernel.ndim()==dimension_ );
  std::vector<unsigned> nneighb=kernel.getSupport( dx_ );
  std::vector<index_t> neighbors=getNeighbors( kernel.getCenter(), nneighb );
  std::vector<double> xx( dimension_ );
  std::vector<std::unique_ptr<Value>> vv( dimension_ );
  std::string str_min, str_max;
  for(unsigned i=0; i<dimension_; ++i) {
    vv[i].reset(new Value());
    if( pbc_[i] ) {
      Tools::convert(min_[i],str_min);
      Tools::convert(max_[i],str_max);
      vv[i]->setDomain( str_min, str_max );
    } else {
      vv[i]->setNotPeriodic();
    }
  }

// vv_ptr contains plain pointers obtained from vv.
// this is the simplest way to replace a unique_ptr here.
// perhaps the interface of kernel.evaluate() should be changed
// in order to accept a std::vector<std::unique_ptr<Value>>
  auto vv_ptr=Tools::unique2raw(vv);

  std::vector<double> der( dimension_ );
  for(unsigned i=0; i<neighbors.size(); ++i) {
    index_t ineigh=neighbors[i];
    getPoint( ineigh, xx );
    for(unsigned j=0; j<dimension_; ++j) vv[j]->set(xx[j]);
    double newval = kernel.evaluate( vv_ptr, der, usederiv_ );
    if( usederiv_ ) addValueAndDerivatives( ineigh, newval, der );
    else addValue( ineigh, newval );
  }
}


double GridBase::getValue(const vector<unsigned> & indices) const {
  return getValue(getIndex(indices));
}

double GridBase::getValue(const vector<double> & x) const {
  if(!dospline_) {
    return getValue(getIndex(x));
  } else {
    vector<double> der(dimension_);
    return getValueAndDerivatives(x,der);
  }
}

double GridBase::getValueAndDerivatives
(const vector<unsigned> & indices, vector<double>& der) const {
  return getValueAndDerivatives(getIndex(indices),der);
}

double GridBase::getValueAndDerivatives
(const vector<double> & x, vector<double>& der) const {
  plumed_dbg_assert(der.size()==dimension_ && usederiv_);

  if(dospline_) {
    double X,X2,X3,value;
    std::array<double,maxdim> fd, C, D;
    std::vector<double> dder(dimension_);
// reset
    value=0.0;
    for(unsigned int i=0; i<dimension_; ++i) der[i]=0.0;

    vector<unsigned> indices(dimension_);
    getIndices(x, indices);
    vector<double> xfloor(dimension_);
    getPoint(indices, xfloor);
    vector<index_t> neigh; unsigned nneigh; getSplineNeighbors(indices, neigh, nneigh);

// loop over neighbors
    vector<unsigned> nindices;
    for(unsigned int ipoint=0; ipoint<nneigh; ++ipoint) {
      double grid=getValueAndDerivatives(neigh[ipoint],dder);
      getIndices(neigh[ipoint], nindices);
      double ff=1.0;

      for(unsigned j=0; j<dimension_; ++j) {
        int x0=1;
        if(nindices[j]==indices[j]) x0=0;
        double dx=getDx(j);
        X=fabs((x[j]-xfloor[j])/dx-(double)x0);
        X2=X*X;
        X3=X2*X;
        double yy;
        if(fabs(grid)<0.0000001) yy=0.0;
        else yy=-dder[j]/grid;
        C[j]=(1.0-3.0*X2+2.0*X3) - (x0?-1.0:1.0)*yy*(X-2.0*X2+X3)*dx;
        D[j]=( -6.0*X +6.0*X2) - (x0?-1.0:1.0)*yy*(1.0-4.0*X +3.0*X2)*dx;
        D[j]*=(x0?-1.0:1.0)/dx;
        ff*=C[j];
      }
      for(unsigned j=0; j<dimension_; ++j) {
        fd[j]=D[j];
        for(unsigned i=0; i<dimension_; ++i) if(i!=j) fd[j]*=C[i];
      }
      value+=grid*ff;
      for(unsigned j=0; j<dimension_; ++j) der[j]+=grid*fd[j];
    }
    return value;
  } else {
    return getValueAndDerivatives(getIndex(x),der);
  }
}

void GridBase::setValue(const vector<unsigned> & indices, double value) {
  setValue(getIndex(indices),value);
}

void GridBase::setValueAndDerivatives
(const vector<unsigned> & indices, double value, vector<double>& der) {
  setValueAndDerivatives(getIndex(indices),value,der);
}

void GridBase::addValue(const vector<unsigned> & indices, double value) {
  addValue(getIndex(indices),value);
}

void GridBase::addValueAndDerivatives
(const vector<unsigned> & indices, double value, vector<double>& der) {
  addValueAndDerivatives(getIndex(indices),value,der);
}

void GridBase::writeHeader(OFile& ofile) {
  for(unsigned i=0; i<dimension_; ++i) {
    ofile.addConstantField("min_" + argnames[i]);
    ofile.addConstantField("max_" + argnames[i]);
    ofile.addConstantField("nbins_" + argnames[i]);
    ofile.addConstantField("periodic_" + argnames[i]);
  }
}

void Grid::clear() {
  grid_.assign(maxsize_,0.0);
  if(usederiv_) der_.assign(maxsize_*dimension_,0.0);
}

void Grid::writeToFile(OFile& ofile) {
  vector<double> xx(dimension_);
  vector<double> der(dimension_);
  double f;
  writeHeader(ofile);
  for(index_t i=0; i<getSize(); ++i) {
    xx=getPoint(i);
    if(usederiv_) {f=getValueAndDerivatives(i,der);}
    else {f=getValue(i);}
    if(i>0 && dimension_>1 && getIndices(i)[dimension_-2]==0) ofile.printf("\n");
    for(unsigned j=0; j<dimension_; ++j) {
      ofile.printField("min_" + argnames[j], str_min_[j] );
      ofile.printField("max_" + argnames[j], str_max_[j] );
      ofile.printField("nbins_" + argnames[j], static_cast<int>(nbin_[j]) );
      if( pbc_[j] ) ofile.printField("periodic_" + argnames[j], "true" );
      else          ofile.printField("periodic_" + argnames[j], "false" );
    }
    for(unsigned j=0; j<dimension_; ++j) { ofile.fmtField(" "+fmt_); ofile.printField(argnames[j],xx[j]); }
    ofile.fmtField(" "+fmt_); ofile.printField(funcname,f);
    if(usederiv_) for(unsigned j=0; j<dimension_; ++j) { ofile.fmtField(" "+fmt_); ofile.printField("der_" + argnames[j],der[j]); }
    ofile.printField();
  }
}

void GridBase::writeCubeFile(OFile& ofile, const double& lunit) {
  plumed_assert( dimension_==3 );
  ofile.printf("PLUMED CUBE FILE\n");
  ofile.printf("OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n");
  // Number of atoms followed by position of origin (origin set so that center of grid is in center of cell)
  ofile.printf("%d %f %f %f\n",1,-0.5*lunit*(max_[0]-min_[0]),-0.5*lunit*(max_[1]-min_[1]),-0.5*lunit*(max_[2]-min_[2]));
  ofile.printf("%u %f %f %f\n",nbin_[0],lunit*dx_[0],0.0,0.0);  // Number of bins in each direction followed by
  ofile.printf("%u %f %f %f\n",nbin_[1],0.0,lunit*dx_[1],0.0);  // shape of voxel
  ofile.printf("%u %f %f %f\n",nbin_[2],0.0,0.0,lunit*dx_[2]);
  ofile.printf("%d %f %f %f\n",1,0.0,0.0,0.0); // Fake atom otherwise VMD doesn't work
  std::vector<unsigned> pp(3);
  for(pp[0]=0; pp[0]<nbin_[0]; ++pp[0]) {
    for(pp[1]=0; pp[1]<nbin_[1]; ++pp[1]) {
      for(pp[2]=0; pp[2]<nbin_[2]; ++pp[2]) {
        ofile.printf("%f ",getValue(pp) );
        if(pp[2]%6==5) ofile.printf("\n");
      }
      ofile.printf("\n");
    }
  }
}

std::unique_ptr<GridBase> GridBase::create(const std::string& funcl, const std::vector<Value*> & args, IFile& ifile,
    const vector<std::string> & gmin,const vector<std::string> & gmax,
    const vector<unsigned> & nbin,bool dosparse, bool dospline, bool doder) {
  std::unique_ptr<GridBase> grid=GridBase::create(funcl,args,ifile,dosparse,dospline,doder);
  std::vector<unsigned> cbin( grid->getNbin() );
  std::vector<std::string> cmin( grid->getMin() ), cmax( grid->getMax() );
  for(unsigned i=0; i<args.size(); ++i) {
    plumed_massert( cmin[i]==gmin[i], "mismatched grid min" );
    plumed_massert( cmax[i]==gmax[i], "mismatched grid max" );
    if( args[i]->isPeriodic() ) {
      plumed_massert( cbin[i]==nbin[i], "mismatched grid nbins" );
    } else {
      plumed_massert( (cbin[i]-1)==nbin[i], "mismatched grid nbins");
    }
  }
  return grid;
}

std::unique_ptr<GridBase> GridBase::create(const std::string& funcl, const std::vector<Value*> & args, IFile& ifile, bool dosparse, bool dospline, bool doder)
{
  std::unique_ptr<GridBase> grid;
  unsigned nvar=args.size(); bool hasder=false; std::string pstring;
  std::vector<int> gbin1(nvar); std::vector<unsigned> gbin(nvar);
  std::vector<std::string> labels(nvar),gmin(nvar),gmax(nvar);
  std::vector<std::string> fieldnames; ifile.scanFieldList( fieldnames );
// Retrieve names for fields
  for(unsigned i=0; i<args.size(); ++i) labels[i]=args[i]->getName();
// And read the stuff from the header
  plumed_massert( ifile.FieldExist( funcl ), "no column labelled " + funcl + " in in grid input");
  for(unsigned i=0; i<args.size(); ++i) {
    ifile.scanField( "min_" + labels[i], gmin[i]);
    ifile.scanField( "max_" + labels[i], gmax[i]);
    ifile.scanField( "periodic_" + labels[i], pstring );
    ifile.scanField( "nbins_" + labels[i], gbin1[i]);
    plumed_assert( gbin1[i]>0 );
    if( args[i]->isPeriodic() ) {
      plumed_massert( pstring=="true", "input value is periodic but grid is not");
      std::string pmin, pmax;
      args[i]->getDomain( pmin, pmax ); gbin[i]=gbin1[i];
      if( pmin!=gmin[i] || pmax!=gmax[i] ) plumed_merror("mismatch between grid boundaries and periods of values");
    } else {
      gbin[i]=gbin1[i]-1;  // Note header in grid file indicates one more bin that there should be when data is not periodic
      plumed_massert( pstring=="false", "input value is not periodic but grid is");
    }
    hasder=ifile.FieldExist( "der_" + args[i]->getName() );
    if( doder && !hasder ) plumed_merror("missing derivatives from grid file");
    for(unsigned j=0; j<fieldnames.size(); ++j) {
      for(unsigned k=i+1; k<args.size(); ++k) {
        if( fieldnames[j]==labels[k] ) plumed_merror("arguments in input are not in same order as in grid file");
      }
      if( fieldnames[j]==labels[i] ) break;
    }
  }

  if(!dosparse) {grid.reset(new Grid(funcl,args,gmin,gmax,gbin,dospline,doder));}
  else {grid.reset(new SparseGrid(funcl,args,gmin,gmax,gbin,dospline,doder));}

  vector<double> xx(nvar),dder(nvar);
  vector<double> dx=grid->getDx();
  double f,x;
  while( ifile.scanField(funcl,f) ) {
    for(unsigned i=0; i<nvar; ++i) {
      ifile.scanField(labels[i],x); xx[i]=x+dx[i]/2.0;
      ifile.scanField( "min_" + labels[i], gmin[i]);
      ifile.scanField( "max_" + labels[i], gmax[i]);
      ifile.scanField( "nbins_" + labels[i], gbin1[i]);
      ifile.scanField( "periodic_" + labels[i], pstring );
    }
    if(hasder) { for(unsigned i=0; i<nvar; ++i) { ifile.scanField( "der_" + args[i]->getName(), dder[i] ); } }
    index_t index=grid->getIndex(xx);
    if(doder) {grid->setValueAndDerivatives(index,f,dder);}
    else {grid->setValue(index,f);}
    ifile.scanField();
  }
  return grid;
}

double Grid::getMinValue() const {
  double minval;
  minval=DBL_MAX;
  for(index_t i=0; i<grid_.size(); ++i) {
    if(grid_[i]<minval)minval=grid_[i];
  }
  return minval;
}

double Grid::getMaxValue() const {
  double maxval;
  maxval=DBL_MIN;
  for(index_t i=0; i<grid_.size(); ++i) {
    if(grid_[i]>maxval)maxval=grid_[i];
  }
  return maxval;
}

void Grid::scaleAllValuesAndDerivatives( const double& scalef ) {
  if(usederiv_) {
    for(index_t i=0; i<grid_.size(); ++i) {
      grid_[i]*=scalef;
      for(unsigned j=0; j<dimension_; ++j) der_[i*dimension_+j]*=scalef;
    }
  } else {
    for(index_t i=0; i<grid_.size(); ++i) grid_[i]*=scalef;
  }
}

void Grid::logAllValuesAndDerivatives( const double& scalef ) {
  if(usederiv_) {
    for(index_t i=0; i<grid_.size(); ++i) {
      grid_[i] = scalef*log(grid_[i]);
      for(unsigned j=0; j<dimension_; ++j) der_[i*dimension_+j] = scalef/der_[i*dimension_+j];
    }
  } else {
    for(index_t i=0; i<grid_.size(); ++i) grid_[i] = scalef*log(grid_[i]);
  }
}

void Grid::setMinToZero() {
  double min=grid_[0];
  for(index_t i=1; i<grid_.size(); ++i) if(grid_[i]<min) min=grid_[i];
  for(index_t i=0; i<grid_.size(); ++i) grid_[i] -= min;
}

void Grid::applyFunctionAllValuesAndDerivatives( double (*func)(double val), double (*funcder)(double valder) ) {
  if(usederiv_) {
    for(index_t i=0; i<grid_.size(); ++i) {
      grid_[i]=func(grid_[i]);
      for(unsigned j=0; j<dimension_; ++j) der_[i*dimension_+j]=funcder(der_[i*dimension_+j]);
    }
  } else {
    for(index_t i=0; i<grid_.size(); ++i) grid_[i]=func(grid_[i]);
  }
}

double Grid::getDifferenceFromContour( const std::vector<double>& x, std::vector<double>& der ) const {
  return getValueAndDerivatives( x, der ) - contour_location;
}

void Grid::findSetOfPointsOnContour(const double& target, const std::vector<bool>& nosearch,
                                    unsigned& npoints, std::vector<std::vector<double> >& points ) {
// Set contour location for function
  contour_location=target;
// Resize points to maximum possible value
  points.resize( dimension_*maxsize_ );

// Two points for search
  std::vector<unsigned> ind(dimension_);
  std::vector<double> direction( dimension_, 0 );

// Run over whole grid
  npoints=0; RootFindingBase<Grid> mymin( this );
  for(unsigned i=0; i<maxsize_; ++i) {
    for(unsigned j=0; j<dimension_; ++j) ind[j]=getIndices(i)[j];

    // Get the value of a point on the grid
    double val1=getValue(i) - target;

    // Now search for contour in each direction
    bool edge=false;
    for(unsigned j=0; j<dimension_; ++j) {
      if( nosearch[j] ) continue ;
      // Make sure we don't search at the edge of the grid
      if( !pbc_[j] && (ind[j]+1)==nbin_[j] ) continue;
      else if( (ind[j]+1)==nbin_[j] ) { edge=true; ind[j]=0; }
      else ind[j]+=1;
      double val2=getValue(ind) - target;
      if( val1*val2<0 ) {
        // Use initial point location as first guess for search
        points[npoints].resize(dimension_); for(unsigned k=0; k<dimension_; ++k) points[npoints][k]=getPoint(i)[k];
        // Setup direction vector
        direction[j]=0.999999999*dx_[j];
        // And do proper search for contour point
        mymin.linesearch( direction, points[npoints], &Grid::getDifferenceFromContour );
        direction[j]=0.0; npoints++;
      }
      if( pbc_[j] && edge ) { edge=false; ind[j]=nbin_[j]-1; }
      else ind[j]-=1;
    }
  }
}

/// OVERRIDES ARE BELOW

Grid::index_t Grid::getSize() const {
  return maxsize_;
}

double Grid::getValue(index_t index) const {
  plumed_dbg_assert(index<maxsize_);
  return grid_[index];
}

double Grid::getValueAndDerivatives
(index_t index, vector<double>& der) const {
  plumed_dbg_assert(index<maxsize_ && usederiv_ && der.size()==dimension_);
  der.resize(dimension_);
  for(unsigned i=0; i<dimension_; i++) der[i]=der_[dimension_*index+i];
  return grid_[index];
}

void Grid::setValue(index_t index, double value) {
  plumed_dbg_assert(index<maxsize_ && !usederiv_);
  grid_[index]=value;
}

void Grid::setValueAndDerivatives
(index_t index, double value, vector<double>& der) {
  plumed_dbg_assert(index<maxsize_ && usederiv_ && der.size()==dimension_);
  grid_[index]=value;
  for(unsigned i=0; i<dimension_; i++) der_[dimension_*index+i]=der[i];
}

void Grid::addValue(index_t index, double value) {
  plumed_dbg_assert(index<maxsize_ && !usederiv_);
  grid_[index]+=value;
}

void Grid::addValueAndDerivatives
(index_t index, double value, vector<double>& der) {
  plumed_dbg_assert(index<maxsize_ && usederiv_ && der.size()==dimension_);
  grid_[index]+=value;
  for(unsigned int i=0; i<dimension_; ++i) der_[index*dimension_+i]+=der[i];
}

Grid::index_t SparseGrid::getSize() const {
  return map_.size();
}

Grid::index_t SparseGrid::getMaxSize() const {
  return maxsize_;
}

double SparseGrid::getValue(index_t index)const {
  plumed_assert(index<maxsize_);
  double value=0.0;
  const auto it=map_.find(index);
  if(it!=map_.end()) value=it->second;
  return value;
}

double SparseGrid::getValueAndDerivatives
(index_t index, vector<double>& der)const {
  plumed_assert(index<maxsize_ && usederiv_ && der.size()==dimension_);
  double value=0.0;
  for(unsigned int i=0; i<dimension_; ++i) der[i]=0.0;
  const auto it=map_.find(index);
  if(it!=map_.end()) value=it->second;
  const auto itder=der_.find(index);
  if(itder!=der_.end()) der=itder->second;
  return value;
}

void SparseGrid::setValue(index_t index, double value) {
  plumed_assert(index<maxsize_ && !usederiv_);
  map_[index]=value;
}

void SparseGrid::setValueAndDerivatives
(index_t index, double value, vector<double>& der) {
  plumed_assert(index<maxsize_ && usederiv_ && der.size()==dimension_);
  map_[index]=value;
  der_[index]=der;
}

void SparseGrid::addValue(index_t index, double value) {
  plumed_assert(index<maxsize_ && !usederiv_);
  map_[index]+=value;
}

void SparseGrid::addValueAndDerivatives
(index_t index, double value, vector<double>& der) {
  plumed_assert(index<maxsize_ && usederiv_ && der.size()==dimension_);
  map_[index]+=value;
  der_[index].resize(dimension_);
  for(unsigned int i=0; i<dimension_; ++i) der_[index][i]+=der[i];
}

void SparseGrid::writeToFile(OFile& ofile) {
  vector<double> xx(dimension_);
  vector<double> der(dimension_);
  double f;
  writeHeader(ofile);
  ofile.fmtField(" "+fmt_);
  for(const auto & it : map_) {
    index_t i=it.first;
    xx=getPoint(i);
    if(usederiv_) {f=getValueAndDerivatives(i,der);}
    else {f=getValue(i);}
    if(i>0 && dimension_>1 && getIndices(i)[dimension_-2]==0) ofile.printf("\n");
    for(unsigned j=0; j<dimension_; ++j) {
      ofile.printField("min_" + argnames[j], str_min_[j] );
      ofile.printField("max_" + argnames[j], str_max_[j] );
      ofile.printField("nbins_" + argnames[j], static_cast<int>(nbin_[j]) );
      if( pbc_[j] ) ofile.printField("periodic_" + argnames[j], "true" );
      else          ofile.printField("periodic_" + argnames[j], "false" );
    }
    for(unsigned j=0; j<dimension_; ++j) ofile.printField(argnames[j],xx[j]);
    ofile.printField(funcname, f);
    if(usederiv_) { for(unsigned j=0; j<dimension_; ++j) ofile.printField("der_" + argnames[j],der[j]); }
    ofile.printField();
  }
}

double SparseGrid::getMinValue() const {
  double minval;
  minval=0.0;
  for(auto const & i : map_) {
    if(i.second<minval) minval=i.second;
  }
  return minval;
}

double SparseGrid::getMaxValue() const {
  double maxval;
  maxval=0.0;
  for(auto const & i : map_) {
    if(i.second>maxval) maxval=i.second;
  }
  return maxval;
}

void Grid::projectOnLowDimension(double &val, std::vector<int> &vHigh, WeightBase * ptr2obj ) {
  unsigned i=0;
  for(i=0; i<vHigh.size(); i++) {
    if(vHigh[i]<0) { // this bin needs to be integrated out
      // parallelize here???
      for(unsigned j=0; j<(getNbin())[i]; j++) {
        vHigh[i]=int(j);
        projectOnLowDimension(val,vHigh,ptr2obj); // recursive function: this is the core of the mechanism
        vHigh[i]=-1;
      }
      return; //
    }
  }
  // when there are no more bin to dig in then retrieve the value
  if(i==vHigh.size()) {
    //std::cerr<<"POINT: ";
    //for(unsigned j=0;j<vHigh.size();j++){
    //   std::cerr<<vHigh[j]<<" ";
    //}
    std::vector<unsigned> vv(vHigh.size());
    for(unsigned j=0; j<vHigh.size(); j++)vv[j]=unsigned(vHigh[j]);
    //
    // this is the real assignment !!!!! (hack this to have bias or other stuff)
    //

    // this case: produce fes
    //val+=exp(beta*getValue(vv)) ;
    double myv=getValue(vv);
    val=ptr2obj->projectInnerLoop(val,myv) ;
    // to be added: bias (same as before without negative sign)
    //std::cerr<<" VAL: "<<val <<endl;
  }
}

Grid Grid::project(const std::vector<std::string> & proj, WeightBase *ptr2obj ) {
  // find extrema only for the projection
  vector<string>   smallMin,smallMax;
  vector<unsigned> smallBin;
  vector<unsigned> dimMapping;
  vector<bool> smallIsPeriodic;
  vector<string> smallName;

  // check if the two key methods are there
  WeightBase* pp = dynamic_cast<WeightBase*>(ptr2obj);
  if (!pp)plumed_merror("This WeightBase is not complete: you need a projectInnerLoop and projectOuterLoop ");

  for(unsigned j=0; j<proj.size(); j++) {
    for(unsigned i=0; i<getArgNames().size(); i++) {
      if(proj[j]==getArgNames()[i]) {
        unsigned offset;
        // note that at sizetime the non periodic dimension get a bin more
        if(getIsPeriodic()[i]) {offset=0;} else {offset=1;}
        smallMax.push_back(getMax()[i]);
        smallMin.push_back(getMin()[i]);
        smallBin.push_back(getNbin()[i]-offset);
        smallIsPeriodic.push_back(getIsPeriodic()[i]);
        dimMapping.push_back(i);
        smallName.push_back(getArgNames()[i]);
        break;
      }
    }
  }
  Grid smallgrid("projection",smallName,smallMin,smallMax,smallBin,false,false,smallIsPeriodic,smallMin,smallMax);
  // check that the two grids are commensurate
  for(unsigned i=0; i<dimMapping.size(); i++) {
    plumed_massert(  (smallgrid.getMax())[i] == (getMax())[dimMapping[i]],  "the two input grids are not compatible in max"   );
    plumed_massert(  (smallgrid.getMin())[i] == (getMin())[dimMapping[i]],  "the two input grids are not compatible in min"   );
    plumed_massert(  (smallgrid.getNbin())[i]== (getNbin())[dimMapping[i]], "the two input grids are not compatible in bin"   );
  }
  vector<unsigned> toBeIntegrated;
  for(unsigned i=0; i<getArgNames().size(); i++) {
    bool doappend=true;
    for(unsigned j=0; j<dimMapping.size(); j++) {
      if(dimMapping[j]==i) {doappend=false; break;}
    }
    if(doappend)toBeIntegrated.push_back(i);
  }
  //for(unsigned i=0;i<dimMapping.size();i++ ){
  //     cerr<<"Dimension to preserve "<<dimMapping[i]<<endl;
  //}
  //for(unsigned i=0;i<toBeIntegrated.size();i++ ){
  //     cerr<<"Dimension to integrate "<<toBeIntegrated[i]<<endl;
  //}

  // loop over all the points in the Grid, find the corresponding fixed index, rotate over all the other ones
  for(unsigned i=0; i<smallgrid.getSize(); i++) {
    std::vector<unsigned> v;
    v=smallgrid.getIndices(i);
    std::vector<int> vHigh((getArgNames()).size(),-1);
    for(unsigned j=0; j<dimMapping.size(); j++)vHigh[dimMapping[j]]=int(v[j]);
    // the vector vhigh now contains at the beginning the index of the low dimension and -1 for the values that need to be integrated
    double initval=0.;
    projectOnLowDimension(initval,vHigh, ptr2obj);
    smallgrid.setValue(i,initval);
  }
  // reset to zero just for biasing (this option can be evtl enabled in a future...)
  //double vmin;vmin=-smallgrid.getMinValue()+1;
  for(unsigned i=0; i<smallgrid.getSize(); i++) {
    //         //if(dynamic_cast<BiasWeight*>(ptr2obj)){
    //         //        smallgrid.addValue(i,vmin);// go to 1
    //         //}
    double vv=smallgrid.getValue(i);
    smallgrid.setValue(i,ptr2obj->projectOuterLoop(vv));
    //         //if(dynamic_cast<BiasWeight*>(ptr2obj)){
    //         //        smallgrid.addValue(i,-vmin);// bring back to the value
    //         //}
  }

  return smallgrid;
}

double Grid::integrate( std::vector<unsigned>& npoints ) {
  plumed_dbg_assert( npoints.size()==dimension_ ); plumed_assert( dospline_ );

  unsigned ntotgrid=1; double box_vol=1.0;
  std::vector<double> ispacing( npoints.size() );
  for(unsigned j=0; j<dimension_; ++j) {
    if( !pbc_[j] ) {
      ispacing[j] = ( max_[j] - dx_[j] - min_[j] ) / static_cast<double>( npoints[j] );
      npoints[j]+=1;
    } else {
      ispacing[j] = ( max_[j] - min_[j] ) / static_cast<double>( npoints[j] );
    }
    ntotgrid*=npoints[j]; box_vol*=ispacing[j];
  }

  std::vector<double> vals( dimension_ );
  std::vector<unsigned> t_index( dimension_ ); double integral=0.0;
  for(unsigned i=0; i<ntotgrid; ++i) {
    t_index[0]=(i%npoints[0]);
    unsigned kk=i;
    for(unsigned j=1; j<dimension_-1; ++j) { kk=(kk-t_index[j-1])/npoints[j-1]; t_index[j]=(kk%npoints[j]); }
    if( dimension_>=2 ) t_index[dimension_-1]=((kk-t_index[dimension_-2])/npoints[dimension_-2]);

    for(unsigned j=0; j<dimension_; ++j) vals[j]=min_[j] + t_index[j]*ispacing[j];

    integral += getValue( vals );
  }

  return box_vol*integral;
}

void Grid::mpiSumValuesAndDerivatives( Communicator& comm ) {
  comm.Sum( grid_ ); for(unsigned i=0; i<der_.size(); ++i) comm.Sum( der_[i] );
}


bool indexed_lt(pair<Grid::index_t, double> const &x, pair<Grid::index_t, double> const   &y) {
  return x.second < y.second;
}

double GridBase::findMaximalPathMinimum(const std::vector<double> &source, const std::vector<double> &sink) {
  plumed_dbg_assert(source.size() == dimension_);
  plumed_dbg_assert(sink.size() == dimension_);
  // Start and end indices
  index_t source_idx = getIndex(source);
  index_t sink_idx = getIndex(sink);
  // Path cost
  double maximal_minimum = 0;
  // In one dimension, path searching is very easy--either go one way if it's not periodic,
  // or go both ways if it is periodic. There's no reason to pay the cost of Dijkstra.
  if (dimension_ == 1) {
    // Do a search from the grid source to grid sink that does not
    // cross the grid boundary.
    double curr_min_bias = getValue(source_idx);
    // Either search from a high source to a low sink.
    if (source_idx > sink_idx) {
      for (index_t i = source_idx; i >= sink_idx; i--) {
        if (curr_min_bias == 0.0) {
          break;
        }
        curr_min_bias = fmin(curr_min_bias, getValue(i));
      }
      // Or search from a low source to a high sink.
    } else if (source_idx < sink_idx) {
      for (index_t i = source_idx; i <= sink_idx; i++) {
        if (curr_min_bias == 0.0) {
          break;
        }
        curr_min_bias = fmin(curr_min_bias, getValue(i));
      }
    }
    maximal_minimum = curr_min_bias;
    // If the grid is periodic, also do the search that crosses
    // the grid boundary.
    if (pbc_[0]) {
      double curr_min_bias = getValue(source_idx);
      // Either go from a high source to the upper boundary and
      // then from the bottom boundary to the sink
      if (source_idx > sink_idx) {
        for (index_t i = source_idx; i < maxsize_; i++) {
          if (curr_min_bias == 0.0) {
            break;
          }
          curr_min_bias = fmin(curr_min_bias, getValue(i));
        }
        for (index_t i = 0; i <= sink_idx; i++) {
          if (curr_min_bias == 0.0) {
            break;
          }
          curr_min_bias = fmin(curr_min_bias, getValue(i));
        }
        // Or go from a low source to the bottom boundary and
        // then from the high boundary to the sink
      } else if (source_idx < sink_idx) {
        for (index_t i = source_idx; i > 0; i--) {
          if (curr_min_bias == 0.0) {
            break;
          }
          curr_min_bias = fmin(curr_min_bias, getValue(i));
        }
        curr_min_bias = fmin(curr_min_bias, getValue(0));
        for (index_t i = maxsize_ - 1; i <= sink_idx; i--) {
          if (curr_min_bias == 0.0) {
            break;
          }
          curr_min_bias = fmin(curr_min_bias, getValue(i));
        }
      }
      // If the boundary crossing paths was more biased, it's
      // minimal bias replaces the non-boundary-crossing path's
      // minimum.
      maximal_minimum = fmax(maximal_minimum, curr_min_bias);
    }
    // The one dimensional path search is complete.
    return maximal_minimum;
    // In two or more dimensions, path searching isn't trivial and we really
    // do need to use a path search algorithm. Dijkstra is the simplest decent
    // one. Using it we've never found the path search to be performance
    // limiting in any solvated biomolecule test system, but faster options are
    // easy to imagine if they become necessary. NB-In this case, we're actually
    // using a greedy variant of Dijkstra's algorithm where the first possible
    // path to a point always controls the path cost to that point. The structure
    // of the cost function in this case guarantees that the calculated costs will
    // be correct using this variant even though fine details of the paths may not
    // match a normal Dijkstra search.
  } else if (dimension_ > 1) {
    // Prepare calculation temporaries for Dijkstra's algorithm.
    // Minimal path costs from source to a given grid point
    vector<double> mins_from_source = vector<double>(maxsize_, -1.0);
    // Heap for tracking available steps, steps are recorded as std::pairs of
    // an index and a value.
    vector< pair<index_t, double> > next_steps;
    pair<index_t, double> curr_indexed_val;
    make_heap(next_steps.begin(), next_steps.end(), indexed_lt);
    // The search begins at the source index.
    next_steps.push_back(pair<index_t, double>(source_idx, getValue(source_idx)));
    push_heap(next_steps.begin(), next_steps.end(), indexed_lt);
    // At first no points have been examined and the optimal path has not been found.
    index_t n_examined = 0;
    bool path_not_found = true;
    // Until a path is found,
    while (path_not_found) {
      // Examine the grid point currently most accessible from
      // the set of all previously explored grid points by popping
      // it from the top of the heap.
      pop_heap(next_steps.begin(), next_steps.end(), indexed_lt);
      curr_indexed_val = next_steps.back();
      next_steps.pop_back();
      n_examined++;
      // Check if this point is the sink point, and if so
      // finish the loop.
      if (curr_indexed_val.first == sink_idx) {
        path_not_found = false;
        maximal_minimum = curr_indexed_val.second;
        break;
        // Check if this point has reached the worst possible
        // value, and if so stop looking for paths.
      } else if (curr_indexed_val.second == 0.0) {
        maximal_minimum = 0.0;
        break;
      }
      // If the search is not over, add this grid point's neighbors to the
      // possible next points to search for the sink.
      vector<index_t> neighs = getNearestNeighbors(curr_indexed_val.first);
      for (unsigned k = 0; k < neighs.size(); k++) {
        index_t i = neighs[k];
        // If the neighbor has not already been added to the list of possible next steps,
        if (mins_from_source[i] == -1.0) {
          // Set the cost to reach it via a path through the current point being examined.
          mins_from_source[i] = fmin(curr_indexed_val.second, getValue(i));
          // Add the neighboring point to the heap of potential next steps.
          next_steps.push_back(pair<index_t, double>(i, mins_from_source[i]));
          push_heap(next_steps.begin(), next_steps.end(), indexed_lt);
        }
      }
      // Move on to the next best looking step along any of the paths
      // growing from the source.
    }
    // The multidimensional path search is now complete.
    return maximal_minimum;
  }
  return 0.0;
}

}
