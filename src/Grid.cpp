#include <vector>
#include <cassert>
#include <cmath>
#include "Grid.h"
#include <iostream>

using namespace std;
using namespace PLMD;

Grid::Grid(vector<double> gmin, vector<double> gmax, vector<unsigned> nbin,
           vector<bool> pbc, bool dospline, bool usederiv, bool doclear){
// various checks
 assert(gmax.size()==gmin.size());
 assert(gmax.size()==nbin.size());
 assert(gmax.size()==pbc.size());
 for(unsigned int i=0;i<gmax.size();++i){
  assert(gmax[i]>gmin[i]);
  assert(nbin[i]>0);
 }
 dimension_=gmax.size();
 min_=gmin;
 max_=gmax;
 nbin_=nbin;
 pbc_=pbc;
 dospline_=dospline;
 usederiv_=usederiv;
 if(dospline_) assert(dospline_==usederiv_);
 maxsize_=1;
 for(unsigned int i=0;i<dimension_;++i){
  dx_.push_back((max_[i]-min_[i])/(double)nbin_[i]);
  maxsize_*=nbin_[i];
 }
 if(doclear) clear();
}

void Grid::clear(){
 grid_.resize(maxsize_);
 if(usederiv_) der_.resize(maxsize_);
 for(unsigned int i=0;i<maxsize_;++i){
  grid_[i]=0.0;
  if(usederiv_){
   (der_[i]).resize(dimension_); 
   for(unsigned int j=0;j<dimension_;++j) der_[i][j]=0.0;
  }
 }
}

vector<double> Grid::getMin() const {
 return min_;
}

vector<double> Grid::getMax() const {
 return max_;
}

vector<double> Grid::getDx() const {
 return dx_;
}

vector<unsigned> Grid::getNbin() const {
 return nbin_;
}

unsigned Grid::getSize() const {
 return maxsize_;
}

unsigned Grid::getDimension() const {
 return dimension_;
}

// we are flattening arrays using a column-major order
unsigned Grid::getIndex(vector<unsigned> indices) const {
 assert(indices.size()==dimension_);
 unsigned index=indices[dimension_-1];
 for(unsigned int i=dimension_-1;i>0;--i){
  index=index*nbin_[i-1]+indices[i-1];
 }
 return index;
}

unsigned Grid::getIndex(vector<double> x) const {
 assert(x.size()==dimension_);
 return getIndex(getIndices(x));
}

// we are flattening arrays using a column-major order
vector<unsigned> Grid::getIndices(unsigned index) const {
 vector<unsigned> indices;
 unsigned kk=index;
 indices.push_back(index%nbin_[0]);
 for(unsigned int i=1;i<dimension_-1;++i){
  kk=(kk-indices[i-1])/nbin_[i-1];
  indices.push_back(kk%nbin_[i]);
 }
 if(dimension_>=2){
  indices.push_back((kk-indices[dimension_-2])/nbin_[dimension_-2]);
 }
 return indices;
}

vector<unsigned> Grid::getIndices(vector<double> x) const {
 assert(x.size()==dimension_);
 vector<unsigned> indices;
 for(unsigned int i=0;i<dimension_;++i){
   indices.push_back(floor((x[i]-min_[i])/dx_[i]));
 }
 return indices;
}

vector<double> Grid::getPoint(vector<unsigned> indices) const {
 assert(indices.size()==dimension_);
 vector<double> x;
 for(unsigned int i=0;i<dimension_;++i){
  x.push_back(min_[i]+(double)(indices[i])*dx_[i]);
 }
 return x;
}

vector<double> Grid::getPoint(unsigned index) const {
 assert(index>=0 && index<maxsize_);
 return getPoint(getIndices(index));
}

vector<double> Grid::getPoint(vector<double> x) const {
 assert(x.size()==dimension_);
 return getPoint(getIndices(x));
}

vector<unsigned> Grid::getNeighbors
 (vector<unsigned> indices,vector<unsigned> nneigh){
 assert(indices.size()==dimension_ && nneigh.size()==dimension_);

 vector<unsigned> neighbors, small_bin; 
 unsigned small_nbin=1;
 for(unsigned j=0;j<dimension_;++j){
  small_bin.push_back(2*nneigh[j]+1);
  small_nbin*=small_bin[j];
 }
 
 for(unsigned index=0;index<small_nbin;++index){
  vector<unsigned> small_indices;
  unsigned kk=index;
  small_indices.push_back(index%small_bin[0]);
  for(unsigned i=1;i<dimension_-1;++i){
   kk=(kk-small_indices[i-1])/small_bin[i-1];
   small_indices.push_back(kk%small_bin[i]);
  }
  if(dimension_>=2){
   small_indices.push_back((kk-small_indices[dimension_-2])/small_bin[dimension_-2]);
  }
  vector<unsigned> tmp_indices;
  for(unsigned i=0;i<dimension_;++i){
   int i0=small_indices[i]-nneigh[i]+indices[i];
   if(!pbc_[i] && i0<0)         continue;
   if(!pbc_[i] && i0>=nbin_[i]) continue;
   if( pbc_[i] && i0<0)         i0+=nbin_[i];
   if( pbc_[i] && i0>=nbin_[i]) i0-=nbin_[i];
   tmp_indices.push_back((unsigned)i0);
  }
  if(tmp_indices.size()==dimension_){neighbors.push_back(getIndex(tmp_indices));}
 } 
 return neighbors;
}
 
vector<unsigned> Grid::getNeighbors
 (vector<double> x,vector<unsigned> nneigh){
 assert(x.size()==dimension_ && nneigh.size()==dimension_);
 return getNeighbors(getIndices(x),nneigh);
}

vector<unsigned> Grid::getNeighbors
 (unsigned index,vector<unsigned> nneigh){
 assert(index>=0 && index<maxsize_ && nneigh.size()==dimension_);
 return getNeighbors(getIndices(index),nneigh);
}

vector<unsigned> Grid::getSplineNeighbors(vector<unsigned> indices){
 assert(indices.size()==dimension_);
 vector<unsigned> neighbors;
 unsigned nneigh=pow(2.0,int(dimension_));
 
 for(unsigned int i=0;i<nneigh;++i){
  unsigned tmp=i;
  vector<unsigned> nindices;
  for(unsigned int j=0;j<dimension_;++j){
   unsigned i0=tmp%2+indices[j];
   tmp/=2;
   if(!pbc_[j] && i0==nbin_[j]) continue;
   if( pbc_[j] && i0==nbin_[j]) i0=0;
   nindices.push_back(i0);
  }
  if(nindices.size()==dimension_) neighbors.push_back(getIndex(nindices));
 }
 return neighbors;
}
 
// all the methods below should be overwritten by a class 
// that inherits from Grid
double Grid::getValue(unsigned index) {
 assert(index>=0 && index<maxsize_);
 return grid_[index];
}

double Grid::getValue(vector<unsigned> indices) {
 return getValue(getIndex(indices));
}

double Grid::getValue(vector<double> x) {
 if(!dospline_){
  return getValue(getIndex(x));
 } else {
  vector<double> der;
  der.resize(dimension_);
  return getValueAndDerivatives(x,der);
 }
}

double Grid::getValueAndDerivatives
 (unsigned index, vector<double>& der){
 assert(index>=0 && index<maxsize_ && usederiv_ && der.size()==dimension_);
 der=der_[index];
 return grid_[index];
}

double Grid::getValueAndDerivatives
 (vector<unsigned> indices, vector<double>& der){
 return getValueAndDerivatives(getIndex(indices),der);
}

double Grid::getValueAndDerivatives
(vector<double> x, vector<double>& der) {
 assert(der.size()==dimension_ && usederiv_);
 
 if(dospline_){
  double X,X2,X3,value;
  double* fd=new double[dimension_];
  double* C=new double[dimension_];
  double* D=new double[dimension_];
  vector<double> dder;
  dder.resize(dimension_);
// reset
  value=0.0;
  for(unsigned int i=0;i<dimension_;++i) der[i]=0.0;

  vector<unsigned> indices=getIndices(x);
  vector<unsigned> neigh=getSplineNeighbors(indices);
  vector<double>   xfloor=getPoint(x);

// loop over neighbors
  for(unsigned int ipoint=0;ipoint<neigh.size();++ipoint){
   double grid=getValueAndDerivatives(neigh[ipoint],dder);
   vector<unsigned> nindices=getIndices(neigh[ipoint]);
   double ff=1.0;

   for(unsigned j=0;j<dimension_;++j){
    int x0=1;
    if(nindices[j]==indices[j]) x0=0;
    double dx=getDx()[j];
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
   for(unsigned j=0;j<dimension_;++j){
    fd[j]=D[j];
    for(unsigned i=0;i<dimension_;++i) if(i!=j) fd[j]*=C[i];
   }
   value+=grid*ff;
   for(unsigned j=0;j<dimension_;++j) der[j]+=grid*fd[j];
  }
  delete(fd); delete(C); delete(D);
  return value;
 }else{
  return getValueAndDerivatives(getIndex(x),der);
 }
}

void Grid::setValue(unsigned index, double value){
 assert(index>=0 && index<maxsize_ && !usederiv_);
 grid_[index]=value;
}

void Grid::setValue(vector<unsigned> indices, double value){
 setValue(getIndex(indices),value); 
}

void Grid::setValueAndDerivatives
 (unsigned index, double value, vector<double>& der){
 assert(index>=0 && index<maxsize_ && usederiv_ && der.size()==dimension_);
 grid_[index]=value;
 der_[index]=der;
}

void Grid::setValueAndDerivatives
 (vector<unsigned> indices, double value, vector<double>& der){
 setValueAndDerivatives(getIndex(indices),value,der); 
}

void Grid::addValue(unsigned index, double value){
 assert(index>=0 && index<maxsize_ && !usederiv_);
 grid_[index]+=value;
}

void Grid::addValue(vector<unsigned> indices, double value){
 addValue(getIndex(indices),value);
}

void Grid::addValueAndDerivatives
 (unsigned index, double value, vector<double>& der){
 assert(index>=0 && index<maxsize_ && usederiv_ && der.size()==dimension_);
 grid_[index]+=value;
 for(unsigned int i=0;i<dimension_;++i) der_[index][i]+=der[i];
}

void Grid::addValueAndDerivatives
 (vector<unsigned> indices, double value, vector<double>& der){
 addValueAndDerivatives(getIndex(indices),value,der);
}

// Sparse version of grid with map
void SparseGrid::clear(){
 map_.clear();
}

unsigned SparseGrid::getSize() const{
 return maxsize_;
}

double SparseGrid::getUsedSize() const {
 return map_.size();
}

double SparseGrid::getValue(unsigned index) {
 assert(index>=0 && index<maxsize_);
 double value=0.0;
 it_=map_.find(index);
 if(it_!=map_.end()) value=it_->second;
 return value;
}

double SparseGrid::getValue(vector<unsigned> indices) {
 return getValue(getIndex(indices));
}

double SparseGrid::getValue(vector<double> x) {
 if(!dospline_){
  return getValue(getIndex(x));
 } else {
  vector<double> der;
  der.resize(dimension_);
  return getValueAndDerivatives(x,der);
 }
}

double SparseGrid::getValueAndDerivatives
 (unsigned index, vector<double>& der){
 assert(index>=0 && index<maxsize_ && usederiv_ && der.size()==dimension_);
 double value=0.0;
 for(unsigned int i=0;i<dimension_;++i) der[i]=0.0;
 it_=map_.find(index);
 if(it_!=map_.end()) value=it_->second;
 itder_=der_.find(index);
 if(itder_!=der_.end()) der=itder_->second;
 return value;
}

double SparseGrid::getValueAndDerivatives
 (vector<unsigned> indices, vector<double>& der){
 return getValueAndDerivatives(getIndex(indices),der);
}

double SparseGrid::getValueAndDerivatives
 (vector<double> x, vector<double>& der){
 assert(der.size()==dimension_ && usederiv_);


if(dospline_){
  double X,X2,X3,value;
  double* fd=new double[dimension_];
  double* C=new double[dimension_];
  double* D=new double[dimension_];
  vector<double> dder;
  dder.resize(dimension_);
// reset
  value=0.0;
  for(unsigned int i=0;i<dimension_;++i) der[i]=0.0;

  vector<unsigned> indices=getIndices(x);
  vector<unsigned> neigh=getSplineNeighbors(indices);
  vector<double>   xfloor=getPoint(x);

// loop over neighbors
  for(unsigned int ipoint=0;ipoint<neigh.size();++ipoint){
   double grid=getValueAndDerivatives(neigh[ipoint],dder);
   vector<unsigned> nindices=getIndices(neigh[ipoint]);
   double ff=1.0;

   for(unsigned j=0;j<dimension_;++j){
    int x0=1;
    if(nindices[j]==indices[j]) x0=0;
    double dx=getDx()[j];
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
   for(unsigned j=0;j<dimension_;++j){
    fd[j]=D[j];
    for(unsigned i=0;i<dimension_;++i) if(i!=j) fd[j]*=C[i];
   }
   value+=grid*ff;
   for(unsigned j=0;j<dimension_;++j) der[j]+=grid*fd[j];
  }
  delete(fd); delete(C); delete(D);
  return value;
 }else{
  return getValueAndDerivatives(getIndex(x),der);
 }
}

void SparseGrid::setValue(unsigned index, double value){
 assert(index>=0 && index<maxsize_ && !usederiv_);
 map_[index]=value;
}

void SparseGrid::setValue(vector<unsigned> indices, double value){
 setValue(getIndex(indices),value); 
}

void SparseGrid::setValueAndDerivatives
 (unsigned index, double value, vector<double>& der){
 assert(index>=0 && index<maxsize_ && usederiv_ && der.size()==dimension_);
 map_[index]=value;
 der_[index]=der;
}

void SparseGrid::setValueAndDerivatives
(vector<unsigned> indices, double value, vector<double>& der){
 setValueAndDerivatives(getIndex(indices),value,der); 
}

void SparseGrid::addValue(unsigned index, double value){
 assert(index>=0 && index<maxsize_ && !usederiv_);
 map_[index]+=value;
}

void SparseGrid::addValue(vector<unsigned> indices, double value){
 addValue(getIndex(indices),value);
}

void SparseGrid::addValueAndDerivatives
 (unsigned index, double value, vector<double>& der){
 assert(index>=0 && index<maxsize_ && usederiv_ && der.size()==dimension_);
 map_[index]+=value;
 der_[index].resize(dimension_);
 for(unsigned int i=0;i<dimension_;++i) der_[index][i]+=der[i]; 
}

void SparseGrid::addValueAndDerivatives
(vector<unsigned> indices, double value, vector<double>& der){
 addValueAndDerivatives(getIndex(indices),value,der);
}
