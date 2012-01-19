#include <vector>
#include <cmath>
#include <iostream>
#include <cstdio>

#include "Grid.h"
#include "PlumedException.h"

using namespace std;
using namespace PLMD;

Grid::Grid(const vector<double> & gmin, const vector<double> & gmax, const vector<unsigned> & nbin,
           const vector<bool> & pbc, bool dospline, bool usederiv, bool doclear){
// various checks
 plumed_assert(gmax.size()==gmin.size());
 plumed_assert(gmax.size()==nbin.size());
 plumed_assert(gmax.size()==pbc.size());
 for(unsigned int i=0;i<gmax.size();++i){
  plumed_assert(gmax[i]>gmin[i]);
  plumed_assert(nbin[i]>0);
 }
 dimension_=gmax.size();
 min_=gmin;
 max_=gmax;
 nbin_=nbin;
 pbc_=pbc;
 dospline_=dospline;
 usederiv_=usederiv;
 if(dospline_) plumed_assert(dospline_==usederiv_);
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

vector<bool> Grid::getIsPeriodic() const {
 return pbc_;
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
unsigned Grid::getIndex(const vector<unsigned> & indices) const {
 plumed_assert(indices.size()==dimension_);
 unsigned index=indices[dimension_-1];
 for(unsigned int i=dimension_-1;i>0;--i){
  index=index*nbin_[i-1]+indices[i-1];
 }
 return index;
}

unsigned Grid::getIndex(const vector<double> & x) const {
 plumed_assert(x.size()==dimension_);
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

vector<unsigned> Grid::getIndices(const vector<double> & x) const {
 plumed_assert(x.size()==dimension_);
 vector<unsigned> indices;
 for(unsigned int i=0;i<dimension_;++i){
   indices.push_back(unsigned(floor((x[i]-min_[i])/dx_[i])));
 }
 return indices;
}

vector<double> Grid::getPoint(const vector<unsigned> & indices) const {
 plumed_assert(indices.size()==dimension_);
 vector<double> x;
 for(unsigned int i=0;i<dimension_;++i){
  x.push_back(min_[i]+(double)(indices[i])*dx_[i]);
 }
 return x;
}

vector<double> Grid::getPoint(unsigned index) const {
 plumed_assert(index<maxsize_);
 return getPoint(getIndices(index));
}

vector<double> Grid::getPoint(const vector<double> & x) const {
 plumed_assert(x.size()==dimension_);
 return getPoint(getIndices(x));
}

void Grid::getPoint(unsigned index,std::vector<double> & point) const{
 plumed_assert(index<maxsize_);
 getPoint(getIndices(index),point);
}

void Grid::getPoint(const std::vector<unsigned> & indices,std::vector<double> & point) const{
 plumed_assert(indices.size()==dimension_);
 plumed_assert(point.size()==dimension_);
 for(unsigned int i=0;i<dimension_;++i){
  point[i]=(min_[i]+(double)(indices[i])*dx_[i]);
 }
}

void Grid::getPoint(const std::vector<double> & x,std::vector<double> & point) const{
 plumed_assert(x.size()==dimension_);
 getPoint(getIndices(x),point);
}


vector<unsigned> Grid::getNeighbors
 (const vector<unsigned> &indices,const vector<unsigned> &nneigh)const{
 plumed_assert(indices.size()==dimension_ && nneigh.size()==dimension_);

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
 (const vector<double> & x,const vector<unsigned> & nneigh)const{
 plumed_assert(x.size()==dimension_ && nneigh.size()==dimension_);
 return getNeighbors(getIndices(x),nneigh);
}

vector<unsigned> Grid::getNeighbors
 (unsigned index,const vector<unsigned> & nneigh)const{
 plumed_assert(index<maxsize_ && nneigh.size()==dimension_);
 return getNeighbors(getIndices(index),nneigh);
}

vector<unsigned> Grid::getSplineNeighbors(const vector<unsigned> & indices)const{
 plumed_assert(indices.size()==dimension_);
 vector<unsigned> neighbors;
 unsigned nneigh=unsigned(pow(2.0,int(dimension_)));
 
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

double Grid::getValue(unsigned index) const {
 plumed_assert(index<maxsize_);
 return grid_[index];
}

double Grid::getValue(const vector<unsigned> & indices) const {
 return getValue(getIndex(indices));
}

double Grid::getValue(const vector<double> & x) const {
 if(!dospline_){
  return getValue(getIndex(x));
 } else {
  vector<double> der(dimension_);
  return getValueAndDerivatives(x,der);
 }
}

double Grid::getValueAndDerivatives
 (unsigned index, vector<double>& der) const{
 plumed_assert(index<maxsize_ && usederiv_ && der.size()==dimension_);
 der=der_[index];
 return grid_[index];
}

double Grid::getValueAndDerivatives
 (const vector<unsigned> & indices, vector<double>& der) const{
 return getValueAndDerivatives(getIndex(indices),der);
}

double Grid::getValueAndDerivatives
(const vector<double> & x, vector<double>& der) const {
 plumed_assert(der.size()==dimension_ && usederiv_);
 
 if(dospline_){
  double X,X2,X3,value;
  vector<double> fd(dimension_);
  vector<double> C(dimension_);
  vector<double> D(dimension_);
  vector<double> dder(dimension_);
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
  return value;
 }else{
  return getValueAndDerivatives(getIndex(x),der);
 }
}

void Grid::setValue(unsigned index, double value){
 plumed_assert(index<maxsize_ && !usederiv_);
 grid_[index]=value;
}

void Grid::setValue(const vector<unsigned> & indices, double value){
 setValue(getIndex(indices),value); 
}

void Grid::setValueAndDerivatives
 (unsigned index, double value, vector<double>& der){
 plumed_assert(index<maxsize_ && usederiv_ && der.size()==dimension_);
 grid_[index]=value;
 der_[index]=der;
}

void Grid::setValueAndDerivatives
 (const vector<unsigned> & indices, double value, vector<double>& der){
 setValueAndDerivatives(getIndex(indices),value,der); 
}

void Grid::addValue(unsigned index, double value){
 plumed_assert(index<maxsize_ && !usederiv_);
 grid_[index]+=value;
}

void Grid::addValue(const vector<unsigned> & indices, double value){
 addValue(getIndex(indices),value);
}

void Grid::addValueAndDerivatives
 (unsigned index, double value, vector<double>& der){
 plumed_assert(index<maxsize_ && usederiv_ && der.size()==dimension_);
 grid_[index]+=value;
 for(unsigned int i=0;i<dimension_;++i) der_[index][i]+=der[i];
}

void Grid::addValueAndDerivatives
 (const vector<unsigned> & indices, double value, vector<double>& der){
 addValueAndDerivatives(getIndex(indices),value,der);
}

void Grid::writeHeader(FILE* file){
 fprintf(file,"#! DERIVATIVE %d\n",int(usederiv_));
 fprintf(file,"#! NVAR      %2u\n",dimension_);
 fprintf(file,"#! BIN"); 
 for(unsigned i=0;i<dimension_;++i){fprintf(file," %14u",nbin_[i]);}
 fprintf(file,"\n");
 fprintf(file,"#! MIN");
 for(unsigned i=0;i<dimension_;++i){fprintf(file," %14.9f",min_[i]);}                        
 fprintf(file,"\n");
 fprintf(file,"#! MAX");
 for(unsigned i=0;i<dimension_;++i){fprintf(file," %14.9f",max_[i]);} 
 fprintf(file,"\n");
 fprintf(file,"#! PBC");
 for(unsigned i=0;i<dimension_;++i){fprintf(file," %14d",int(pbc_[i]));}
 fprintf(file,"\n");
}

void Grid::writeToFile(FILE* file){
 vector<double> xx(dimension_);
 vector<double> der(dimension_);
 double f;
 writeHeader(file);
 for(unsigned i=0;i<getSize();++i){
   xx=getPoint(i);
   if(usederiv_){f=getValueAndDerivatives(i,der);} 
   else{f=getValue(i);}
   if(dimension_>1 && getIndices(i)[dimension_-2]==0){fprintf(file,"\n");} 
   for(unsigned j=0;j<dimension_;++j){fprintf(file,"%14.9f ",xx[j]);}
   fprintf(file,"  %14.9f  ",f);
   if(usederiv_){for(unsigned j=0;j<dimension_;++j){fprintf(file,"%14.9f ",der[j]);}}
   fprintf(file,"\n");
 }
}

Grid* Grid::create(FILE* file, bool dosparse, bool dospline, bool doder)
{
 Grid* grid=NULL;
 unsigned nvar,ibool;
 char str1[50],str2[50];
 fscanf(file,"%s %s %u",str1,str2,&ibool);
 bool hasder=bool(ibool);
 if(doder){plumed_assert(doder==hasder);}
 fscanf(file,"%s %s %u",str1,str2,&nvar);

 vector<unsigned> gbin(nvar);
 vector<double>   gmin(nvar),gmax(nvar);
 vector<bool>     gpbc(nvar);
 fscanf(file,"%s %s",str1,str2);
 for(unsigned i=0;i<nvar;++i){fscanf(file,"%u",&gbin[i]);}
 fscanf(file,"%s %s",str1,str2);
 for(unsigned i=0;i<nvar;++i){fscanf(file,"%lf",&gmin[i]);}
 fscanf(file,"%s %s",str1,str2);
 for(unsigned i=0;i<nvar;++i){fscanf(file,"%lf",&gmax[i]);}
 fscanf(file,"%s %s",str1,str2);
 for(unsigned i=0;i<nvar;++i){fscanf(file,"%u",&ibool);gpbc[i]=bool(ibool);}

 if(!dosparse){grid=new Grid(gmin,gmax,gbin,gpbc,dospline,doder);}
 else{grid=new SparseGrid(gmin,gmax,gbin,gpbc,dospline,doder);}

 vector<double> xx(nvar),dder(nvar);
 vector<double> dx=grid->getDx();
 double f,x;
 while(1){
  int nread;
  for(unsigned i=0;i<nvar;++i){nread=fscanf(file,"%lf",&x);xx[i]=x+dx[i]/2.0;}
  if(nread<1){break;}
  fscanf(file,"%lf",&f);
  if(hasder){for(unsigned i=0;i<nvar;++i){fscanf(file,"%lf",&dder[i]);}}
  unsigned index=grid->getIndex(xx);
  if(doder){grid->setValueAndDerivatives(index,f,dder);}
  else{grid->setValue(index,f);}
 }
 return grid;
}


// Sparse version of grid with map
void SparseGrid::clear(){
 map_.clear();
}

unsigned SparseGrid::getSize() const{
 return map_.size(); 
}

unsigned SparseGrid::getMaxSize() const {
 return maxsize_; 
}

double SparseGrid::getValue(unsigned index)const{
 plumed_assert(index<maxsize_);
 double value=0.0;
 iterator it=map_.find(index);
 if(it!=map_.end()) value=it->second;
 return value;
}

double SparseGrid::getValueAndDerivatives
 (unsigned index, vector<double>& der)const{
 plumed_assert(index<maxsize_ && usederiv_ && der.size()==dimension_);
 double value=0.0;
 for(unsigned int i=0;i<dimension_;++i) der[i]=0.0;
 iterator it=map_.find(index);
 if(it!=map_.end()) value=it->second;
 iterator_der itder=der_.find(index);
 if(itder!=der_.end()) der=itder->second;
 return value;
}

void SparseGrid::setValue(unsigned index, double value){
 plumed_assert(index<maxsize_ && !usederiv_);
 map_[index]=value;
}

void SparseGrid::setValueAndDerivatives
 (unsigned index, double value, vector<double>& der){
 plumed_assert(index<maxsize_ && usederiv_ && der.size()==dimension_);
 map_[index]=value;
 der_[index]=der;
}

void SparseGrid::addValue(unsigned index, double value){
 plumed_assert(index<maxsize_ && !usederiv_);
 map_[index]+=value;
}

void SparseGrid::addValueAndDerivatives
 (unsigned index, double value, vector<double>& der){
 plumed_assert(index<maxsize_ && usederiv_ && der.size()==dimension_);
 map_[index]+=value;
 der_[index].resize(dimension_);
 for(unsigned int i=0;i<dimension_;++i) der_[index][i]+=der[i]; 
}

void SparseGrid::writeToFile(FILE* file){
 vector<double> xx(dimension_);
 vector<double> der(dimension_);
 double f;
 writeHeader(file);
 for(iterator it=map_.begin();it!=map_.end();it++){
   unsigned i=(*it).first;
   xx=getPoint(i);
   if(usederiv_){f=getValueAndDerivatives(i,der);} 
   else{f=getValue(i);}
   if(dimension_>1 && getIndices(i)[dimension_-2]==0){fprintf(file,"\n");}
   for(unsigned j=0;j<dimension_;++j){fprintf(file,"%14.9f ",xx[j]);}
   fprintf(file,"  %14.9f  ",f);
   if(usederiv_){for(unsigned j=0;j<dimension_;++j){fprintf(file,"%14.9f ",der[j]);}}
   fprintf(file,"\n");
 }
}
