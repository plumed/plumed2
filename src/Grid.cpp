#include <vector>
#include <cassert>
#include <math.h>
#include "Grid.h"
#include <iostream>

using namespace std;


Grid::Grid(vector<double> gmin,vector<double> gmax,vector<unsigned> nbin,
           bool dospline, bool doclear):
           dospline_(dospline) {
// various checks
 assert(gmax.size()==gmin.size());
 assert(gmax.size()==nbin.size());
 for(unsigned int i=0;i<gmax.size();++i){
  assert(gmax[i]>gmin[i]);
  assert(nbin[i]>0);
 }
 min_=gmin;
 max_=gmax;
 nbin_=nbin;
 maxsize_=1;
 for(unsigned int i=0;i<max_.size();++i){
  dx_.push_back((max_[i]-min_[i])/(double)(nbin_[i]-1));
  maxsize_*=nbin_[i];
 }
 if(doclear) {clear();}
}

Grid::Grid(vector<double> gmin,vector<double> gmax,vector<double> dx,
           bool dospline, bool doclear):
           dospline_(dospline) {
// various checks
 assert(gmax.size()==gmin.size());
 assert(gmax.size()==dx.size());
 for(unsigned int i=0;i<gmax.size();++i){
  assert(gmax[i]>gmin[i]);
  assert(dx[i]>0.0);
 }
 min_=gmin;
 max_=gmax;
 dx_=dx;
 maxsize_=1;
 for(unsigned int i=0;i<max_.size();++i){
  nbin_.push_back(floor((max_[i]-min_[i])/dx_[i])+1);
  maxsize_*=nbin_[i];
 }
 if(doclear) {clear();}
}

void Grid::clear(){
 grid_.resize(maxsize_);
 for(unsigned int i=0;i<size();++i){grid_[i]=0.0;}
}

vector<double> Grid::getMin() const {
 return min_;
}

vector<double> Grid::getMax() const {
 return max_;
}

vector<double> Grid::getSide() const {
 return dx_;
}

unsigned Grid::size() const {
 return grid_.size();
}

unsigned Grid::dimension() const {
 return max_.size();
}

// we are flattening arrays using a column-major order
unsigned Grid::getIndex(vector<unsigned> indices) const {
 assert(indices.size()==dimension());
 unsigned index=indices[dimension()-1];
 for(unsigned int i=dimension()-1;i>0;--i){
  index=index*nbin_[i-1]+indices[i-1];
 }
 return index;
}

unsigned Grid::getIndex(vector<double> x) const {
 assert(x.size()==dimension());
 return getIndex(getIndices(x));
}

// we are flattening arrays using a column-major order
vector<unsigned> Grid::getIndices(unsigned index) const {
 vector<unsigned> indices;
 unsigned N=dimension();
 unsigned kk=index;
 indices.push_back(index%nbin_[0]);
 for(unsigned int i=1;i<N-1;++i){
  kk=(kk-indices[i-1])/nbin_[i-1];
  indices.push_back(kk%nbin_[i]);
 }
 if(N>=2) indices.push_back((kk-indices[N-2])/nbin_[N-2]);
 return indices;
}

vector<unsigned> Grid::getIndices(vector<double> x) const {
 assert(x.size()==dimension());
 vector<unsigned> indices;
 for(unsigned int i=0;i<dimension();++i){
   indices.push_back(floor((x[i]-min_[i])/dx_[i]));
 }
 return indices;
}

vector<double> Grid::getPoint(vector<unsigned> indices) const {
 assert(indices.size()==dimension());
 vector<double> x;
 for(unsigned int i=0;i<dimension();++i){
  x.push_back(min_[i]+(double)(indices[i])*dx_[i]);
 }
 return x;
}

vector<double> Grid::getPoint(unsigned index) const {
 assert(index>=0 && index<maxsize_);
 return getPoint(getIndices(index));
}

vector<double> Grid::getPoint(vector<double> x) const {
 assert(x.size()==dimension());
 return getPoint(getIndices(x));
}

double Grid::getValue(unsigned index) {
 assert(index>=0 && index<size());
 return grid_[index];
}

double Grid::getValue(vector<unsigned> indices) {
 return getValue(getIndex(indices));
}

double Grid::getValue(vector<double> x) {
 if(dospline_){
// TO DO 
 } else {return getValue(getIndices(x));}
}

void Grid::setValue(unsigned index, double value){
 assert(index>=0 && index<size());
 grid_[index]=value;
}

void Grid::setValue(vector<unsigned> indices, double value){
 setValue(getIndex(indices),value); 
}

void Grid::addValue(unsigned index, double value){
 assert(index>=0 && index<size());
 grid_[index]+=value;
}

void Grid::addValue(vector<unsigned> indices, double value){
 addValue(getIndex(indices),value);
}

// Sparse version of grid with map
void SparseGrid::clear(){
 map_.clear();
}

unsigned SparseGrid::size() const{
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
 if(dospline_){
// TO DO 
 } else {return getValue(getIndices(x));}
}

void SparseGrid::setValue(unsigned index, double value){
 assert(index>=0 && index<maxsize_);
 map_[index]=value;
}

void SparseGrid::setValue(vector<unsigned> indices, double value){
 setValue(getIndex(indices),value); 
}

void SparseGrid::addValue(unsigned index, double value){
 assert(index>=0 && index<maxsize_);
 map_[index]+=value;
}

void SparseGrid::addValue(vector<unsigned> indices, double value){
 addValue(getIndex(indices),value);
}

double SparseGrid::getMaxSize() const {
 return maxsize_;
}
