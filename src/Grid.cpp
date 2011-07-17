#include <vector>
#include <cassert>
#include <math.h>
#include "Grid.h"

using namespace std;


Grid::Grid(vector<double> gmin,vector<double> gmax,vector<unsigned> nbin,
           bool dospline, bool doclear):
           dospline_(dospline) {
// various checks
 assert(gmax.size()==gmin.size());
 assert(gmax.size()==nbin.size());
 for(unsigned int i=0;i<gmax.size();++i){
  assert(gmax[i]>gmin[i]);
  assert(nbin[i]>0.0);
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

void Grid::clear(){
 grid_.resize(maxsize_);
 for(unsigned int i=0;i<size();++i){grid_[i]=0.0;}
}

double Grid::getMin(unsigned i) const {
 return min_[i];
}

double Grid::getMax(unsigned i) const {
 return max_[i];
}

double Grid::getSide(unsigned i) const {
 return dx_[i];
}

unsigned Grid::size() const {
 return grid_.size();
}

unsigned Grid::dimension() const {
 return max_.size();
}

unsigned Grid::getIndex(vector<unsigned> indices) const {
 assert(indices.size()==dimension());
 unsigned index;
 // TO DO
 return index;
}

unsigned Grid::getIndex(vector<double> x) const {
 assert(x.size()==dimension());
 vector<unsigned> indices=getIndices(x);
 return getIndex(indices);
}

vector<unsigned> Grid::getIndices(unsigned index) const {
 vector<unsigned> indices;
 // TO DO
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
 vector<unsigned> indices=getIndices(index);
 return getPoint(indices);
}

double Grid::getValue(unsigned index) const { 
 return grid_[index];
}

double Grid::getValue(vector<unsigned> indices) const {
 return grid_[getIndex(indices)];
}

double Grid::getValue(vector<double> x) const{
 assert(x.size()==dimension());
 vector<unsigned> indices;
 double value;
 if(dospline_){
// TO DO 
 } else {
  value=getValue(getIndices(x));
 }
 return value;
}

void Grid::setValue(unsigned index, double value){
 grid_[index]=value;
}

void Grid::setValue(vector<unsigned> indices, double value){
 grid_[getIndex(indices)]=value; 
}

void Grid::addValue(unsigned index, double value){
 grid_[index]+=value;
}

void Grid::addValue(vector<unsigned> indices, double value){
 grid_[getIndex(indices)]+=value;
}

// Sparse version of grid with map
void SparseGrid::clear(){
 map_.clear();
}

unsigned SparseGrid::size() const{
 return map_.size();
}