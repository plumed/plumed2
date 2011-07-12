#include <vector>
#include <algorithm>
#include "Vector.h"
#include "Pbc.h"
#include "NeighborList.h"

using namespace PLMD;
using namespace std;

NeighborList::NeighborList(vector<unsigned> list0, vector<unsigned> list1,
                           double distance, unsigned stride, Pbc *pbc=NULL):  
                           list0_(list0), list1_(list1),
                           distance_(distance), pbc_(pbc), stride_(stride)
{
// find maximum index atom
 unsigned imax0= *max_element(list0_.begin(),list0_.end());
 unsigned imax1= *max_element(list1_.begin(),list1_.end());
// initialize indexes
 imax_=max(imax0,imax1)+1;
 indexes_.resize(imax_);
 newindexes_.resize(imax_);
// initialize neighbor list with all the atoms
// NB: the neighbor list contains the atomic index of the atoms
// not the index in the array of positions
 neighbors_.resize(list0_.size());
 for(unsigned int i=0;i<list0_.size();++i){
  for(unsigned int j=0;j<list1_.size();++j){
    neighbors_[i].push_back(list1_[j]);
  }
 }
}

// indexes[i] is the index in the positions array (the one requested from
// the MD code) of the atom with "atomic" index i 
vector<unsigned> NeighborList::createIndexes(vector<unsigned> request)
{
 vector<unsigned> indexes;
 indexes.resize(imax_);
 for(unsigned int i=0;i<request.size();++i){
  indexes[request[i]]=i;
 }
 return indexes;
}

// this method should be called at the end of the step
// before the update of the neighbor list
vector<unsigned> NeighborList::getFullList()
{
 vector<unsigned> request=list0_;
 request.insert(request.end(),list1_.begin(),list1_.end());
 indexes_=createIndexes(request);
 return request;
}

// this method should be called at the beginning of the step
// when the neighbor list should be updated
vector<unsigned> NeighborList::getUpdatedList(vector<Vector> positions)
{
 vector<unsigned> request=list0_;
// clean neighbors list
 neighbors_.clear();
 neighbors_.resize(list0_.size());
// check if positions array has the correct length
 if(positions.size()!=list0_.size()+list1_.size()){
  // i need to access to log here. still don't know how to do it
  //log.printf("ERROR you need to request all the atoms one step before updating the NL\n"); 
 }
// update neighbors
 for(unsigned int i=0;i<list0_.size();++i){
  for(unsigned int j=list0_.size();j<list0_.size()+list1_.size();++j){
   Vector distance;
   if(pbc_!=NULL){
    distance=pbc_->distance(positions[i],positions[j]);
   } else {
    distance=delta(positions[i],positions[j]);
   }
   double value=distance.modulo();
   if(value<=distance_) {
    neighbors_[i].push_back(list1_[j-list0_.size()]);
   } 
  }
  request.insert(request.end(),neighbors_[i].begin(),neighbors_[i].end());
 }
 newindexes_=createIndexes(request);
 return request;
}

// this method should be called at the end of the step
// when nl must be updated, before getFullList()
void NeighborList::update()
{
 indexes_=newindexes_;
}

unsigned NeighborList::getStride() const
{
 return stride_;
}

vector<unsigned> NeighborList::getNeighbors(unsigned index)
{
 vector<unsigned> neigh;
 for(unsigned int i=0;i<getNumberOfNeighbors(index);++i){
  unsigned iatom=neighbors_[index][i];
  neigh.push_back(indexes_[iatom]);
 }
 return neigh; 
}

unsigned NeighborList::getNumberOfNeighbors(unsigned index)
{
 return neighbors_[index].size();
}

unsigned NeighborList::getNumberOfAtoms() const
{
 return neighbors_.size();
}
