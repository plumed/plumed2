#include <vector>
#include <algorithm>
#include "Vector.h"
#include "Pbc.h"
#include "NeighborList.h"

using namespace PLMD;
using namespace std;

NeighborList::NeighborList(vector<int> list0, vector<int> list1,
                           double distance, unsigned stride, Pbc *pbc=NULL):  
                           list0_(list0), list1_(list1),
                           distance_(distance), pbc_(pbc), stride_(stride)
{
 if(stride_<2){
  // i need to access to log here. still don't know how to do it
  //log.printf("stride should be greater or equal to 2")
 }
// find maximum index atom
 int imax0= *max_element(list0_.begin(),list0_.end());
 int imax1= *max_element(list1_.begin(),list1_.end());
 indexes_.resize(max(imax0,imax1)+1);
 neighbors_.resize(list0_.size());
// initialize neighbor list with all the atoms
 for(int i=0;i<list0_.size();++i){
  for(int j=0;j<list1_.size();++j){
    neighbors_[i].push_back(list1_[j]);
  }
 }
}

vector<int> NeighborList::getFullList()
{
 indexes_.clear();
 vector<int> request=list0_;
 request.insert(request.end(),list1_.begin(),list1_.end());
 for(int i=0;i<request.size();++i) indexes_[request[i]]=i; 
 return request;
}

vector<int> NeighborList::update(vector<Vector> positions)
{
 vector<int> request=list0_;
 neighbors_.clear();
 neighbors_.resize(list0_.size());
// check if positions array has the correct length
 if(positions.size()!=list0_.size()+list1_.size()){
  // i need to access to log here. still don't know how to do it
  //log.printf("ERROR you need to request all the atoms before updating NL\n"); 
 }
// update neighbors
 for(int i=0;i<list0_.size();++i){
  for(int j=list0_.size();j<list0_.size()+list1_.size();++j){
   Vector distance;
   if(pbc_!=NULL){
    distance=pbc_->distance(positions[i],positions[j]);
   } else {
    distance=delta(positions[i],positions[j]);
   }
   double value=distance.modulo();
   if(value>0.0 && value<=distance_) {
    neighbors_[i].push_back(list1_[j-list0_.size()]);
   } 
  }
 request.insert(request.end(),neighbors_[i].begin(),neighbors_[i].end());
 }
 indexes_.clear();
 for(int i=0;i<request.size();++i) indexes_[request[i]]=i;
 return request;
}

int NeighborList::getStride() const
{
 return stride_;
}

vector<int> NeighborList::getNeighbors(int index)
{
 vector<int> neigh;
 for(int i=0;i<getNumberOfNeighbors(index);++i){
  int iatom=neighbors_[index][i];
  neigh.push_back(indexes_[iatom]);
 }
 return neigh; 
}

int NeighborList::getNumberOfNeighbors(int index)
{
 return neighbors_[index].size();
}

int NeighborList::getNumberOfAtoms() const
{
 return neighbors_.size();
}
