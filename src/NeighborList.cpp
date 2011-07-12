#include <vector>
#include <algorithm>
#include "Vector.h"
#include "Pbc.h"
#include "NeighborList.h"

using namespace PLMD;
using namespace std;

NeighborList::NeighborList(vector<unsigned> list0, vector<unsigned> list1, bool do_pair,
                           Pbc *pbc=NULL, double distance=1.0e+20, unsigned stride=0):  
                           do_pair_(do_pair), pbc_(pbc), distance_(distance), stride_(stride)
{
// store full list of atoms needed
 fullatomlist_=list0;
 fullatomlist_.insert(fullatomlist_.end(),list1.begin(),list1.end());
 nlist0_=list0.size();
 nlist1_=list1.size();
 twolists_=true;
 if(!do_pair){
  nallpairs_=nlist0_*nlist1_;
 }else{
  // if(nlist0_!=nlist1_) error
  nallpairs_=nlist0_;
 }
}

NeighborList::NeighborList(vector<unsigned> list0, Pbc *pbc=NULL,
                           double distance=1.0e+20, unsigned stride=0):  
                           pbc_(pbc), distance_(distance), stride_(stride)
{
// store full list of atoms needed
 fullatomlist_=list0;
 nlist0_=list0.size();
 twolists_=false;
 nallpairs_=nlist0_*(nlist0_-1)/2;
}

// this method should be called at the end of the step
// before the update of the neighbor list or in the prep stage
vector<unsigned> NeighborList::getFullAtomList() const
{
 return fullatomlist_;
}

// this method return the pair of indexes in the array of positions
// of a given pair "ipair" from the complete list
pair<unsigned,unsigned> NeighborList::getIndexPair(unsigned ipair)
{
 pair<unsigned,unsigned> index;
 if(twolists_ && do_pair_){
  index=pair<unsigned,unsigned>(ipair,ipair+nlist0_);
 }else if (twolists_ && !do_pair_){
  index=pair<unsigned,unsigned>(ipair/nlist1_,ipair%nlist1_+nlist0_);
 }else if (!twolists_){
// still to fix the map between mono and multi index
  index=pair<unsigned,unsigned>(ipair/nlist1_,ipair%nlist1_+nlist0_);
 }
 return index;
}

// this method should be called at the beginning of the step
// when the neighbor list is updated
vector<unsigned> NeighborList::getReducedAtomList(vector<Vector> positions)
{
// clean neighbors list
 neighbors_.clear();
// check if positions array has the correct length
 if(positions.size()!=fullatomlist_.size()){
  // i need to access to log here. still don't know how to do it
  //log.printf("ERROR you need to request all the atoms one step before updating the NL\n"); 
 }
// update neighbors
 for(unsigned int i=0;i<nallpairs_;++i){
   pair<unsigned,unsigned> index=getIndexPair(i);
   unsigned index0=index.first;
   unsigned index1=index.second;
   Vector distance;
   if(pbc_!=NULL){
    distance=pbc_->distance(positions[index0],positions[index1]);
   } else {
    distance=delta(positions[index0],positions[index1]);
   }
   double value=distance.modulo();
   if(value<=distance_) {
    neighbors_.push_back(index);
   } 
 }
// I need the new list of atoms to request
 requestlist_=getRequestList();
 return requestlist_;
}

template<typename T>
void NeighborList::removeDuplicates(std::vector<T>& vec)
{
    std::sort(vec.begin(), vec.end());
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}

vector<unsigned> NeighborList::getRequestList()
{
 vector<unsigned> request;
 for(unsigned int i=0;i<size();++i){
  request.push_back(fullatomlist_[neighbors_[i].first]);
  request.push_back(fullatomlist_[neighbors_[i].second]);
 }
 removeDuplicates<unsigned>(request);
 return request;
}

// this method should be called at the end of the step when nl is updated
void NeighborList::update()
{
 std::vector< pair<unsigned,unsigned> > newneighbors;
 for(unsigned int i=0;i<size();++i){
  unsigned newindex0,newindex1;
  unsigned index0=fullatomlist_[neighbors_[i].first];
  unsigned index1=fullatomlist_[neighbors_[i].second];
  for(unsigned j=0;j<requestlist_.size();++j){
   if(requestlist_[j]==index0) newindex0=j;
   if(requestlist_[j]==index1) newindex1=j;
  }
  newneighbors.push_back(pair<unsigned,unsigned>(newindex0,newindex1));
 }
 neighbors_.clear();
 neighbors_=newneighbors;
}

unsigned NeighborList::getStride() const
{
 return stride_;
}

unsigned NeighborList::size() const
{
 return neighbors_.size();
}

const pair<unsigned,unsigned> & NeighborList::operator[] (unsigned i) const
{
 return neighbors_[i];
}

vector<unsigned> NeighborList::getNeighbors(unsigned index)
{
 vector<unsigned> neighbors;
 for(unsigned int i=0;i<size();++i){
  if(neighbors_[i].first==index)  neighbors.push_back(neighbors_[i].second);
  if(neighbors_[i].second==index) neighbors.push_back(neighbors_[i].first);
 }
 return neighbors;
}