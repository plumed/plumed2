#include <vector>
#include <algorithm>
#include "Vector.h"
#include "Pbc.h"
#include "AtomNumber.h"
#include "NeighborList.h"

using namespace PLMD;
using namespace std;

namespace PLMD{

NeighborList::NeighborList(const vector<AtomNumber>& list0, const vector<AtomNumber>& list1,
                           const bool& do_pair, const bool& do_pbc, const Pbc& pbc,
                           const double& distance, const unsigned& stride):
                           do_pair_(do_pair), do_pbc_(do_pbc), pbc_(pbc),
                           distance_(distance), stride_(stride)
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
 initialize();
}


NeighborList::NeighborList(const vector<AtomNumber>& list0, const bool& do_pbc,
                           const Pbc& pbc, const double& distance,
                           const unsigned& stride):
                           do_pbc_(do_pbc), pbc_(pbc),
                           distance_(distance), stride_(stride)
{
// store full list of atoms needed
 fullatomlist_=list0;
 nlist0_=list0.size();
 twolists_=false;
 nallpairs_=nlist0_*(nlist0_-1)/2;
 initialize();
}

// initialize neighborlist with all possible pairs
void NeighborList::initialize()
{
 neighbors_.clear();
 for(unsigned int i=0;i<nallpairs_;++i){
   neighbors_.push_back(getIndexPair(i));
 }
}

// this methods returns the list of all atoms
// it should be called at the end of the step
// before the update of the neighbor list or in the prep stage
vector<AtomNumber>& NeighborList::getFullAtomList()
{
 return fullatomlist_;
}

// this method returns the pair of indexes of the positions array
// of pair "ipair" among the list of all possible pairs
pair<unsigned,unsigned> NeighborList::getIndexPair(unsigned ipair)
{
 pair<unsigned,unsigned> index;
 if(twolists_ && do_pair_){
  index=pair<unsigned,unsigned>(ipair,ipair+nlist0_);
 }else if (twolists_ && !do_pair_){
  index=pair<unsigned,unsigned>(ipair/nlist1_,ipair%nlist1_+nlist0_);
 }else if (!twolists_){
// still to fix the map between mono and multi index here
  index=pair<unsigned,unsigned>(ipair/nlist1_,ipair%nlist1_+nlist0_);
 }
 return index;
}

// this method updates the neighbor list and prepare the new
// list of atoms to request to the MD code
void NeighborList::update(const vector<Vector>& positions)
{
// clean neighbors list
 neighbors_.clear();
// check if positions array has the correct length
 if(positions.size()!=fullatomlist_.size()){
  // i need to access to log here. still don't know how to do it
  //log.printf("ERROR you need to request all the atoms one step before updating the NL\n"); 
 }
// update neighbors. cycle on all possible pairs
 for(unsigned int i=0;i<nallpairs_;++i){
   pair<unsigned,unsigned> index=getIndexPair(i);
   unsigned index0=index.first;
   unsigned index1=index.second;
   Vector distance;
   if(do_pbc_){
    distance=pbc_.distance(positions[index0],positions[index1]);
   } else {
    distance=delta(positions[index0],positions[index1]);
   }
   double value=distance.modulo();
   if(value<=distance_) {
    neighbors_.push_back(index);
   } 
 }
 setRequestList();
}

// template to remove duplicates from a list of types <T>
template<typename T>
void NeighborList::removeDuplicates(std::vector<T>& vec)
{
    std::sort(vec.begin(), vec.end());
    vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
}

// extract the list of atoms from the list of close pairs
void NeighborList::setRequestList()
{
 requestlist_.clear();
 for(unsigned int i=0;i<size();++i){
  requestlist_.push_back(fullatomlist_[neighbors_[i].first]);
  requestlist_.push_back(fullatomlist_[neighbors_[i].second]);
 }
// is it going to work when T=AtomNumber?
 //removeDuplicates<AtomNumber>(request);
}

// this method should be called at the end of the step when nl is updated
// it updates the indexes in the neighbor list to match the new
// ordering in the future positions array at the next iteration
// and returns the new list of atoms to request to the MD code
vector<AtomNumber>& NeighborList::getReducedAtomList()
{
 std::vector< pair<unsigned,unsigned> > newneighbors;
 for(unsigned int i=0;i<size();++i){
  unsigned newindex0,newindex1;
  unsigned index0=fullatomlist_[neighbors_[i].first].index();
  unsigned index1=fullatomlist_[neighbors_[i].second].index();
  for(unsigned j=0;j<requestlist_.size();++j){
   if(requestlist_[j].index()==index0) newindex0=j;
   if(requestlist_[j].index()==index1) newindex1=j;
  }
  newneighbors.push_back(pair<unsigned,unsigned>(newindex0,newindex1));
 }
 neighbors_.clear();
 neighbors_=newneighbors;
 return requestlist_;
}

// get the update stride of the Neighbor List
unsigned NeighborList::getStride() const
{
 return stride_;
}

// size of the Neighbor List
unsigned NeighborList::size() const
{
 return neighbors_.size();
}

// this method returns the pair of indexes in the positions array of
// the close pair "i"
const pair<unsigned,unsigned> & NeighborList::operator[] (unsigned i) const
{
 return neighbors_[i];
}

// this method returns the list of indexes in the positions array of the
// neighbors of the atom "index", where "index" refers to the index in the
// positions array
vector<unsigned> NeighborList::getNeighbors(unsigned index)
{
 vector<unsigned> neighbors;
 for(unsigned int i=0;i<size();++i){
  if(neighbors_[i].first==index)  neighbors.push_back(neighbors_[i].second);
  if(neighbors_[i].second==index) neighbors.push_back(neighbors_[i].first);
 }
 return neighbors;
}

}