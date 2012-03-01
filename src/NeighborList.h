#ifndef __PLUMED_NeighborList_h
#define __PLUMED_NeighborList_h

#include <vector>
#include <algorithm>
#include "Tools.h"

namespace PLMD{

//+DEVELDOC TOOLBOX NeighborList
/**
A class for creating neighbour lists of atoms.
*/
//+ENDDEVELDOC


/// A class that implements neighbor lists from two lists or a single list of atoms
template <typename T, typename D, typename I>
class NeighborList  
{
  bool do_pair_,do_pbc_,twolists_;
  const D& pbc_;
  std::vector<I> fullatomlist_,requestlist_;
  std::vector<std::pair<unsigned,unsigned> > neighbors_;
  double distance_;
  unsigned stride_,nlist0_,nlist1_,nallpairs_,lastupdate_;
/// Initialize the neighbor list with all possible pairs
  void initialize();
/// Return the pair of indexes in the positions array
/// of the two atoms forming the i-th pair among all possible pairs  
  std::pair<unsigned,unsigned> getIndexPair(unsigned i);
/// Extract the list of atoms from the current list of close pairs
  void setRequestList();
public:
  NeighborList(const std::vector<I>& list0,
               const std::vector<I>& list1,
               const bool& do_pair, const bool& do_pbc, const D& pbc, 
               const double& distance=1.0e+30, const unsigned& stride=0);
  NeighborList(const std::vector<I>& list0, const bool& do_pbc,
               const D& pbc, const double& distance=1.0e+30,
               const unsigned& stride=0);
/// Return the list of all atoms. These are needed to rebuild the neighbor list.                         
  std::vector<I>& getFullAtomList();
/// Update the indexes in the neighbor list to match the
/// ordering in the new positions array
/// and return the new list of atoms that must be requested to the main code
  std::vector<I>& getReducedAtomList();
/// Update the neighbor list and prepare the new
/// list of atoms that will be requested to the main code  
  void update(const std::vector<T>& positions);
/// Get the update stride of the neighbor list
  unsigned getStride() const;
/// Get the last step in which the neighbor list was updated  
  unsigned getLastUpdate() const;
/// Set the step of the last update
  void setLastUpdate(unsigned step);
/// Get the size of the neighbor list
  unsigned size() const;
/// Get the i-th pair of the neighbor list
  std::pair<unsigned,unsigned> getClosePair(unsigned i) const;
/// Get the list of neighbors of the i-th atom
  std::vector<unsigned> getNeighbors(unsigned i);
  ~NeighborList(){};
};


template <typename T, typename D, typename I>
NeighborList<T,D,I>::NeighborList(const std::vector<I>& list0, const std::vector<I>& list1,
                           const bool& do_pair, const bool& do_pbc, const D& pbc,
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
  plumed_assert(nlist0_==nlist1_);
  nallpairs_=nlist0_;
 }
 initialize();
 lastupdate_=0;
}

template <typename T, typename D, typename I>
NeighborList<T,D,I>::NeighborList(const std::vector<I>& list0, const bool& do_pbc,
                           const D& pbc, const double& distance,
                           const unsigned& stride):
                           do_pbc_(do_pbc), pbc_(pbc),
                           distance_(distance), stride_(stride){
 fullatomlist_=list0;
 nlist0_=list0.size();
 twolists_=false;
 nallpairs_=nlist0_*(nlist0_-1)/2;
 initialize();
 lastupdate_=0;
}

template <typename T, typename D, typename I>
void NeighborList<T,D,I>::initialize() {
 neighbors_.clear();
 for(unsigned int i=0;i<nallpairs_;++i){
   neighbors_.push_back(getIndexPair(i));
 }
}

template <typename T, typename D, typename I>
std::vector<I>& NeighborList<T,D,I>::getFullAtomList() {
 return fullatomlist_;
}

template <typename T, typename D, typename I>
std::pair<unsigned,unsigned> NeighborList<T,D,I>::getIndexPair(unsigned ipair) {
 std::pair<unsigned,unsigned> index;
 if(twolists_ && do_pair_){
  index=std::pair<unsigned,unsigned>(ipair,ipair+nlist0_);
 }else if (twolists_ && !do_pair_){
  index=std::pair<unsigned,unsigned>(ipair/nlist1_,ipair%nlist1_+nlist0_);
 }else if (!twolists_){
  unsigned ii = nallpairs_-1-ipair;
  unsigned  K = unsigned(floor((sqrt(double(8*ii+1))+1)/2));
  unsigned jj = ii-K*(K-1)/2;
  index=std::pair<unsigned,unsigned>(nlist0_-1-K,nlist0_-1-jj);
 }
 return index;
}

template <typename T, typename D, typename I>
void NeighborList<T,D,I>::update(const std::vector<T>& positions) {
 neighbors_.clear();
// check if positions array has the correct length 
 plumed_assert(positions.size()==fullatomlist_.size());
 for(unsigned int i=0;i<nallpairs_;++i){
   std::pair<unsigned,unsigned> index=getIndexPair(i);
   unsigned index0=index.first;
   unsigned index1=index.second;
   T distance;
   if(do_pbc_){
    distance=pbc_.distance(positions[index0],positions[index1]);
   } else {
    distance=delta(positions[index0],positions[index1]);
   }
   double value=distance.modulo();
   if(value<=distance_) {neighbors_.push_back(index);} 
 }
 setRequestList();
}

template <typename T, typename D, typename I>
void NeighborList<T,D,I>::setRequestList() {
 requestlist_.clear();
 for(unsigned int i=0;i<size();++i){
  requestlist_.push_back(fullatomlist_[neighbors_[i].first]);
  requestlist_.push_back(fullatomlist_[neighbors_[i].second]);
 }
 Tools::removeDuplicates(requestlist_);
}

template <typename T, typename D, typename I>
std::vector<I>& NeighborList<T,D,I>::getReducedAtomList() {
 std::vector< std::pair<unsigned,unsigned> > newneighbors;
 for(unsigned int i=0;i<size();++i){
  unsigned newindex0,newindex1;
  I index0=fullatomlist_[neighbors_[i].first];
  I index1=fullatomlist_[neighbors_[i].second];
  for(unsigned j=0;j<requestlist_.size();++j){
   if(requestlist_[j]==index0) newindex0=j;
   if(requestlist_[j]==index1) newindex1=j;
  }
  newneighbors.push_back(std::pair<unsigned,unsigned>(newindex0,newindex1));
 }
 neighbors_.clear();
 neighbors_=newneighbors;
 return requestlist_;
}

template <typename T, typename D, typename I>
unsigned NeighborList<T,D,I>::getStride() const {
 return stride_;
}

template <typename T, typename D, typename I>
unsigned NeighborList<T,D,I>::getLastUpdate() const {
 return lastupdate_;
}

template <typename T, typename D, typename I>
void NeighborList<T,D,I>::setLastUpdate(unsigned step) {
 lastupdate_=step;
}

template <typename T, typename D, typename I>
unsigned NeighborList<T,D,I>::size() const {
 return neighbors_.size();
}

template <typename T, typename D, typename I>
std::pair<unsigned,unsigned> NeighborList<T,D,I>::getClosePair(unsigned i) const {
 return neighbors_[i];
}

template <typename T, typename D, typename I>
std::vector<unsigned> NeighborList<T,D,I>::getNeighbors(unsigned index) {
 std::vector<unsigned> neighbors;
 for(unsigned int i=0;i<size();++i){
  if(neighbors_[i].first==index)  neighbors.push_back(neighbors_[i].second);
  if(neighbors_[i].second==index) neighbors.push_back(neighbors_[i].first);
 }
 return neighbors;
}


}

#endif
