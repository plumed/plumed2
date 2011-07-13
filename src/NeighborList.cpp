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
 lastupdate_=-1;
}

NeighborList::NeighborList(const vector<AtomNumber>& list0, const bool& do_pbc,
                           const Pbc& pbc, const double& distance,
                           const unsigned& stride):
                           do_pbc_(do_pbc), pbc_(pbc),
                           distance_(distance), stride_(stride)
{
 fullatomlist_=list0;
 nlist0_=list0.size();
 twolists_=false;
 nallpairs_=nlist0_*(nlist0_-1)/2;
 initialize();
 lastupdate_=-1;
}

void NeighborList::initialize()
{
 neighbors_.clear();
 for(unsigned int i=0;i<nallpairs_;++i){
   neighbors_.push_back(getIndexPair(i));
 }
}

vector<AtomNumber>& NeighborList::getFullAtomList()
{
 return fullatomlist_;
}

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

void NeighborList::update(const vector<Vector>& positions, int step)
{
 lastupdate_=step;
 neighbors_.clear();
// check if positions array has the correct length
 if(positions.size()!=fullatomlist_.size()){
  // i need to access to log here. still don't know how to do it
  //log.printf("ERROR you need to request all the atoms one step before updating the NL\n"); 
 }
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

void NeighborList::setRequestList()
{
 requestlist_.clear();
 for(unsigned int i=0;i<size();++i){
  requestlist_.push_back(fullatomlist_[neighbors_[i].first]);
  requestlist_.push_back(fullatomlist_[neighbors_[i].second]);
 }
}

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

unsigned NeighborList::getStride() const
{
 return stride_;
}

bool NeighborList::doUpdate(int step)
{
 bool doit;
 if(stride_>0 && step!=lastupdate_ && (step-lastupdate_)%stride_==0){
  doit=true;
 }else{
  doit=false;
 }
 return doit;
}

unsigned NeighborList::size() const
{
 return neighbors_.size();
}

 pair<unsigned,unsigned> NeighborList::getClosePair(unsigned i) const
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

}