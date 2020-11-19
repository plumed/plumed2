/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "NeighborList.h"
#include "Vector.h"
#include "Pbc.h"
#include "AtomNumber.h"
#include "Communicator.h"
#include "OpenMP.h"
#include "Tools.h"
#include <vector>
#include <algorithm>
#include <numeric>

namespace PLMD {

NeighborList::NeighborList(const std::vector<AtomNumber>& list0, const std::vector<AtomNumber>& list1,
                           const bool& serial, const bool& do_pair, const bool& do_pbc, const Pbc& pbc, Communicator& cm,
                           const double& distance, const unsigned& stride): reduced(false),
  serial_(serial), do_pair_(do_pair), do_pbc_(do_pbc), pbc_(&pbc), comm(cm),
  distance_(distance), stride_(stride)
{
// store full list of atoms needed
  fullatomlist_=list0;
  fullatomlist_.insert(fullatomlist_.end(),list1.begin(),list1.end());
  nlist0_=list0.size();
  nlist1_=list1.size();
  twolists_=true;
  if(!do_pair) {
    nallpairs_=nlist0_*nlist1_;
  } else {
    plumed_assert(nlist0_==nlist1_) << "when using PAIR option, the two groups should have the same number of elements\n"
                                    << "the groups you specified have size "<<nlist0_<<" and "<<nlist1_;
    nallpairs_=nlist0_;
  }
  initialize();
  lastupdate_=0;
}

NeighborList::NeighborList(const std::vector<AtomNumber>& list0, const bool& serial, const bool& do_pbc,
                           const Pbc& pbc, Communicator& cm, const double& distance,
                           const unsigned& stride): reduced(false),
  serial_(serial), do_pbc_(do_pbc), pbc_(&pbc), comm(cm),
  distance_(distance), stride_(stride) {
  fullatomlist_=list0;
  nlist0_=list0.size();
  twolists_=false;
  nallpairs_=nlist0_*(nlist0_-1)/2;
  initialize();
  lastupdate_=0;
}

void NeighborList::initialize() {
  neighbors_.clear();
  for(unsigned int i=0; i<nallpairs_; ++i) {
    neighbors_.push_back(getIndexPair(i));
  }
}

std::vector<AtomNumber>& NeighborList::getFullAtomList() {
  return fullatomlist_;
}

std::pair<unsigned,unsigned> NeighborList::getIndexPair(unsigned ipair) {
  std::pair<unsigned,unsigned> index;
  if(twolists_ && do_pair_) {
    index=std::pair<unsigned,unsigned>(ipair,ipair+nlist0_);
  } else if (twolists_ && !do_pair_) {
    index=std::pair<unsigned,unsigned>(ipair/nlist1_,ipair%nlist1_+nlist0_);
  } else if (!twolists_) {
    unsigned ii = nallpairs_-1-ipair;
    unsigned  K = unsigned(std::floor((std::sqrt(double(8*ii+1))+1)/2));
    unsigned jj = ii-K*(K-1)/2;
    index=std::pair<unsigned,unsigned>(nlist0_-1-K,nlist0_-1-jj);
  }
  return index;
}

void NeighborList::update(const std::vector<Vector>& positions) {
  neighbors_.clear();
  const double d2=distance_*distance_;
  // check if positions array has the correct length
  plumed_assert(positions.size()==fullatomlist_.size());

  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  unsigned nt=OpenMP::getNumThreads();
  if(serial_) {
    stride=1;
    rank=0;
    nt=1;
  }
  std::vector<unsigned> local_flat_nl;

  #pragma omp parallel num_threads(nt)
  {
    std::vector<unsigned> private_flat_nl;
    #pragma omp for nowait
    for(unsigned int i=rank; i<nallpairs_; i+=stride) {
      std::pair<unsigned,unsigned> index=getIndexPair(i);
      unsigned index0=index.first;
      unsigned index1=index.second;
      Vector distance;
      if(do_pbc_) {
        distance=pbc_->distance(positions[index0],positions[index1]);
      } else {
        distance=delta(positions[index0],positions[index1]);
      }
      double value=modulo2(distance);
      if(value<=d2) {
        private_flat_nl.push_back(index0);
        private_flat_nl.push_back(index1);
      }
    }
    #pragma omp critical
    local_flat_nl.insert(local_flat_nl.end(), private_flat_nl.begin(), private_flat_nl.end());
  }

  // find total dimension of neighborlist
  std::vector <int> local_nl_size(stride, 0);
  local_nl_size[rank] = local_flat_nl.size();
  if(!serial_) comm.Sum(&local_nl_size[0], stride);
  int tot_size = std::accumulate(local_nl_size.begin(), local_nl_size.end(), 0);
  if(tot_size==0) {setRequestList(); return;}
  // merge
  std::vector<unsigned> merge_nl(tot_size, 0);
  // calculate vector of displacement
  std::vector<int> disp(stride);
  disp[0] = 0;
  int rank_size = 0;
  for(unsigned i=0; i<stride-1; ++i) {
    rank_size += local_nl_size[i];
    disp[i+1] = rank_size;
  }
  // Allgather neighbor list
  if(comm.initialized()&&!serial_) comm.Allgatherv((!local_flat_nl.empty()?&local_flat_nl[0]:NULL), local_nl_size[rank], &merge_nl[0], &local_nl_size[0], &disp[0]);
  else merge_nl = local_flat_nl;
  // resize neighbor stuff
  neighbors_.resize(tot_size/2);
  for(unsigned i=0; i<tot_size/2; i++) {
    unsigned j=2*i;
    neighbors_[i] = std::make_pair(merge_nl[j],merge_nl[j+1]);
  }

  setRequestList();
}

void NeighborList::setRequestList() {
  requestlist_.clear();
  for(unsigned int i=0; i<size(); ++i) {
    requestlist_.push_back(fullatomlist_[neighbors_[i].first]);
    requestlist_.push_back(fullatomlist_[neighbors_[i].second]);
  }
  Tools::removeDuplicates(requestlist_);
  reduced=false;
}

std::vector<AtomNumber>& NeighborList::getReducedAtomList() {
  if(!reduced)for(unsigned int i=0; i<size(); ++i) {
      unsigned newindex0=0,newindex1=0;
      AtomNumber index0=fullatomlist_[neighbors_[i].first];
      AtomNumber index1=fullatomlist_[neighbors_[i].second];
// I exploit the fact that requestlist_ is an ordered vector
      auto p = std::find(requestlist_.begin(), requestlist_.end(), index0); plumed_dbg_assert(p!=requestlist_.end()); newindex0=p-requestlist_.begin();
      p = std::find(requestlist_.begin(), requestlist_.end(), index1); plumed_dbg_assert(p!=requestlist_.end()); newindex1=p-requestlist_.begin();
      neighbors_[i]=std::pair<unsigned,unsigned>(newindex0,newindex1);
    }
  reduced=true;
  return requestlist_;
}

unsigned NeighborList::getStride() const {
  return stride_;
}

unsigned NeighborList::getLastUpdate() const {
  return lastupdate_;
}

void NeighborList::setLastUpdate(unsigned step) {
  lastupdate_=step;
}

unsigned NeighborList::size() const {
  return neighbors_.size();
}

std::pair<unsigned,unsigned> NeighborList::getClosePair(unsigned i) const {
  return neighbors_[i];
}

std::pair<AtomNumber,AtomNumber> NeighborList::getClosePairAtomNumber(unsigned i) const {
  std::pair<AtomNumber,AtomNumber> Aneigh;
  Aneigh=std::pair<AtomNumber,AtomNumber>(fullatomlist_[neighbors_[i].first],fullatomlist_[neighbors_[i].second]);
  return Aneigh;
}

std::vector<unsigned> NeighborList::getNeighbors(unsigned index) {
  std::vector<unsigned> neighbors;
  for(unsigned int i=0; i<size(); ++i) {
    if(neighbors_[i].first==index)  neighbors.push_back(neighbors_[i].second);
    if(neighbors_[i].second==index) neighbors.push_back(neighbors_[i].first);
  }
  return neighbors;
}

}
