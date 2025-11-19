/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "Exception.h"
#include "Vector.h"
#include "Pbc.h"
#include "AtomNumber.h"
#include "Communicator.h"
#include "OpenMP.h"
#include "Tools.h"
#include <string>
#include "LinkCells.h"
#include "View.h"
#include <vector>
#include <algorithm>
#include <numeric>
#include <cstdlib>

#pragma GCC diagnostic error "-Wswitch"

namespace PLMD {

NeighborList::NeighborList(const std::vector<AtomNumber>& list0,
                           const std::vector<AtomNumber>& list1,
                           const bool serial,
                           const bool do_pair,
                           const bool do_pbc,
                           const Pbc& pbc,
                           Communicator& cm,
                           const double distance,
                           const unsigned stride)
  : serial_(serial),
    do_pbc_(do_pbc),
    style_(do_pair ? NNStyle::Pair : NNStyle::TwoList),
    pbc_(&pbc),
    comm(cm),
    //copy-initialize fullatomlist_
    fullatomlist_(list0),
    distance_(distance),
    nlist0_(list0.size()),
    nlist1_(list1.size()),
    stride_(stride) {
  // store the rest of the atoms into fullatomlist_
  fullatomlist_.insert(fullatomlist_.end(),list1.begin(),list1.end());
  if(style_ != NNStyle::Pair) {
    nallpairs_=nlist0_*nlist1_;
  } else {
    plumed_assert(nlist0_==nlist1_)
        << "when using PAIR option, the two groups should have the same number"
        " of elements\n" << "the groups you specified have size "
        <<nlist0_<<" and "<<nlist1_;
    nallpairs_=nlist0_;
  }
  initialize();
}

NeighborList::NeighborList(const std::vector<AtomNumber>& list0,
                           const bool serial,
                           const bool do_pbc,
                           const Pbc& pbc,
                           Communicator& cm,
                           const double distance,
                           const unsigned stride)
  : serial_(serial),
    do_pbc_(do_pbc),
    style_(NNStyle::SingleList),
    pbc_(&pbc),
    comm(cm),
    //copy-initialize fullatomlist_
    fullatomlist_(list0),
    distance_(distance),
    nlist0_(list0.size()),
    nallpairs_(nlist0_*(nlist0_-1)/2),
    stride_(stride) {
  initialize();
}

NeighborList::~NeighborList()=default;

void NeighborList::initialize() {
  constexpr const char* envKey="PLUMED_IGNORE_NL_MEMORY_ERROR";
  if(!std::getenv(envKey)) {
    //blocking memory allocation on more than 10 GB of memory
    //A single list of more than 50000 atoms
    //two different lists of more than 35355 atoms each (lista*listb < max, see below)
    //that is more than 1250000000 pairs
    //each pairIDs occupies 64 bit (where unsigned are 32bit integers)
    //4294967296 is max(uint32)+1 and is more than 34 GB (correspond to a system of 65536 atoms)
    if(nallpairs_ > 1250000000 ) {
      const unsigned GB = sizeof(decltype(neighbors_)::value_type) * nallpairs_ / 1000000000;
      plumed_error() << "A NeighborList is trying to allocate "
                     + std::to_string( GB ) +" GB of data for the list of neighbors\n"
                     "You can skip this error by exporting \""+envKey+"\"";
    }
  }
  try {
    neighbors_.resize(nallpairs_);
  } catch (...) {
    plumed_error_nested() << "An error happened while allocating the neighbor "
                          "list, please decrease the number of atoms used";
  }
  //TODO: test if this is feasible for accelerating the loop
  //#pragma omp parallel for default(shared)
  for(unsigned int i=0; i<nallpairs_; ++i) {
    neighbors_[i]=getIndexPair(i);
  }
}

std::vector<AtomNumber>& NeighborList::getFullAtomList() {
  return fullatomlist_;
}

NeighborList::pairIDs NeighborList::getIndexPair(const unsigned ipair) {
  pairIDs index;
  switch (style_) {
  case NNStyle::Pair : {
    index=pairIDs(ipair,ipair+nlist0_);
  }
  break;
  case NNStyle::TwoList : {
    index=pairIDs(ipair/nlist1_,ipair%nlist1_+nlist0_);
  }
  break;
  case NNStyle::SingleList: {
    unsigned ii = nallpairs_-1-ipair;
    unsigned  K = unsigned(std::floor((std::sqrt(double(8*ii+1))+1)/2));
    unsigned jj = ii-K*(K-1)/2;
    index=pairIDs(nlist0_-1-K,nlist0_-1-jj);
  }
  }
  return index;
}

void NeighborList::update(const std::vector<Vector>& positions) {
  neighbors_.clear();
  // check if positions array has the correct length
  plumed_assert(positions.size()==fullatomlist_.size());

#ifdef _OPENMP
  //nt is unused if openmp is not declared
  const unsigned nt=(serial_)? 1 : OpenMP::getNumThreads();
#endif //_OPENMP
  if(useCellList_) {
    switch (style_) {
    case NNStyle::TwoList: {
      LinkCells cells(comm);
      cells.setCutoff(distance_);

      std::vector<unsigned> indexesForCells(fullatomlist_.size());
      std::iota(indexesForCells.begin(),indexesForCells.end(),0);
      cells.setupCells(make_const_view(positions),*pbc_);
      auto listA = cells.getCollection(View{positions.data(),nlist0_},
                                       View<const unsigned> {indexesForCells.data(),nlist0_});
      auto listB = cells.getCollection(View{positions.data()+nlist0_,nlist1_},
                                       View<const unsigned> {indexesForCells.data()+nlist0_,nlist1_});
      plumed_assert((listA.lcell_lists.size()+listB.lcell_lists.size()) == positions.size())
          << listA.lcell_lists.size()<<"+"<<listB.lcell_lists.size() <<"==" <<positions.size();
      //#pragma omp parallel num_threads(nt)
      std::vector<unsigned> cells_required(27);
      for(unsigned c =0; c < cells.getNumberOfCells(); ++c) {
        auto atomsInC= listA.getCellIndexes(c);
        if(atomsInC.size()>0) {
          auto cell = cells.findMyCell(c);
          unsigned ncells_required=0;
          cells.addRequiredCells(cell,ncells_required, cells_required,do_pbc_);
          for (auto A : atomsInC) {
            for (unsigned cb=0; cb <ncells_required ; ++cb) {
              for (auto B : listB.getCellIndexes(cells_required[cb])) {
                neighbors_.push_back({A,B});
              }
            }
          }
        }


      }
    }
    break;//*/
    case NNStyle::SingleList: {
      LinkCells cells(comm);
      cells.setCutoff(distance_);

      std::vector<unsigned> indexesForCells(fullatomlist_.size());
      std::iota(indexesForCells.begin(),indexesForCells.end(),0);
      cells.setupCells(make_const_view(positions),*pbc_);
      auto listA = cells.getCollection(positions,indexesForCells);
      //#pragma omp parallel num_threads(nt)
      std::vector<unsigned> cells_required(27);
      for(unsigned c =0; c < cells.getNumberOfCells(); ++c) {
        auto atomsInC= listA.getCellIndexes(c);
        if(atomsInC.size()>0) {
          auto cell = cells.findMyCell(c);
          unsigned ncells_required=0;
          cells.addRequiredCells(cell,ncells_required, cells_required,do_pbc_);
          for (auto A : atomsInC) {
            for (unsigned cb=0; cb <ncells_required ; ++cb) {
              for (auto B : listA.getCellIndexes(cells_required[cb])) {
                if (B>A) {
                  neighbors_.push_back({A,B});
                }
              }
            }
          }
        }


      }
    }
    break;
    case NNStyle::Pair:
      plumed_error() << "Cell list should not be active with a Pair NL";
    }
  } else {
    const double d2=distance_*distance_;
    // check if positions array has the correct length
    plumed_assert(positions.size()==fullatomlist_.size());

    const unsigned stride=(serial_)? 1 : comm.Get_size();
    const unsigned rank  =(serial_)? 0 : comm.Get_rank();
    const unsigned elementsPerRank = std::ceil(double(nallpairs_)/stride);
    const unsigned int start= rank*elementsPerRank;
    const unsigned int end = ((start + elementsPerRank)< nallpairs_)?(start + elementsPerRank): nallpairs_;
    std::vector<unsigned> local_flat_nl;

    #pragma omp parallel num_threads(nt)
    {
      std::vector<unsigned> private_flat_nl;
      #pragma omp for nowait
      for(unsigned int i=start; i<end; ++i) {
        auto [index0, index1 ] = getIndexPair(i);
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
      local_flat_nl.insert(local_flat_nl.end(),
                           private_flat_nl.begin(),
                           private_flat_nl.end());
    }

    // find total dimension of neighborlist
    std::vector <int> local_nl_size(stride, 0);
    local_nl_size[rank] = local_flat_nl.size();
    if(!serial_) {
      comm.Sum(&local_nl_size[0], stride);
    }
    int tot_size = std::accumulate(local_nl_size.begin(), local_nl_size.end(), 0);
    if(tot_size!=0) {
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
      if(comm.initialized()&&!serial_) {
        comm.Allgatherv((!local_flat_nl.empty()?&local_flat_nl[0]:NULL),
                        local_nl_size[rank],
                        &merge_nl[0],
                        &local_nl_size[0],
                        &disp[0]);
      } else {
        merge_nl = local_flat_nl;
      }
      // resize neighbor stuff
      neighbors_.resize(tot_size/2);
      for(int i=0; i<tot_size/2; i++) {
        unsigned j=2*i;
        neighbors_[i] = std::make_pair(merge_nl[j],merge_nl[j+1]);
      }
    }
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
  if(!reduced) {
    for(unsigned int i=0; i<size(); ++i) {
      AtomNumber index0=fullatomlist_[neighbors_[i].first];
      AtomNumber index1=fullatomlist_[neighbors_[i].second];
// I exploit the fact that requestlist_ is an ordered vector
// And I assume that index0 and index1 actually exists in the requestlist_ (see setRequestList())
// so I can use lower_bond that uses binary seach instead of find
      plumed_dbg_assert(std::is_sorted(requestlist_.begin(),requestlist_.end()));
      auto p = std::lower_bound(requestlist_.begin(), requestlist_.end(), index0);
      plumed_dbg_assert(p!=requestlist_.end());
      unsigned newindex0=p-requestlist_.begin();
      p = std::lower_bound(requestlist_.begin(), requestlist_.end(), index1);
      plumed_dbg_assert(p!=requestlist_.end());
      unsigned newindex1=p-requestlist_.begin();
      neighbors_[i]=pairIDs(newindex0,newindex1);
    }
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

void NeighborList::setLastUpdate(const unsigned step) {
  lastupdate_=step;
}

unsigned NeighborList::size() const {
  return neighbors_.size();
}

NeighborList::pairIDs NeighborList::getClosePair(const unsigned i) const {
  return neighbors_[i];
}

NeighborList::pairAtomNumbers
NeighborList::getClosePairAtomNumber(const unsigned i) const {
  pairAtomNumbers Aneigh=pairAtomNumbers(
                           fullatomlist_[neighbors_[i].first],
                           fullatomlist_[neighbors_[i].second]);
  return Aneigh;
}

std::vector<unsigned> NeighborList::getNeighbors(const unsigned index) const {
  std::vector<unsigned> neighbors;
  for(unsigned int i=0; i<size(); ++i) {
    if(neighbors_[i].first==index) {
      neighbors.push_back(neighbors_[i].second);
    }
    if(neighbors_[i].second==index) {
      neighbors.push_back(neighbors_[i].first);
    }
  }
  return neighbors;
}

} // namespace PLMD
