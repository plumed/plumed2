#ifndef __PLUMED_NeighborList_h
#define __PLUMED_NeighborList_h

#include <vector>
#include "Vector.h"
#include "Pbc.h"
#include "AtomNumber.h"

namespace PLMD{

/// NeighborList 
class NeighborList  
{
  bool do_pair_,do_pbc_,twolists_;
  PLMD::Pbc pbc_;
  std::vector<PLMD::AtomNumber> fullatomlist_,requestlist_;
  std::vector< std::pair<unsigned,unsigned> > neighbors_;
  double distance_;
  unsigned stride_,nlist0_,nlist1_,nallpairs_;
  void initialize();
  std::pair<unsigned,unsigned> getIndexPair(unsigned ipair);
  template<typename T> 
   void removeDuplicates(std::vector<T>& vec);
  void setRequestList();
public:
  NeighborList(const std::vector<PLMD::AtomNumber>& list0,
               const std::vector<PLMD::AtomNumber>& list1,
               const bool& do_pair, const bool& do_pbc, const PLMD::Pbc& pbc, 
               const double& distance=1.0e+30, const unsigned& stride=0);
  NeighborList(const std::vector<PLMD::AtomNumber>& list0, const bool& do_pbc,
               const PLMD::Pbc& pbc, const double& distance=1.0e+30,
               const unsigned& stride=0);
                             
  std::vector<PLMD::AtomNumber>& getFullAtomList();
  std::vector<PLMD::AtomNumber>& getReducedAtomList();
  void update(const std::vector<PLMD::Vector>& positions);
  unsigned getStride() const;
  unsigned size() const;
  const std::pair<unsigned,unsigned> & operator[] (unsigned i)const;
  std::vector<unsigned> getNeighbors(unsigned i);
  ~NeighborList(){};
};

}

#endif