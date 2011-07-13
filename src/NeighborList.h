#ifndef __PLUMED_NeighborList_h
#define __PLUMED_NeighborList_h

#include <vector>
#include "Vector.h"
#include "Pbc.h"
#include "AtomNumber.h"

/// NeighborList 
class NeighborList  
{
  bool do_pair_,twolists_;
  PLMD::Pbc *pbc_;
  std::vector<PLMD::AtomNumber> fullatomlist_,requestlist_;
  std::vector< std::pair<unsigned,unsigned> > neighbors_;
  double distance_;
  unsigned stride_,nlist0_,nlist1_,nallpairs_;
  void initialize();
  std::pair<unsigned,unsigned> getIndexPair(unsigned ipair);
  template<typename T> 
   void removeDuplicates(std::vector<T>& vec);
  std::vector<unsigned> getRequestList();
public:
  NeighborList(std::vector<unsigned> list0, std::vector<unsigned> list1,
               bool do_pair, PLMD::Pbc *pbc, double distance, unsigned stride);
  NeighborList(std::vector<unsigned> list0, PLMD::Pbc *pbc,
               double distance, unsigned stride);
                             
  std::vector<PLMD::AtomNumber> & getFullAtomList() const;
  std::vector<PLMD::AtomNumber> & getReducedAtomList(std::vector<PLMD::Vector> positions);
  void update();
  unsigned getStride() const;
  unsigned size() const;
  const std::pair<unsigned,unsigned> & operator[] (unsigned i)const;
  std::vector<unsigned> getNeighbors(unsigned i);
  ~NeighborList(){};
};

#endif