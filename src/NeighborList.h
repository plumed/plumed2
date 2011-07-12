#ifndef __PLUMED_NeighborList_h
#define __PLUMED_NeighborList_h

#include <vector>
#include "Vector.h"
#include "Pbc.h"

/// NeighborList 
class NeighborList  
{
  std::vector<unsigned> list0_,list1_,indexes_,newindexes_; 
  std::vector< std::vector<unsigned> > neighbors_;
  double distance_;
  PLMD::Pbc *pbc_;
  unsigned stride_, imax_;
  std::vector<unsigned> createIndexes(std::vector<unsigned> request);
public:
  NeighborList(std::vector<unsigned> list0, std::vector<unsigned> list1,
               double distance, unsigned stride, PLMD::Pbc *pbc);
  std::vector<unsigned> getFullList();
  std::vector<unsigned> prepareUpdate(std::vector<PLMD::Vector> positions);
  void finishUpdate();
  std::vector<unsigned> getNeighbors(unsigned index);  
  unsigned getStride() const;
  unsigned getNumberOfAtoms() const;
  unsigned getNumberOfNeighbors(unsigned index);
  ~NeighborList(){};
};

#endif

