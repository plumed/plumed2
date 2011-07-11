#ifndef __PLUMED_NeighborList_h
#define __PLUMED_NeighborList_h

#include <vector>
#include "Vector.h"
#include "Pbc.h"

/// NeighborList 
class NeighborList  
{
  std::vector<int> list0_,list1_,indexes_; 
  std::vector< std::vector<int> > neighbors_;
  double distance_;
  PLMD::Pbc *pbc_;
  int stride_;
public:
  NeighborList(std::vector<int> list0, std::vector<int> list1,
               double distance, unsigned stride, PLMD::Pbc *pbc);
  std::vector<int> getFullList();
  std::vector<int> update(std::vector<PLMD::Vector> positions);
  std::vector<int> getNeighbors(int index);  
  int getStride() const;
  int getNumberOfAtoms() const;
  int getNumberOfNeighbors(int index);
  ~NeighborList(){};
};

#endif

