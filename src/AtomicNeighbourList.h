#ifndef __PLUMED_AtomicNeighborList_h
#define __PLUMED_AtomicNeighborList_h

#include "Vector.h"
#include <cassert>
#include <vector>

namespace PLMD {

/// An atomic neighbour list serves two purposes:
/// (1) It stores a list of atoms in a colvar from which a colvar is created (all_atoms)
/// (2) It stores a list of pairs of atoms (neighbours) which must either be within a certain distance (rcut)

class AtomicNeighbourList {
friend class ActionAtomistic;
private:
/// Pointer to the action for which this is a neighbour list  
  ActionAtomistic* action;
/// The style of the neighbour list 1=none, 2=don't calculate colvar unless all pairs are within rcut and 3=ignore atoms that are only part of pairs that are more than rcut appart
  unsigned style;
/// Are we calculating the colvar for this action
  bool active;
/// The cutoff
  double rcut;
/// The number of atoms that are currently active
  unsigned nactive;
/// The list of all atoms
  std::vector<unsigned> all_atoms;
/// Pairs we are skipping
  std::vector<unsigned> skipto;
/// The list of pairs of atoms of interest
  std::vector< std::pair<unsigned,unsigned> > neighbours;
public:
  AtomicNeighbourList(ActionAtomistic* act);
/// Clear everything in the neighbour list (used during input)
  void clear();
/// Add an atom to the all_atoms list
  void addAtom( const unsigned& atom1 );
/// Add a pair of atoms to neighbours
  void addPair( const unsigned& atom1, const unsigned& atom2 );
/// Complete the setup of the neighbour list
  void completeSetup( const unsigned& ltype, const double& nl_cut );
/// Update the neighbour list
  void update( std::vector<bool>& atom_skips );
/// Are we calculating the neighbour list at this step
  bool mustCalculate() const; 
};

inline
bool AtomicNeighbourList::mustCalculate() const {
  return active;
}

}

#endif
