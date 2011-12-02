#ifndef __PLUMED_ActionAtomistic_h
#define __PLUMED_ActionAtomistic_h

#include "ActionWithExternalArguments.h"
#include "Atoms.h"
#include "AtomicNeighbourList.h"
#include "PlumedMain.h"
#include <vector>
#include <set>
#include "Pbc.h"

namespace PLMD {

/// Action which can access to atomistic data
class ActionAtomistic : public ActionWithExternalArguments {

  friend class Atoms;
  friend class AtomicNeighbourList;
  friend class ColvarModifier;
  std::string atomGroupName;
  bool pbcOn;
  FILE* forcefile;

// Stuff for atoms
  std::vector<bool>     skips;
  std::vector<Vector>   positions;        // positions of the needed atoms
  Tensor                box;
  Pbc                   pbc;
  Tensor                virial;
  std::vector<double>   masses;
  std::vector<double>   charges;
// Stuff for neighbour lists
  unsigned nliststyle;
  int updateFreq;
  unsigned lastUpdate;
  double   nl_cut;
  std::vector<AtomicNeighbourList> nlists;

  bool                  lockRequestAtoms; // forbid changes to request atoms
protected:
/// Request an array of atoms.
/// This method is used to ask for a list of atoms. Atoms
/// should be asked for by number. If this routine is called
/// during the simulation, atoms will be available at the next step
/// MAYBE WE HAVE TO FIND SOMETHING MORE CLEAR FOR DYNAMIC
/// LISTS OF ATOMS
/// Read in actionAtomistics keywords
  void readActionAtomistic();
//  void readActionAtomistic( int& maxatoms, unsigned& maxgroups );
/// Recover the backbone atoms from the backbone keyword
  bool readBackboneAtoms( const std::string& type, std::vector< std::vector<unsigned> >& backbone );
/// Get position of i-th atom
  const Vector & getPositions(int)const;
/// Get the separation between two atoms
  Vector getSeparation(unsigned i, unsigned j) const;
/// Get the separation between two atoms (specified as vectors)
  Vector getSeparation(const Vector& v1, const Vector& v2) const;
/// Get the box
  const Tensor & getBox() const;
/// Get mass of i-th atom
  double getMasses(int i) const;
/// Get charge of i-th atom
  double getCharges(int i) const;
/// Get a reference to forces array
//  std::vector<Vector> & modifyForces();
/// Get a reference to virial array
//  Tensor & modifyVirial();
/// Get number of available atoms
  unsigned getNumberOfAtoms() const; 
/// Get the number of neighbour lists
  unsigned getNumberOfColvars() const;
/// Get the absolute index of an atom
  AtomNumber getAbsoluteIndex(int i)const;
/// Parse a list of atoms
  bool parseAtomList(const std::string&key, const unsigned num );
/// Print the full list of atoms to the log file
  void printAllAtoms( const std::string report );
/// Apply forces to the atoms
  void applyForces( const std::vector<Vector>& forces, const Tensor& virial );
/// Retrieve the set of skips for this quantity
  void retrieveSkips( std::vector<bool>& s ) const ;
/// Setup the style for the neighbour lists
  void setNeighbourListStyle( const std::string style );
/// Add a new neighbour list to the set of neighbour lists ( i.e. add a new colvar )
  void addColvar( AtomicNeighbourList& newlist );
/// Setup the neighbour lists
  void setupNeighbourList( const std::vector< std::pair<unsigned, unsigned> >& nlist_template );
/// This checks that a neighbour list is of a specific type
  bool checkNeighbourListType( const unsigned& req ) const ;
/// This checks that the neighbour lists have been setup
  void checkNeighbourLists() const ;
/// This checks that we are not declaring neighbour lists when they are not allowed
  void checkForBadNeighbourLists() const ;
/// This retrieves that atoms from the neighbour list that are currently close enough together to matter
  bool getColvarAtoms( const unsigned& ilist, std::vector<unsigned>& indexes ) const ;
public:

// virtual functions:

  ActionAtomistic(const ActionOptions&ao);
  ~ActionAtomistic();

//  void clearOutputForces();
  void retrieveData();

  void prepare();
  void updateNeighbourLists();
  void calculateNumericalDerivatives();
  void lockRequests();
  void unlockRequests();
  void resetSkips();
  void requireAtoms( const std::vector<bool>& required_atoms );
};

inline
void ActionAtomistic::resetSkips(){
  for(unsigned i=0;i<skips.size();++i) skips[i]=false;
  plumed.getAtoms().updateSkipsForGroup( atomGroupName, skips );
} 

inline
void ActionAtomistic::requireAtoms( const std::vector<bool>& required_atoms ) {
  assert( required_atoms.size()==skips.size() );
  for(unsigned i=0;i<skips.size();++i){ if( !required_atoms[i] ) skips[i]=true; }
  plumed.getAtoms().updateSkipsForGroup( atomGroupName, skips );
}

inline
const Tensor & ActionAtomistic::getBox() const {
  return box;
}

inline
unsigned ActionAtomistic::getNumberOfAtoms() const {
  return positions.size();
}

inline
unsigned ActionAtomistic::getNumberOfColvars() const {
  return nlists.size();
}

inline
const Vector & ActionAtomistic::getPositions(int i) const {
  assert(!skips[i]);
  return positions[i];
}

inline
Vector ActionAtomistic::getSeparation(unsigned i, unsigned j) const {
  assert( !skips[i] && !skips[j] );
  if ( pbcOn ) return pbc.distance( positions[i], positions[j] );
  return delta( positions[i], positions[j] );
}

inline
Vector ActionAtomistic::getSeparation(const Vector& v1, const Vector& v2) const {
  if ( pbcOn ) return pbc.distance( v1, v2 );
  return delta( v1, v2 );
}

inline
double ActionAtomistic::getMasses(int i) const {
  assert(!skips[i]);
  return masses[i];
}

inline
double ActionAtomistic::getCharges(int i) const {
  assert(!skips[i]);
  return charges[i];
}

inline
bool ActionAtomistic::getColvarAtoms( const unsigned& ilist, std::vector<unsigned>& indexes ) const {
  assert( nlists[ilist].active );

  if( indexes.size()!=nlists[ilist].nactive ) indexes.resize( nlists[ilist].nactive );
  for(unsigned i=0;i<nlists[ilist].all_atoms.size();++i){
      if( !skips[nlists[ilist].all_atoms[i]] ) indexes[i]=nlists[ilist].all_atoms[i];
  }
  return true;
}

inline
bool ActionAtomistic::checkNeighbourListType( const unsigned& req ) const {
  return ( nliststyle==req );
}

//inline
//std::vector<Vector> & ActionAtomistic::modifyForces() {
//  return forces;
//}

//inline
//Tensor & ActionAtomistic::modifyVirial() {
//  return virial;
//}

//inline
//void ActionAtomistic::clearOutputForces() {
//  for(unsigned i=0;i<forces.size();++i) forces[i].clear();
//}

inline
void ActionAtomistic::lockRequests() {
  lockRequestAtoms=true;
}

inline
void ActionAtomistic::unlockRequests() {
  lockRequestAtoms=false;
}

}
#endif
