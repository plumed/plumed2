#ifndef __PLUMED_ActionAtomistic_h
#define __PLUMED_ActionAtomistic_h

#include "ActionWithExternalArguments.h"
#include "Atoms.h"
#include "PlumedMain.h"
#include <vector>
#include <set>
#include "Pbc.h"

namespace PLMD {

/// Action which can access to atomistic data
class ActionAtomistic : public ActionWithExternalArguments {

  friend class Atoms;

  std::string atomGroupName;
  bool pbcOn;

// Stuff for atoms
  std::vector<bool>     skips;
  std::vector<Vector>   positions;        // positions of the needed atoms
  Tensor                box;
  Pbc                   pbc;
  Tensor                virial;
  std::vector<double>   masses;
  std::vector<double>   charges;

  std::vector<Vector>   forces;           // forces on the needed atoms

  bool                  lockRequestAtoms; // forbid changes to request atoms

// Stuff for dynamic content
  unsigned              updateFreq;
  unsigned              lastUpdate;
  double                nl_cut;
protected:
/// Request an array of atoms.
/// This method is used to ask for a list of atoms. Atoms
/// should be asked for by number. If this routine is called
/// during the simulation, atoms will be available at the next step
/// MAYBE WE HAVE TO FIND SOMETHING MORE CLEAR FOR DYNAMIC
/// LISTS OF ATOMS
/// Read in actionAtomistics keywords
  void readActionAtomistic( int& maxatoms, unsigned& maxgroups );
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
  std::vector<Vector> & modifyForces();
/// Get a reference to virial array
  Tensor & modifyVirial();
/// Get number of available atoms
  unsigned getNumberOfAtoms() const; 
/// Get the absolute index of an atom
  AtomNumber getAbsoluteIndex(int i)const;
/// Parse a list of atoms
  void parseAtomList(const std::string&key,std::vector<AtomNumber> &t);
/// Apply forces to the atoms
  void applyForces( const std::vector<Vector>& forces, const Tensor& virial );
/// Is it time to update any dynamic content in the "colvar"
  bool updateTime() const ;
/// Update the requests for the dynamic atom lists
  void updateDynamicAtoms();
/// Check if we are using updates
  bool updateIsOn() const;
/// Check if we have dynamic groups
  bool usingDynamicGroups() const;
public:

// virtual functions:

  ActionAtomistic(const ActionOptions&ao);
  ~ActionAtomistic();

  virtual void clearOutputForces();
  virtual void retrieveData();

  virtual void interpretGroupsKeyword( const unsigned& natoms, const std::string& atomGroupName, const std::vector<std::vector<unsigned> >& groups )=0;
  virtual void interpretAtomsKeyword( const std::vector<std::vector<unsigned> >& flist )=0;
  virtual void updateNeighbourList( const double& cutoff, std::vector<bool>& at_skips  )=0; 
  void prepare();
  void calculateNumericalDerivatives();
  void lockRequests();
  void unlockRequests();
};

inline
bool ActionAtomistic::updateTime() const {
  return ( updateFreq>0 && (getStep()-lastUpdate)>=updateFreq );
}

inline
unsigned ActionAtomistic::getNumberOfAtoms() const {
  return positions.size();
}

inline
const Vector & ActionAtomistic::getPositions(int i) const {
  assert(!skips[i]);
  return positions[i];
}

inline
Vector ActionAtomistic::getSeparation(unsigned i, unsigned j) const {
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
std::vector<Vector> & ActionAtomistic::modifyForces() {
  return forces;
}

inline
Tensor & ActionAtomistic::modifyVirial() {
  return virial;
}

inline
void ActionAtomistic::clearOutputForces() {
  for(unsigned i=0;i<forces.size();++i) forces[i].clear();
}

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
