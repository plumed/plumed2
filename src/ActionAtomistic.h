#ifndef __PLUMED_ActionAtomistic_h
#define __PLUMED_ActionAtomistic_h

#include "ActionWithValue.h"
#include "Atoms.h"
#include "PlumedMain.h"
#include <vector>
#include <set>
#include "Pbc.h"

namespace PLMD {

/// Action which can access to atomistic data
class ActionAtomistic : public ActionWithValue {

  friend class Atoms;

  std::vector<int>      indexes;         // the set of needed atoms
  std::set<int>         unique;
  std::vector<Vector>   positions;       // positions of the needed atoms
  double                energy;
  Tensor                box;
  Pbc                   pbc;
  Tensor                virial;
  std::vector<double>   masses;
  std::vector<double>   charges;

  std::vector<Vector>   forces;          // forces on the needed atoms
  double                forceOnEnergy;

  bool                  lockRequestAtoms; // forbid changes to request atoms

protected:
/// Request an array of atoms.
/// This method is used to ask for a list of atoms. Atoms
/// should be asked for by number. If this routine is called
/// during the simulation, atoms will be available at the next step
/// MAYBE WE HAVE TO FIND SOMETHING MORE CLEAR FOR DYNAMIC
/// LISTS OF ATOMS
  void requestAtoms(const std::vector<AtomNumber> & a);
/// Get position of i-th atom
  const Vector & getPositions(int)const;
/// Get position of i-th atom
  const Tensor & getBox()const;
/// Get the array of all positions
  const std::vector<Vector> & getPositions()const;
/// Get energy
  const double & getEnergy()const;
/// Get mass of i-th atom
  double getMasses(int i)const;
/// Get charge of i-th atom
  double getCharges(int i)const;
/// Get a reference to forces array
  std::vector<Vector> & modifyForces();
/// Get a reference to virial array
  Tensor & modifyVirial();
/// Get a reference to force on energy
  double & modifyForceOnEnergy();
/// Get number of available atoms
  unsigned getNatoms()const{return indexes.size();};
/// Compute the pbc distance between two positions
  Vector pbcDistance(const Vector&,const Vector&)const;
/// Get the absolute index of an atom
  AtomNumber getAbsoluteIndex(int i)const;
/// Parse a list of atoms
  void parseAtomList(const std::string&key,std::vector<AtomNumber> &t);
/// Get reference to Pbc
  const Pbc & getPbc() const;
public:

// virtual functions:

  ActionAtomistic(const ActionOptions&ao);
  ~ActionAtomistic();


  void clearOutputForces();

  void   calculateNumericalDerivatives();

  void retrieveAtoms();
  void applyForces();
  void lockRequests();
  void unlockRequests();
};

inline
const Vector & ActionAtomistic::getPositions(int i)const{
  return positions[i];
}

inline
double ActionAtomistic::getMasses(int i)const{
  return masses[i];
}

inline
double ActionAtomistic::getCharges(int i)const{
  return charges[i];
}

inline
AtomNumber ActionAtomistic::getAbsoluteIndex(int i)const{
  return AtomNumber::index(indexes[i]);
}

inline
const std::vector<Vector> & ActionAtomistic::getPositions()const{
  return positions;
}

inline
const double & ActionAtomistic::getEnergy()const{
  return energy;
}

inline
const Tensor & ActionAtomistic::getBox()const{
  return box;
}

inline
std::vector<Vector> & ActionAtomistic::modifyForces(){
  return forces;
}

inline
Tensor & ActionAtomistic::modifyVirial(){
  return virial;
}

inline
void ActionAtomistic::clearOutputForces(){
  for(unsigned i=0;i<forces.size();++i)forces[i].clear();
  forceOnEnergy=0.0;
}


inline
double & ActionAtomistic::modifyForceOnEnergy(){
  return forceOnEnergy;
}

inline
const Pbc & ActionAtomistic::getPbc() const{
 return pbc;
}

inline
void ActionAtomistic::lockRequests(){
  lockRequestAtoms=true;
}

inline
void ActionAtomistic::unlockRequests(){
  lockRequestAtoms=false;
}

}

#endif
