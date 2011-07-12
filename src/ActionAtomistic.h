#ifndef __PLUMED_ActionAtomistic_h
#define __PLUMED_ActionAtomistic_h

#include <vector>
#include <string>
#include "Action.h"
#include "Atoms.h"

namespace PLMD {

/// Action which can access to atomistic data
class ActionAtomistic :
  virtual public Action
  {

  Atoms::Request*       atomRequest;     // handler for request of atoms
  std::vector<int> indexes;         // the set of needed atoms
  std::vector<Vector>   positions;       // positions of the needed atoms
  Tensor                box;
  Tensor                virial;
  std::vector<Vector>   forces;          // forces on the needed atoms
  std::vector<double>   masses;
  std::vector<double>   charges;

protected:
/// Request an array of atoms.
/// This method is used to ask for a list of atoms. Atoms
/// should be asked for by number. If this routine is called
/// during the simulation, atoms will be available at the next step
/// MAYBE WE HAVE TO FIND SOMETHING MORE CLEAR FOR DYNAMIC
/// LISTS OF ATOMS
  void requestAtoms(const std::vector<int> & a);
/// Get position of i-th atom
  const Vector & getPositions(int)const;
/// Get position of i-th atom
  const Tensor & getBox()const;
/// Get the array of all positions
  const std::vector<Vector> & getPositions()const;
/// Get mass of i-th atom
  double getMasses(int i)const;
/// Get charge of i-th atom
  double getCharges(int i)const;
/// Get a reference to forces array
  std::vector<Vector> & modifyForces();
/// Get a reference to virial array
  Tensor & modifyVirial();
/// Get number of available atoms
  int getNatoms()const{return indexes.size();};
/// Compute the pbc distance between two positions
  Vector pbcDistance(const Vector&,const Vector&)const;
/// Get the absolute index of an atom
  int getAbsoluteIndex(int i)const;

public:

// virtual functions:

  ActionAtomistic(const ActionOptions&ao);
  ~ActionAtomistic();

  void activate();
  void deactivate();

  void clearOutputForces();

  void   calculateNumericalDerivatives();
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
int ActionAtomistic::getAbsoluteIndex(int i)const{
  return indexes[i];
}

inline
const std::vector<Vector> & ActionAtomistic::getPositions()const{
  return positions;
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
  for(unsigned i=0;i<forces.size();++i){
    forces[i][0]=0.0;
    forces[i][1]=0.0;
    forces[i][2]=0.0;
  }
}



}

#endif
