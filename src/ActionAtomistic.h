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

  double                sigma;           // typical variable scale
  Atoms::Request*       atomRequest;     // handler for request of atoms
  std::vector<int> indexes;         // the set of needed atoms
  std::vector<Vector>   positions;       // positions of the needed atoms
  Tensor                box;
  Tensor                virial;
  std::vector<Vector>   forces;          // forces on the needed atoms
  std::vector<double>   masses;
  std::vector<double>   charges;

protected:
  void requestAtoms(const std::vector<int> & a);
  const Vector & getPositions(int)const;
  const Tensor & getBox()const;
  const std::vector<Vector> & getPositions()const;
  double getMasses(int i)const;
  double getCharges(int i)const;
  std::vector<Vector> & modifyForces();
  Tensor & modifyVirial();
  int getNatoms(){return indexes.size();};
  Vector pbcDistance(const Vector&,const Vector&)const;

public:

// virtual functions:

  ActionAtomistic(const ActionOptions&ao);
  ~ActionAtomistic();

  void activate();
  void deactivate();

  void clearOutputForces();
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
