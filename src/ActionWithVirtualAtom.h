#ifndef __PLUMED_ActionWithVirtualAtom_h
#define __PLUMED_ActionWithVirtualAtom_h

#include "ActionAtomistic.h"
#include "AtomNumber.h"
#include "Vector.h"
#include "Tensor.h"
#include "PlumedMain.h"
#include "Atoms.h"

namespace PLMD{

//+DEVELDOC INHERIT ActionWithVirtualAtom
/**
Inherit from here if you are calculating the position of a virtual atom (eg a center of mass)
*/
//+ENDDEVELDOC

/// Class to add a single virtual atom to the system.
/// (it might be extended to add multiple virtual atoms).
class ActionWithVirtualAtom:
  public ActionAtomistic
{
  AtomNumber index;
  std::vector<Tensor> derivatives;
  std::map<AtomNumber,Tensor> gradients;
  void apply();
protected:
/// Set position of the virtual atom
  void setPosition(const Vector &);
/// Set its mass
  void setMass(double);
/// Set its charge
  void setCharge(double);
/// Request atoms on which the calculation depends
  void requestAtoms(const std::vector<AtomNumber> & a);
/// Set the derivatives of virtual atom coordinate wrt atoms on which it dependes
  void setAtomsDerivatives(const std::vector<Tensor> &d);
public:
  void setGradients();
  const std::map<AtomNumber,Tensor> & getGradients()const;
/// Return the atom id of the corresponding virtual atom
  AtomNumber getIndex()const;
  ActionWithVirtualAtom(const ActionOptions&ao);
  ~ActionWithVirtualAtom();
  static void registerKeywords(Keywords& keys);
  void setGradientsIfNeeded();
};

inline
AtomNumber ActionWithVirtualAtom::getIndex()const{
  return index;
}

inline
void ActionWithVirtualAtom::setPosition(const Vector & pos){
  plumed.getAtoms().positions[index.index()]=pos;
}

inline
void ActionWithVirtualAtom::setMass(double m){
  plumed.getAtoms().masses[index.index()]=m;
}

inline
void ActionWithVirtualAtom::setCharge(double c){
  plumed.getAtoms().charges[index.index()]=c;
}

inline
void ActionWithVirtualAtom::setAtomsDerivatives(const std::vector<Tensor> &d){
  derivatives=d;
}

inline
const std::map<AtomNumber,Tensor> & ActionWithVirtualAtom::getGradients()const{
  return gradients;
}

}

#endif
