#ifndef __PLUMED_ActionWithVirtualAtom_h
#define __PLUMED_ActionWithVirtualAtom_h

#include "ActionAtomistic.h"
#include "AtomNumber.h"
#include "Vector.h"
#include "Tensor.h"

namespace PLMD{

/// Class to add a single virtual atom to the system.
/// (it might be extended to add multiple virtual atoms).
class ActionWithVirtualAtom : public ActionAtomistic {
private:
/// The index of the virtual atom in the Atoms class
  unsigned index;
/// The derivatives for the virtual atom wrt to the positions of the atoms
  std::vector<Tensor> derivatives;
/// The forces on each atom
  std::vector<Vector> f;
protected:
/// Set position of the virtual atom
  void setPosition(const Vector &);
/// Set its mass
  void setMass(double);
/// Set its charge
  void setCharge(double);
/// Request atoms on which the calculation depends
//  void requestAtoms(const std::vector<AtomNumber> & a);
/// Set the derivatives of virtual atom coordinate wrt atoms on which it dependes
  void setAtomsDerivatives(const std::vector<Tensor> &d);
/// Read everything in higher levels for virtual atoms
  void readActionWithVirtualAtom();
public:
/// Return the atom id of the corresponding virtual atom
  AtomNumber getIndex() const;
  ActionWithVirtualAtom(const ActionOptions&ao);
  ~ActionWithVirtualAtom();
  void interpretGroupsKeyword( const unsigned& natoms, const std::string& atomGroupName, const std::vector<std::vector<unsigned> >& groups );
  void interpretAtomsKeyword( const std::vector<std::vector<unsigned> >& flist );
  void apply();
};

inline
AtomNumber ActionWithVirtualAtom::getIndex()const{
  return AtomNumber::index(index);
}

inline
void ActionWithVirtualAtom::setPosition(const Vector & pos){
  plumed.getAtoms().positions[index]=pos;
}

inline
void ActionWithVirtualAtom::setMass(double m){
  plumed.getAtoms().masses[index]=m;
}

inline
void ActionWithVirtualAtom::setCharge(double c){
  plumed.getAtoms().charges[index]=c;
}

inline
void ActionWithVirtualAtom::setAtomsDerivatives(const std::vector<Tensor> &d){
  derivatives=d;
}

}

#endif
