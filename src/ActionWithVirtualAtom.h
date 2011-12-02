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
/// Set the derivatives of virtual atom coordinate wrt atoms on which it dependes
  void setAtomsDerivatives(const std::vector<Tensor> &d);
/// Read everything in higher levels for virtual atoms
  void readActionWithVirtualAtom();
public:
/// Return the atom id of the corresponding virtual atom
  AtomNumber getIndex() const;
  ActionWithVirtualAtom(const ActionOptions&ao);
  ~ActionWithVirtualAtom();
  void apply();
};

inline
AtomNumber ActionWithVirtualAtom::getIndex()const{
  return AtomNumber::index(index);
}

inline
void ActionWithVirtualAtom::setPosition(const Vector & pos){
  plumed.getAtoms().positions[index]=pos;
  setValue( 0,pos[0], 1.0 ); setValue(1, pos[1], 1.0 ) ; setValue( 2, pos[2], 1.0 ); 
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
  assert( d.size()==derivatives.size() );

  derivatives=d;
  for(unsigned i=0;i<getNumberOfAtoms();++i) {
     for(unsigned j=0;j<3;++j){
        addDerivative(0, 3*i+0, d[i](0,j) ); 
        addDerivative(1, 3*i+1, d[i](1,j) ); 
        addDerivative(2, 3*i+2, d[i](2,j) );
     }
  }
}

}

#endif
