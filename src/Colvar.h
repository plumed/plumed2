#ifndef __PLUMED_Colvar_h
#define __PLUMED_Colvar_h

#include "ActionAtomistic.h"
#include "ActionWithValue.h"
#include <vector>

#define PLUMED_COLVAR_INIT(ao) Action(ao),Colvar(ao)

namespace PLMD {

/// Action representing a collective variable
class Colvar :
  public ActionAtomistic,
  public ActionWithValue
  {
protected:
  bool isEnergy;
  void requestAtoms(const std::vector<AtomNumber> & a);
// Set the derivatives for a particular atom equal to the input Vector
// This routine is called setAtomsDerivatives because not because you
// are setting the derivative of many atoms but because you are setting
// the derivatives of a particular atom.  The s is an apostrophe s 
// but you can't put apostrophes in function names
  void           setAtomsDerivatives(int,const Vector&);
  void           setAtomsDerivatives(Value*,int,const Vector&);
  void           setBoxDerivatives(const Tensor&);
  void           setBoxDerivatives(Value*,const Tensor&);
  const Tensor & getBoxDerivatives()const;
  const double & getForce()const;
  void apply();
public:
  bool checkIsEnergy(){return isEnergy;};
  Colvar(const ActionOptions&);
  ~Colvar(){};
  static void registerKeywords( Keywords& keys );
};

inline
void Colvar::setAtomsDerivatives(Value*v,int i,const Vector&d){
  v->setDerivatives(3*i+0,d[0]);
  v->setDerivatives(3*i+1,d[1]);
  v->setDerivatives(3*i+2,d[2]);
}


inline
void Colvar::setBoxDerivatives(Value* v,const Tensor&d){
  unsigned nat=getNumberOfAtoms();
  v->setDerivatives(3*nat+0,d(0,0));
  v->setDerivatives(3*nat+1,d(0,1));
  v->setDerivatives(3*nat+2,d(0,2));
  v->setDerivatives(3*nat+3,d(1,0));
  v->setDerivatives(3*nat+4,d(1,1));
  v->setDerivatives(3*nat+5,d(1,2));
  v->setDerivatives(3*nat+6,d(2,0));
  v->setDerivatives(3*nat+7,d(2,1));
  v->setDerivatives(3*nat+8,d(2,2));
}

inline
void Colvar::setAtomsDerivatives(int i,const Vector&d){
  setAtomsDerivatives(getValue(0),i,d);
}

inline
void Colvar::setBoxDerivatives(const Tensor&d){
  setBoxDerivatives(getValue(0),d);
}

}

#endif
