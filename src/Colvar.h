#ifndef __PLUMED_Colvar_h
#define __PLUMED_Colvar_h

#include "ActionAtomistic.h"
#include "ActionWithValue.h"
#include <vector>

namespace PLMD {

/// Action representing a collective variable
class Colvar : public ActionAtomistic {
protected:
//  bool isEnergy;
  void requestAtoms(const std::vector<AtomNumber> & a);
// These are so as to map to 3d vectors for atoms:
  void setAtomsDerivatives(int,const Vector&);
  void setAtomsDerivatives(Value*,int,const Vector&);
  void setBoxDerivatives(const Tensor&);
  void setBoxDerivatives(Value*,const Tensor&);
  void setColvarValue(double);
  void setColvarValue(Value*,double);
public:
  const Tensor                & getBoxDerivatives()const;
  const double                & getForce()const;
  void apply();


public:


//  bool checkIsEnergy(){return isEnergy;};
  Colvar(const ActionOptions&);
  ~Colvar(){};

};


inline
void Colvar::setAtomsDerivatives(int i,const Vector&d){
  getValue(0)->setDerivatives(3*i+0,d[0]);
  getValue(0)->setDerivatives(3*i+1,d[1]);
  getValue(0)->setDerivatives(3*i+2,d[2]);
}

inline
void Colvar::setAtomsDerivatives(Value*v,int i,const Vector&d){
  v->setDerivatives(3*i+0,d[0]);
  v->setDerivatives(3*i+1,d[1]);
  v->setDerivatives(3*i+2,d[2]);
}


inline
void Colvar::setBoxDerivatives(const Tensor&d){
  getValue(0)->setDerivatives(3*getNatoms()+0,d(0,0));
  getValue(0)->setDerivatives(3*getNatoms()+1,d(0,1));
  getValue(0)->setDerivatives(3*getNatoms()+2,d(0,2));
  getValue(0)->setDerivatives(3*getNatoms()+3,d(1,0));
  getValue(0)->setDerivatives(3*getNatoms()+4,d(1,1));
  getValue(0)->setDerivatives(3*getNatoms()+5,d(1,2));
  getValue(0)->setDerivatives(3*getNatoms()+6,d(2,0));
  getValue(0)->setDerivatives(3*getNatoms()+7,d(2,1));
  getValue(0)->setDerivatives(3*getNatoms()+8,d(2,2));
}

inline
void Colvar::setBoxDerivatives(Value* v,const Tensor&d){
  v->setDerivatives(3*getNatoms()+0,d(0,0));
  v->setDerivatives(3*getNatoms()+1,d(0,1));
  v->setDerivatives(3*getNatoms()+2,d(0,2));
  v->setDerivatives(3*getNatoms()+3,d(1,0));
  v->setDerivatives(3*getNatoms()+4,d(1,1));
  v->setDerivatives(3*getNatoms()+5,d(1,2));
  v->setDerivatives(3*getNatoms()+6,d(2,0));
  v->setDerivatives(3*getNatoms()+7,d(2,1));
  v->setDerivatives(3*getNatoms()+8,d(2,2));
}

}

#endif
