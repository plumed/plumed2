#ifndef __PLUMED_ColvarWithoutModifiers_h
#define __PLUMED_ColvarWithoutModifiers_h

#include <string>
#include <cassert>
#include <vector>
#include "Colvar.h"

namespace PLMD {

class ColvarWithoutModifiers : public Colvar {
protected:
  void finishColvarSetup( const unsigned ncv, const double min, const double max );
public:
  ColvarWithoutModifiers(const ActionOptions&);
  void setAtomicDerivative( const unsigned& icv, const unsigned& iatom, const Vector& der ); 
  void setVirial( const unsigned& icv, const Tensor& vir );
};

inline
void ColvarWithoutModifiers::setAtomicDerivative( const unsigned& icv, const unsigned& iatom, const Vector& der ){
  addDerivative( icv, 3*iatom+0, der[0] );
  addDerivative( icv, 3*iatom+1, der[1] );
  addDerivative( icv, 3*iatom+2, der[2] );
}

inline
void ColvarWithoutModifiers::setVirial( const unsigned& icv, const Tensor& vir ){
  unsigned natoms=getNumberOfAtoms();
  addDerivative( icv, 3*natoms+0, vir(0,0) );
  addDerivative( icv, 3*natoms+1, vir(0,1) );
  addDerivative( icv, 3*natoms+2, vir(0,2) );
  addDerivative( icv, 3*natoms+3, vir(1,0) );
  addDerivative( icv, 3*natoms+4, vir(1,1) );
  addDerivative( icv, 3*natoms+5, vir(1,2) );
  addDerivative( icv, 3*natoms+6, vir(2,0) );
  addDerivative( icv, 3*natoms+7, vir(2,1) );
  addDerivative( icv, 3*natoms+8, vir(2,2) );
}

}

#endif
