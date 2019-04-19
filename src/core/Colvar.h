/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_core_Colvar_h
#define __PLUMED_core_Colvar_h

#include "ActionAtomistic.h"
#include "ActionWithValue.h"
#include <vector>

#define PLUMED_COLVAR_INIT(ao) Action(ao),Colvar(ao)

namespace PLMD {

/**
\ingroup INHERIT
This is the abstract base class to use for implementing new collective variables, within it there is
\ref AddingAColvar "information" as to how to go about implementing a new CV.
*/

class Colvar :
  public ActionAtomistic,
  public ActionWithValue
{
private:
protected:
  bool isEnergy;
  bool isExtraCV;
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
  void apply() override;
/// Set box derivatives automatically.
/// It should be called after the setAtomsDerivatives has been used for all
/// single atoms.
/// \warning It only works for collective variable NOT using PBCs!
  void           setBoxDerivativesNoPbc();
  void           setBoxDerivativesNoPbc(Value*);
public:
  bool checkIsEnergy() {return isEnergy;}
  explicit Colvar(const ActionOptions&);
  ~Colvar() {}
  static void registerKeywords( Keywords& keys );
  unsigned getNumberOfDerivatives() override;
};

inline
void Colvar::setAtomsDerivatives(Value*v,int i,const Vector&d) {
  v->addDerivative(3*i+0,d[0]);
  v->addDerivative(3*i+1,d[1]);
  v->addDerivative(3*i+2,d[2]);
}


inline
void Colvar::setBoxDerivatives(Value* v,const Tensor&d) {
  unsigned nat=getNumberOfAtoms();
  v->addDerivative(3*nat+0,d(0,0));
  v->addDerivative(3*nat+1,d(0,1));
  v->addDerivative(3*nat+2,d(0,2));
  v->addDerivative(3*nat+3,d(1,0));
  v->addDerivative(3*nat+4,d(1,1));
  v->addDerivative(3*nat+5,d(1,2));
  v->addDerivative(3*nat+6,d(2,0));
  v->addDerivative(3*nat+7,d(2,1));
  v->addDerivative(3*nat+8,d(2,2));
}

inline
void Colvar::setAtomsDerivatives(int i,const Vector&d) {
  setAtomsDerivatives(getPntrToValue(),i,d);
}

inline
void Colvar::setBoxDerivatives(const Tensor&d) {
  setBoxDerivatives(getPntrToValue(),d);
}

inline
void Colvar::setBoxDerivativesNoPbc() {
  setBoxDerivativesNoPbc(getPntrToValue());
}

inline
unsigned Colvar::getNumberOfDerivatives() {
  return 3*getNumberOfAtoms() + 9;
}


}

#endif
