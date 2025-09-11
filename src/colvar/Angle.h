/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#ifndef __PLUMED_colvar_Angle_h
#define __PLUMED_colvar_Angle_h
#include "Colvar.h"
#include "ColvarInput.h"
#include "tools/Angle.h"

namespace PLMD {
namespace colvar {

template <typename T>
class Angle : public Colvar {
  bool pbc;
  std::vector<double> value;
  std::vector<double> derivs;
public:
  explicit Angle(const ActionOptions&);
// active methods:
  void calculate() override;
  static void registerKeywords( Keywords& keys );
  static void parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa );
  static unsigned getModeAndSetupValues( ActionWithValue* av );
  static void calculateCV( const ColvarInput<T>& cvin, ColvarOutput<T>& cvout );
};

template <typename T>
void Angle<T>::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords(keys);
  keys.setDisplayName("ANGLE");
  keys.add("atoms","ATOMS","the list of atoms involved in this collective variable (either 3 or 4 atoms)");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.setValueDescription("scalar/vector","the ANGLE involving these atoms");
  keys.reset_style("NUMERICAL_DERIVATIVES","hidden");
}

template <typename T>
void Angle<T>::parseAtomList( const int& num, std::vector<AtomNumber>& atoms, ActionAtomistic* aa ) {
  aa->parseAtomList("ATOMS",num,atoms);
  if(atoms.size()==3) {
    aa->log.printf("  between atoms %d %d %d\n",atoms[0].serial(),atoms[1].serial(),atoms[2].serial());
    atoms.resize(4);
    atoms[3]=atoms[2];
    atoms[2]=atoms[1];
  } else if(atoms.size()==4) {
    aa->log.printf("  between lines %d-%d and %d-%d\n",atoms[0].serial(),atoms[1].serial(),atoms[2].serial(),atoms[3].serial());
  } else if( num<0 || atoms.size()>0 ) {
    aa->error("Number of specified atoms should be either 3 or 4");
  }
}

template <typename T>
unsigned Angle<T>::getModeAndSetupValues( ActionWithValue* av ) {
  av->addValueWithDerivatives();
  av->setNotPeriodic();
  return 0;
}

template <typename T>
Angle<T>::Angle(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  value(1) {
  std::vector<AtomNumber> atoms;
  parseAtomList( -1, atoms, this );
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  if(pbc) {
    log.printf("  using periodic boundary conditions\n");
  } else {
    log.printf("  without periodic boundary conditions\n");
  }

  addValueWithDerivatives();
  setNotPeriodic();
  requestAtoms(atoms);
  checkRead();
}

// calculator
template <typename T>
void Angle<T>::calculate() {

  if(pbc) {
    makeWhole();
  }
  ColvarOutput<double> cvout = ColvarOutput<double>::createColvarOutput(value,derivs,this);
  Angle<double>::calculateCV( ColvarInput<double>::createColvarInput( 0, getPositions(), this ), cvout );
  setValue( value[0] );
  for(unsigned i=0; i<getPositions().size(); ++i) {
    setAtomsDerivatives( i, cvout.getAtomDerivatives(0,i) );
  }
  setBoxDerivatives( cvout.virial[0] );
}

template <typename T>
void Angle<T>::calculateCV( const ColvarInput<T>& cvin, ColvarOutput<T>& cvout ) {
  auto dij=delta(cvin.pos[2],cvin.pos[3]);
  auto dik=delta(cvin.pos[1],cvin.pos[0]);
  VectorTyped<T,3> ddij,ddik;

  cvout.values[0]=PLMD::Angle::compute(dij,dik,ddij,ddik);
  cvout.derivs[0][0]=ddik;
  cvout.derivs[0][1]=-ddik;
  cvout.derivs[0][2]=-ddij;
  cvout.derivs[0][3]=ddij;
  ColvarInput<T>::setBoxDerivativesNoPbc( cvin, cvout );
}

} //namespace colvar
} //namespace PLMD

#endif //__PLUMED_colvar_Angle_h
