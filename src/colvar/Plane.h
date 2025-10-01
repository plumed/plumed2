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
#ifndef __PLUMED_colvar_Plane_h
#define __PLUMED_colvar_Plane_h
#include "Colvar.h"
#include "ColvarInput.h"
#include "tools/Pbc.h"

#include <string>
#include <cmath>

namespace PLMD {
namespace colvar {

template <typename T>
class Plane : public Colvar {
private:
  bool pbc;
  std::vector<double> value;
  std::vector<double> derivs;
public:
  using precision=T;
  static void registerKeywords( Keywords& keys );
  explicit Plane(const ActionOptions&);
  static void parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa );
  static unsigned getModeAndSetupValues( ActionWithValue* av );
// active methods:
  void calculate() override;
  static void calculateCV( const ColvarInput<T>& cvin, ColvarOutput<T>& cvout );
};

template <typename T>
void Plane<T>::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.setDisplayName("PLANE");
  keys.add("atoms","ATOMS","the three or four atoms whose plane we are computing");
  keys.addOutputComponent("x","default","scalar/vector","the x-component of the vector that is normal to the plane containing the atoms");
  keys.addOutputComponent("y","default","scalar/vector","the y-component of the vector that is normal to the plane containing the atoms");
  keys.addOutputComponent("z","default","scalar/vector","the z-component of the vector that is normal to the plane containing the atoms");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.reset_style("NUMERICAL_DERIVATIVES","hidden");
}

template <typename T>
void Plane<T>::parseAtomList( const int& num, std::vector<AtomNumber>& atoms, ActionAtomistic* aa ) {
  aa->parseAtomList("ATOMS",num,atoms);
  if(atoms.size()==3) {
    aa->log.printf("  containing atoms %d %d %d\n",atoms[0].serial(),atoms[1].serial(),atoms[2].serial());
    atoms.resize(4);
    atoms[3]=atoms[2];
    atoms[2]=atoms[1];
  } else if(atoms.size()==4) {
    aa->log.printf("  containing lines %d-%d and %d-%d\n",atoms[0].serial(),atoms[1].serial(),atoms[2].serial(),atoms[3].serial());
  } else if( num<0 || atoms.size()>0 ) {
    aa->error("Number of specified atoms should be either 3 or 4");
  }
}

template <typename T>
unsigned Plane<T>::getModeAndSetupValues( ActionWithValue* av ) {
  av->addComponentWithDerivatives("x");
  av->componentIsNotPeriodic("x");
  av->addComponentWithDerivatives("y");
  av->componentIsNotPeriodic("y");
  av->addComponentWithDerivatives("z");
  av->componentIsNotPeriodic("z");
  return 0;
}

template <typename T>
Plane<T>::Plane(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  value(3) {
  std::vector<AtomNumber> atoms;
  parseAtomList(-1,atoms,this);
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  if(pbc) {
    log.printf("  using periodic boundary conditions\n");
  } else {
    log.printf("  without periodic boundary conditions\n");
  }

  /**/getModeAndSetupValues( this );
  requestAtoms(atoms);
  checkRead();
}

template <typename T>
void Plane<T>::calculate() {

  if(pbc) {
    makeWhole();
  }
  auto cvout = ColvarOutput<double>::createColvarOutput(value,derivs,this);
  Plane<double>::calculateCV( ColvarInput<double>::createColvarInput( 0, getPositions(), this ), cvout );
  Value* valuex=getPntrToComponent("x");
  Value* valuey=getPntrToComponent("y");
  Value* valuez=getPntrToComponent("z");

  for(unsigned i=0; i<getPositions().size(); ++i) {
    setAtomsDerivatives( valuex, i, cvout.getAtomDerivatives(0,i) );
  }
  setBoxDerivatives( valuex, cvout.virial[0] );
  valuex->set( value[0] );

  for(unsigned i=0; i<getPositions().size(); ++i) {
    setAtomsDerivatives( valuey, i, cvout.getAtomDerivatives(1,i) );
  }
  setBoxDerivatives( valuey, cvout.virial[1] );
  valuey->set( value[1] );

  for(unsigned i=0; i<getPositions().size(); ++i) {
    setAtomsDerivatives( valuez, i, cvout.getAtomDerivatives(2,i) );
  }
  setBoxDerivatives( valuez, cvout.virial[2] );
  valuez->set( value[2] );
}

template <typename T>
void Plane<T>::calculateCV( const ColvarInput<T>& cvin, ColvarOutput<T>& cvout ) {
  const auto d1=delta( cvin.pos[1], cvin.pos[0] );
  const auto d2=delta( cvin.pos[2], cvin.pos[3] );
  const auto cp = crossProduct( d1, d2 );
  using T33 = TensorTyped<T,3,3>;
  cvout.derivs[0][0] = crossProduct( Versors::xm<T>, d2 );
  cvout.derivs[0][1] = crossProduct( Versors::xp<T>, d2 );
  cvout.derivs[0][2] = crossProduct( Versors::xm<T>, d1 );
  cvout.derivs[0][3] = crossProduct( Versors::xp<T>, d1 );
  cvout.virial.set( 0, T33(d1,crossProduct(Versors::xp<T>, d2))
                    + T33(d2, crossProduct(Versors::xm<T>, d1)));
  cvout.values[0] = cp[0];

  cvout.derivs[1][0] = crossProduct( Versors::ym<T>, d2 );
  cvout.derivs[1][1] = crossProduct( Versors::yp<T>, d2 );
  cvout.derivs[1][2] = crossProduct( Versors::ym<T>, d1 );
  cvout.derivs[1][3] = crossProduct( Versors::yp<T>, d1 );
  cvout.virial.set(1, T33(d1,crossProduct(Versors::yp<T>, d2))
                   + T33(d2, crossProduct(Versors::ym<T>, d1)));
  cvout.values[1] = cp[1];

  cvout.derivs[2][0] = crossProduct( Versors::zm<T>, d2 );
  cvout.derivs[2][1] = crossProduct( Versors::zp<T>, d2 );
  cvout.derivs[2][2] = crossProduct( Versors::zm<T>, d1 );
  cvout.derivs[2][3] = crossProduct( Versors::zp<T>, d1 );
  cvout.virial.set(2, T33(d1,crossProduct(Versors::zp<T>, d2))
                   + T33(d2, crossProduct(Versors::zm<T>, d1)));
  cvout.values[2] = cp[2];
}

}
}

#endif // __PLUMED_colvar_Plane_h
