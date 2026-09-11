/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2023 The plumed team
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
#ifndef __PLUMED_colvar_Position_h
#define __PLUMED_colvar_Position_h
#include "Colvar.h"
#include "tools/Pbc.h"
#include "ColvarInput.h"

namespace PLMD {
namespace colvar {

template <typename T>
class Position : public Colvar {
  enum components:unsigned {scaled=1,standard=0} mode{standard};
  bool pbc;
  std::vector<double> value;
  std::vector<double> derivs;
public:
  using precision=T;
  static void registerKeywords( Keywords& keys );
  explicit Position(const ActionOptions&);
  static void parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa );
  static unsigned getModeAndSetupValues( ActionWithValue* av );
// active methods:
  void calculate() override;
  static void calculateCV( const ColvarInput<T>& cvin, ColvarOutput<T>& cvout );
};

template <typename T>
void Position<T>::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.setDisplayName("POSITION");
  keys.add("atoms","ATOM","the atom number");
  keys.add("atoms","ATOMS","the atom numbers that you would like to use the positions of");
  keys.addFlag("WHOLEMOLECULES",false,"if this is a vector of positions do you want to make the positions into a whole before");
  keys.addFlag("SCALED_COMPONENTS",false,"calculate the a, b and c scaled components of the position separately and store them as label.a, label.b and label.c");
  keys.addOutputComponent("x","default","scalar/vector","the x-component of the atom position");
  keys.addOutputComponent("y","default","scalar/vector","the y-component of the atom position");
  keys.addOutputComponent("z","default","scalar/vector","the z-component of the atom position");
  keys.addOutputComponent("a","SCALED_COMPONENTS","scalar/vector","the normalized projection on the first lattice vector of the atom position");
  keys.addOutputComponent("b","SCALED_COMPONENTS","scalar/vector","the normalized projection on the second lattice vector of the atom position");
  keys.addOutputComponent("c","SCALED_COMPONENTS","scalar/vector","the normalized projection on the third lattice vector of the atom position");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.reset_style("NUMERICAL_DERIVATIVES","hidden");
}

template <typename T>
Position<T>::Position(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  value(3) {
  std::vector<AtomNumber> atoms;
  parseAtomList(-1,atoms,this);
  mode=static_cast<components>(getModeAndSetupValues(this));

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  checkRead();

  if(pbc) {
    log.printf("  using periodic boundary conditions\n");
  } else {
    log.printf("  without periodic boundary conditions\n");
  }

  requestAtoms(atoms);
}

template <typename T>
void Position<T>::parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa ) {
  aa->parseAtomList("ATOM",num,t);
  if( t.size()==1 ) {
    aa->log.printf("  for atom %d\n",t[0].serial());
  } else if( num<0 || t.size()!=0 ) {
    aa->error("Number of specified atoms should be 1");
  }
}

template <typename T>
unsigned Position<T>::getModeAndSetupValues( ActionWithValue* av ) {
  bool sc;
  av->parseFlag("SCALED_COMPONENTS",sc);
  if(sc) {
    av->addComponentWithDerivatives("a");
    av->componentIsPeriodic("a","-0.5","+0.5");
    av->addComponentWithDerivatives("b");
    av->componentIsPeriodic("b","-0.5","+0.5");
    av->addComponentWithDerivatives("c");
    av->componentIsPeriodic("c","-0.5","+0.5");
    return components::scaled;
  }
  av->addComponentWithDerivatives("x");
  av->componentIsNotPeriodic("x");
  av->addComponentWithDerivatives("y");
  av->componentIsNotPeriodic("y");
  av->addComponentWithDerivatives("z");
  av->componentIsNotPeriodic("z");
  av->log<<"  WARNING: components will not have the proper periodicity - see manual\n";
  return components::standard;
}

// calculator
template <typename T>
void Position<T>::calculate() {

  std::vector<Vector> distance(1);
  if(pbc) {
    distance[0]=pbcDistance(Vector(0.0,0.0,0.0),getPosition(0));
  } else {
    distance[0]=delta(Vector(0.0,0.0,0.0),getPosition(0));
  }

  ColvarOutput<double> cvout = ColvarOutput<double>::createColvarOutput(value,derivs,this);
  Position<double>::calculateCV( ColvarInput<double>::createColvarInput( mode, distance, this ), cvout );
  if(mode==components::scaled) {
    Value* valuea=getPntrToComponent("a");
    Value* valueb=getPntrToComponent("b");
    Value* valuec=getPntrToComponent("c");
    setAtomsDerivatives (valuea,0,cvout.getAtomDerivatives(0,0));
    valuea->set(value[0]);
    setAtomsDerivatives (valueb,0,cvout.getAtomDerivatives(1,0));
    valueb->set(value[1]);
    setAtomsDerivatives (valuec,0,cvout.getAtomDerivatives(2,0));
    valuec->set(value[2]);
  } else {
    Value* valuex=getPntrToComponent("x");
    Value* valuey=getPntrToComponent("y");
    Value* valuez=getPntrToComponent("z");

    setAtomsDerivatives (valuex,0,cvout.getAtomDerivatives(0,0));
    setBoxDerivatives   (valuex,cvout.virial[0]);
    valuex->set(value[0]);

    setAtomsDerivatives (valuey,0,cvout.getAtomDerivatives(1,0));
    setBoxDerivatives   (valuey,cvout.virial[1]);
    valuey->set(value[1]);

    setAtomsDerivatives (valuez,0,cvout.getAtomDerivatives(2,0));
    setBoxDerivatives   (valuez,cvout.virial[2]);
    valuez->set(value[2]);
  }
}

template <typename T>
void Position<T>::calculateCV( const ColvarInput<T>& cvin, ColvarOutput<T>& cvout ) {
  if( cvin.mode==components::scaled ) {
    auto d=cvin.pbc.realToScaled(VectorTyped<T,3>(cvin.pos[0][0],cvin.pos[0][1],cvin.pos[0][2]));
    cvout.values[0]=Tools::pbc(d[0]);
    cvout.values[1]=Tools::pbc(d[1]);
    cvout.values[2]=Tools::pbc(d[2]);
    cvout.derivs[0][0]=matmul(cvin.pbc.getInvBox(),Versors::xp<T>);
    cvout.derivs[1][0]=matmul(cvin.pbc.getInvBox(),Versors::yp<T>);
    cvout.derivs[2][0]=matmul(cvin.pbc.getInvBox(),Versors::zp<T>);
  } else {
    for(unsigned i=0; i<3; ++i) {
      cvout.values[i]=cvin.pos[0][i];
    }
    cvout.derivs[0][0]=Versors::xp<T>;
    cvout.derivs[1][0]=Versors::yp<T>;
    cvout.derivs[2][0]=Versors::zp<T>;
    cvout.virial.set(0, TensorTyped<T,3,3>(VectorTyped<T,3>(cvin.pos[0][0],cvin.pos[0][1],cvin.pos[0][2]),
                                           Versors::xm<T>));
    cvout.virial.set(1, TensorTyped<T,3,3>(VectorTyped<T,3>(cvin.pos[0][0],cvin.pos[0][1],cvin.pos[0][2]),
                                           Versors::ym<T>));
    cvout.virial.set(2, TensorTyped<T,3,3>(VectorTyped<T,3>(cvin.pos[0][0],cvin.pos[0][1],cvin.pos[0][2]),
                                           Versors::zm<T>));
  }
}

}
}
#endif //__PLUMED_colvar_Position_h
