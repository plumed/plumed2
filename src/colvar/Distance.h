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
#ifndef __PLUMED_colvar_Distance_h
#define __PLUMED_colvar_Distance_h
#include "Colvar.h"
#include "ColvarInput.h"
#include "MultiColvarTemplate.h"
#include "tools/Pbc.h"

namespace PLMD {
namespace colvar {

template <typename T=double>
class Distance : public Colvar {
  enum components:unsigned {scaled=2,standard=1,scalar=0} mode{scalar};
  bool pbc;

  //not Typep, the float version is only used for the PTM
  std::vector<double> value;
  std::vector<double> derivs;
public:
  using precision = T;
  static void registerKeywords( Keywords& keys );
  explicit Distance(const ActionOptions&);
  static void parseAtomList( const int& num,
                             std::vector<AtomNumber>& t,
                             ActionAtomistic* aa );
  static unsigned getModeAndSetupValues( ActionWithValue* av );
// active methods:
  void calculate() override;
  static void calculateCV( const ColvarInput<T>& cvin,
                           ColvarOutput<T>& cvout );
};

template <typename T>
void  Distance<T>::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.setDisplayName("DISTANCE");
  constexpr auto scalarOrVector = Keywords::componentType::scalar
                                  | Keywords::componentType::vector;
  keys.add("atoms","ATOMS","the pair of atom that we are calculating the distance between");
  keys.addFlag("COMPONENTS",false,"calculate the x, y and z components of the distance separately and store them as label.x, label.y and label.z");
  keys.addFlag("SCALED_COMPONENTS",false,"calculate the a, b and c scaled components of the distance separately and store them as label.a, label.b and label.c");
  keys.addOutputComponent("x","COMPONENTS",scalarOrVector,"the x-component of the vector connecting the two atoms");
  keys.addOutputComponent("y","COMPONENTS",scalarOrVector,"the y-component of the vector connecting the two atoms");
  keys.addOutputComponent("z","COMPONENTS",scalarOrVector,"the z-component of the vector connecting the two atoms");
  keys.addOutputComponent("a","SCALED_COMPONENTS",scalarOrVector,"the normalized projection on the first lattice vector of the vector connecting the two atoms");
  keys.addOutputComponent("b","SCALED_COMPONENTS",scalarOrVector,"the normalized projection on the second lattice vector of the vector connecting the two atoms");
  keys.addOutputComponent("c","SCALED_COMPONENTS",scalarOrVector,"the normalized projection on the third lattice vector of the vector connecting the two atoms");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.setValueDescription(scalarOrVector,"the DISTANCE between this pair of atoms");
  keys.addDOI("10.1007/978-1-4939-9608-7_21");
  keys.reset_style("NUMERICAL_DERIVATIVES","hidden");
}

template <typename T>
Distance<T>::Distance(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  value(1) {
  std::vector<AtomNumber> atoms;
  parseAtomList(-1,atoms,this);
  if(atoms.size()!=2) {
    error("Number of specified atoms should be 2");
  }

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  if(pbc) {
    log.printf("  using periodic boundary conditions\n");
  } else {
    log.printf("  without periodic boundary conditions\n");
  }

  mode = static_cast<components>(getModeAndSetupValues(this));
  if( mode == components::standard || mode == components::scaled ) {
    value.resize(3);
  }
  requestAtoms(atoms);
}

template <typename T>
void  Distance<T>::parseAtomList( const int& num,
                                  std::vector<AtomNumber>& t,
                                  ActionAtomistic* aa ) {
  aa->parseAtomList("ATOMS",num,t);
  if( t.size()==2 ) {
    aa->log.printf("  between atoms %d %d\n",t[0].serial(),t[1].serial());
  }
}

template <typename T>
unsigned  Distance<T>::getModeAndSetupValues( ActionWithValue* av ) {
  bool c;
  av->parseFlag("COMPONENTS",c);
  bool sc;
  av->parseFlag("SCALED_COMPONENTS",sc);
  if( c && sc ) {
    av->error("COMPONENTS and SCALED_COMPONENTS are not compatible");
  }
  if(c) {
    av->addComponentWithDerivatives("x");
    av->componentIsNotPeriodic("x");
    av->addComponentWithDerivatives("y");
    av->componentIsNotPeriodic("y");
    av->addComponentWithDerivatives("z");
    av->componentIsNotPeriodic("z");
    av->log<<"  WARNING: components will not have the proper periodicity - see manual\n";
    return components::standard;
  } else if(sc) {
    av->addComponentWithDerivatives("a");
    av->componentIsPeriodic("a","-0.5","+0.5");
    av->addComponentWithDerivatives("b");
    av->componentIsPeriodic("b","-0.5","+0.5");
    av->addComponentWithDerivatives("c");
    av->componentIsPeriodic("c","-0.5","+0.5");
    return components::scaled;
  }
  av->addValueWithDerivatives();
  av->setNotPeriodic();
  return components::scalar;
}

// calculator
template <typename T>
void  Distance<T>::calculate() {

  if(pbc) {
    makeWhole();
  }
  using CVInput= ColvarInput<double>;
  using CVOutput=ColvarOutput<double>;
  auto cvout =  CVOutput::createColvarOutput(value, derivs, this);
  Distance<double>::calculateCV( CVInput::createColvarInput( mode, getPositions(), this ), cvout );
  if( mode==components::standard ) {
    Value* valuex=getPntrToComponent("x");
    Value* valuey=getPntrToComponent("y");
    Value* valuez=getPntrToComponent("z");

    for(unsigned i=0; i<2; ++i) {
      setAtomsDerivatives(valuex,i,cvout.getAtomDerivatives(0,i) );
    }
    setBoxDerivatives(valuex,cvout.virial[0]);
    valuex->set(value[0]);

    for(unsigned i=0; i<2; ++i) {
      setAtomsDerivatives(valuey,i,cvout.getAtomDerivatives(1,i) );
    }
    setBoxDerivatives(valuey,cvout.virial[1]);
    valuey->set(value[1]);

    for(unsigned i=0; i<2; ++i) {
      setAtomsDerivatives(valuez,i,cvout.getAtomDerivatives(2,i) );
    }
    setBoxDerivatives(valuez,cvout.virial[2]);
    valuez->set(value[2]);
  } else if( mode==components::scaled ) {

    Value* valuea=getPntrToComponent("a");
    Value* valueb=getPntrToComponent("b");
    Value* valuec=getPntrToComponent("c");
    for(unsigned i=0; i<2; ++i) {
      setAtomsDerivatives(valuea,i,cvout.getAtomDerivatives(0,i) );
    }
    valuea->set(value[0]);
    for(unsigned i=0; i<2; ++i) {
      setAtomsDerivatives(valueb,i,cvout.getAtomDerivatives(1,i) );
    }
    valueb->set(value[1]);
    for(unsigned i=0; i<2; ++i) {
      setAtomsDerivatives(valuec,i,cvout.getAtomDerivatives(2,i) );
    }
    valuec->set(value[2]);
  } else  {
    for(unsigned i=0; i<2; ++i) {
      setAtomsDerivatives(i,cvout.getAtomDerivatives(0,i) );
    }
    setBoxDerivatives(cvout.virial[0]);
    setValue           (value[0]);
  }
}

template <typename T>
void  Distance<T>::calculateCV( const ColvarInput<T>& cvin, ColvarOutput<T>& cvout ) {
  auto distance=delta(cvin.pos[0],cvin.pos[1]);
  //these six can be compacted in a scoped `using namespac3 PLMD::Versors` but plumedcheck won't approve (even in a comment, hence the 3)
  using PLMD::Versors::xp;
  using PLMD::Versors::xm;
  using PLMD::Versors::yp;
  using PLMD::Versors::ym;
  using PLMD::Versors::zp;
  using PLMD::Versors::zm;
  if(cvin.mode==components::standard) {

    cvout.derivs[0][0] = xm<T>;
    cvout.derivs[0][1] = xp<T>;
    cvout.values[0] = distance[0];

    cvout.derivs[1][0] = ym<T>;
    cvout.derivs[1][1] = yp<T>;
    cvout.values[1] = distance[1];

    cvout.derivs[2][0] = zm<T>;
    cvout.derivs[2][1] = zp<T>;
    cvout.values[2] = distance[2];
    ColvarInput<T>::setBoxDerivativesNoPbc( cvin, cvout );
  } else if(cvin.mode==components::scaled) {
    auto d=cvin.pbc.realToScaled(distance);
    cvout.derivs[0][0].copyConv(matmul(cvin.pbc.getInvBox(),xm<T>));
    cvout.derivs[0][1].copyConv(matmul(cvin.pbc.getInvBox(),xp<T>));
    cvout.values[0] = Tools::pbc(d[0]);
    cvout.derivs[1][0].copyConv(matmul(cvin.pbc.getInvBox(),ym<T>));
    cvout.derivs[1][1].copyConv(matmul(cvin.pbc.getInvBox(),yp<T>));
    cvout.values[1] = Tools::pbc(d[1]);
    cvout.derivs[2][0].copyConv(matmul(cvin.pbc.getInvBox(),zm<T>));
    cvout.derivs[2][1].copyConv(matmul(cvin.pbc.getInvBox(),zp<T>));
    cvout.values[2] = Tools::pbc(d[2]);
  } else {
    const double value=distance.modulo();
    const double invvalue=1.0/value;
    cvout.derivs[0][0] = -invvalue*distance;
    cvout.derivs[0][1] = invvalue*distance;
    ColvarInput<T>::setBoxDerivativesNoPbc( cvin, cvout );
    cvout.values[0] = value;
  }
}

} //namespace colvar
} //namespace PLMD

#endif //__PLUMED_colvar_Distance_h


