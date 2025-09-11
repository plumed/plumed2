
/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#ifndef __PLUMED_colvar_Dipole_h
#define __PLUMED_colvar_Dipole_h
#include "Colvar.h"
#include "ColvarInput.h"

namespace PLMD {
namespace colvar {

template <typename T>
class Dipole : public Colvar {
  std::vector<AtomNumber> ga_lista;
  enum /*class*/ modes:unsigned {
    components=1,
    scalar=0
  } mode{scalar};
  bool nopbc;
  std::vector<double> value;
  std::vector<double> derivs;
  Value* valuex=nullptr;
  Value* valuey=nullptr;
  Value* valuez=nullptr;
public:
  explicit Dipole(const ActionOptions&);
  static void parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa );
  static unsigned getModeAndSetupValues( ActionWithValue* av );
  void calculate() override;
  static void registerKeywords(Keywords& keys);
  static void calculateCV( const ColvarInput<T>& cvin, ColvarOutput<T>& cvout );
};

template <typename T>
void Dipole<T>::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.setDisplayName("DIPOLE");
  keys.add("atoms","GROUP","the group of atoms we are calculating the dipole moment for");
  keys.addFlag("COMPONENTS",false,"calculate the x, y and z components of the dipole separately and store them as label.x, label.y and label.z");
  keys.addOutputComponent("x","COMPONENTS","scalar/vector","the x-component of the dipole");
  keys.addOutputComponent("y","COMPONENTS","scalar/vector","the y-component of the dipole");
  keys.addOutputComponent("z","COMPONENTS","scalar/vector","the z-component of the dipole");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.setValueDescription("scalar/vector","the DIPOLE for these atoms");
  keys.reset_style("NUMERICAL_DERIVATIVES","hidden");
}

template <typename T>
Dipole<T>::Dipole(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao) {
  parseAtomList(-1,ga_lista,this);
  mode=static_cast<modes>(getModeAndSetupValues(this));
  if( mode==modes::components ) {
    value.resize(3);
    valuex=getPntrToComponent("x");
    valuey=getPntrToComponent("y");
    valuez=getPntrToComponent("z");
  } else {
    value.resize(1);
    derivs.resize(1,ga_lista.size());
  }
  parseFlag("NOPBC",nopbc);
  checkRead();

  if(nopbc) {
    log.printf("  without periodic boundary conditions\n");
  } else {
    log.printf("  using periodic boundary conditions\n");
  }

  requestAtoms(ga_lista);
}

template <typename T>
void Dipole<T>::parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa ) {
  aa->parseAtomList("GROUP",num,t);
  if( t.size()>0 ) {
    aa->log.printf("  of %u atoms\n",static_cast<unsigned>(t.size()));
    for(unsigned int i=0; i<t.size(); ++i) {
      aa->log.printf("  %d", t[i].serial());
    }
    aa->log.printf("  \n");
  }
}

template <typename T>
unsigned Dipole<T>::getModeAndSetupValues( ActionWithValue* av ) {
  bool c;
  av->parseFlag("COMPONENTS",c);
  if( c ) {
    av->addComponentWithDerivatives("x");
    av->componentIsNotPeriodic("x");
    av->addComponentWithDerivatives("y");
    av->componentIsNotPeriodic("y");
    av->addComponentWithDerivatives("z");
    av->componentIsNotPeriodic("z");
    return modes::components;
  }
  av->addValueWithDerivatives();
  av->setNotPeriodic();
  return modes::scalar;
}

// calculator
template <typename T>
void Dipole<T>::calculate() {
  if( !chargesWereSet ) {
    error("charges were not set by MD code");
  }

  if(!nopbc) {
    makeWhole();
  }
  const unsigned N=getNumberOfAtoms();

  auto cvout= ColvarOutput<double>::createColvarOutput(value, derivs, this);
  Dipole<double>::calculateCV( ColvarInput<double>::createColvarInput( mode, getPositions(), this ),
                               cvout );
  if(mode==modes::scalar) {
    for(unsigned i=0; i<N; i++) {
      setAtomsDerivatives(i,cvout.getAtomDerivatives(0,i));
    }
    setBoxDerivatives(cvout.virial[0]);
    setValue(value[0]);
  } else {
    for(unsigned i=0; i<N; i++) {
      setAtomsDerivatives(valuex,i,cvout.getAtomDerivatives(0,i));
      setAtomsDerivatives(valuey,i,cvout.getAtomDerivatives(1,i));
      setAtomsDerivatives(valuez,i,cvout.getAtomDerivatives(2,i));
    }
    setBoxDerivatives(valuex,cvout.virial[0]);
    setBoxDerivatives(valuey,cvout.virial[1]);
    setBoxDerivatives(valuez,cvout.virial[2]);
    valuex->set(value[0]);
    valuey->set(value[1]);
    valuez->set(value[2]);
  }
}

template <typename T>
void Dipole<T>::calculateCV( const ColvarInput<T>& cvin, ColvarOutput<T>& cvout ) {
  using V3 = VectorTyped<T,3>;
  const unsigned N=cvin.pos.size();
  T ctot=0.;
  for(unsigned i=0; i<N; ++i) {
    ctot += cvin.charges[i];
  }
  ctot/=T(N);

  V3 dipje;
  for(unsigned i=0; i<N; ++i) {
    dipje += (cvin.charges[i]-ctot)*V3(cvin.pos[i][0],cvin.pos[i][1],cvin.pos[i][2]);
  }

  if( cvin.mode==modes::components ) {
    for(unsigned i=0; i<N; i++) {
      cvout.derivs[0][i]=(cvin.charges[i]-ctot)*V3(1.0,0.0,0.0);
      cvout.derivs[1][i]=(cvin.charges[i]-ctot)*V3(0.0,1.0,0.0);
      cvout.derivs[2][i]=(cvin.charges[i]-ctot)*V3(0.0,0.0,1.0);
    }
    for(unsigned i=0; i<3; ++i ) {
      cvout.values[i] = dipje[i];
    }
  } else {
    cvout.values[0] = dipje.modulo();
    T idip = 1./cvout.values[0];
    for(unsigned i=0; i<N; i++) {
      T dfunc=(cvin.charges[i]-ctot)*idip;
      cvout.derivs[0][i] = dfunc*dipje;
    }
  }
  ColvarInput<T>::setBoxDerivativesNoPbc( cvin, cvout );
}

}
}
#endif //__PLUMED_colvar_Dipole_h
