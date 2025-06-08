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
#ifdef __PLUMED_HAS_OPENACC
#define __PLUMED_USE_OPENACC 1
#endif //__PLUMED_HAS_OPENACC
#include "Colvar.h"
#include "ColvarShortcut.h"
#include "MultiColvarTemplate.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR DIPOLE
/*
Calculate the dipole moment for a group of atoms.

The following input tells plumed to calculate the dipole for the group of atoms containing
the atoms from 1-10 and print it every 5 steps

```plumed
d: DIPOLE GROUP=1-10
PRINT FILE=output STRIDE=5 ARG=d
```

The output value from this input is a scalar that tells you the total magnitude of the dipole vector.  If you would
like to access the dipole vector directly you can use the command:

```plumed
d: DIPOLE GROUP=1-10 COMPONENTS
PRINT ARG=d.* FILE=output
```

This command will output three values d.x, d.y and d.z, which are the x, y and z components of the dipole respectively.


You can calculate three instinguishable dipoles using a single DIPOLE command by using an input like the one below:

```plumed
d: DIPOLE GROUP1=1-10 GROUP2=11-20 GROUP3=21-30
PRINT ARG=d FILE=output
```

The output, d, here is a three dimensional vector.  The first element of this vector is the magnitude of the dipole for
atoms 1-10, the second is the magnitude of the dipole for atoms 11-20 and the third is the magnitude of the dipole for
atoms 21-30.  You can also obtain vector components for the three dipoles above by using the following input:

```plumed
d: DIPOLE COMPONENTS GROUP1=1-10 GROUP2=11-20 GROUP3=21-30
PRINT ARG=d.x,d.y,d.z FILE=output
```

The output from the DIPOLE command now consists of three three dimensional vectors called d.x, d.y and d.z that contain the
x, y and z components of the three dipoles respectively.

A final important thing to note is that in all the commands above the default is to use a procedure akin to that used in [WHOLEMOLECULES](WHOLEMOLECULES.md) to ensure
that the sets of atoms that are specified to each GROUP keyword are not broken by the periodic
boundary conditions.  If you would like to turn this off for any reason you add the NOPBC in your input file as shown
below:

```plumed
d: DIPOLE GROUP=1-10 NOPBC
PRINT FILE=output STRIDE=5 ARG=d
```

!!! caution "behaviour for non charge neutral groups"

    If the total charge Q of any of the specified groups is non zero, then a charge Q/N will be subtracted from every atom,
    where N is the number of atoms. This implies that the dipole (which for a charged system depends
    on the position) is computed on the geometric center of the group.

*/
//+ENDPLUMEDOC

class Dipole : public Colvar {
  std::vector<AtomNumber> ga_lista;
  bool components;
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
  static void calculateCV( const ColvarInput& cvin, ColvarOutput& cvout );
};

typedef ColvarShortcut<Dipole> DipoleShortcut;
PLUMED_REGISTER_ACTION(DipoleShortcut,"DIPOLE")
PLUMED_REGISTER_ACTION(Dipole,"DIPOLE_SCALAR")
typedef MultiColvarTemplate<Dipole> DipoleMulti;
PLUMED_REGISTER_ACTION(DipoleMulti,"DIPOLE_VECTOR")

void Dipole::registerKeywords(Keywords& keys) {
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

Dipole::Dipole(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  components(false),
  value(1) {
  parseAtomList(-1,ga_lista,this);
  components=(getModeAndSetupValues(this)==1);
  if( components ) {
    value.resize(3);
    valuex=getPntrToComponent("x");
    valuey=getPntrToComponent("y");
    valuez=getPntrToComponent("z");
  } else {
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

void Dipole::parseAtomList( const int& num, std::vector<AtomNumber>& t, ActionAtomistic* aa ) {
  aa->parseAtomList("GROUP",num,t);
  if( t.size()>0 ) {
    aa->log.printf("  of %u atoms\n",static_cast<unsigned>(t.size()));
    for(unsigned int i=0; i<t.size(); ++i) {
      aa->log.printf("  %d", t[i].serial());
    }
    aa->log.printf("  \n");
  }
}

unsigned Dipole::getModeAndSetupValues( ActionWithValue* av ) {
  bool c;
  av->parseFlag("COMPONENTS",c);
  if( c ) {
    av->addComponentWithDerivatives("x");
    av->componentIsNotPeriodic("x");
    av->addComponentWithDerivatives("y");
    av->componentIsNotPeriodic("y");
    av->addComponentWithDerivatives("z");
    av->componentIsNotPeriodic("z");
    return 1;
  }
  av->addValueWithDerivatives();
  av->setNotPeriodic();
  return 0;
}

// calculator
void Dipole::calculate() {
  if( !chargesWereSet ) {
    error("charges were not set by MD code");
  }

  if(!nopbc) {
    makeWhole();
  }
  unsigned N=getNumberOfAtoms();

  if(!components) {
    ColvarOutput cvout( ColvarOutput::createColvarOutput(value, derivs, this) );
    calculateCV( ColvarInput::createColvarInput( 0, getPositions(), this ), cvout );
    for(unsigned i=0; i<N; i++) {
      setAtomsDerivatives(i,cvout.getAtomDerivatives(0,i));
    }
    setBoxDerivatives(cvout.virial[0]);
    setValue(value[0]);
  } else {
    ColvarOutput cvout( ColvarOutput::createColvarOutput(value, derivs, this) );
    calculateCV( ColvarInput::createColvarInput( 1, getPositions(), this ), cvout );
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

void Dipole::calculateCV( const ColvarInput& cvin, ColvarOutput& cvout ) {
  unsigned N=cvin.pos.size();
  double ctot=0.;
  for(unsigned i=0; i<N; ++i) {
    ctot += cvin.charges[i];
  }
  ctot/=(double)N;

  Vector dipje;
  for(unsigned i=0; i<N; ++i) {
    dipje += (cvin.charges[i]-ctot)*Vector(cvin.pos[i][0],cvin.pos[i][1],cvin.pos[i][2]);
  }

  if( cvin.mode==1 ) {
    for(unsigned i=0; i<N; i++) {
      cvout.derivs[0][i]=(cvin.charges[i]-ctot)*Vector(1.0,0.0,0.0);
      cvout.derivs[1][i]=(cvin.charges[i]-ctot)*Vector(0.0,1.0,0.0);
      cvout.derivs[2][i]=(cvin.charges[i]-ctot)*Vector(0.0,0.0,1.0);
    }
    for(unsigned i=0; i<3; ++i ) {
      cvout.values[i] = dipje[i];
    }
  } else {
    cvout.values[0] = dipje.modulo();
    double idip = 1./cvout.values[0];
    for(unsigned i=0; i<N; i++) {
      double dfunc=(cvin.charges[i]-ctot)*idip;
      cvout.derivs[0][i] = dfunc*dipje;
    }
  }
  ColvarInput::setBoxDerivativesNoPbc( cvin, cvout );
}

}
}
