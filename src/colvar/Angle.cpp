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
#ifdef __PLUMED_HAS_OPENACC
#define __PLUMED_USE_OPENACC 1
#endif //__PLUMED_HAS_OPENACC
#include "Colvar.h"
#include "ColvarShortcut.h"
#include "MultiColvarTemplate.h"
#include "core/ActionRegister.h"
#include "tools/Angle.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR ANGLE
/*
Calculate one or multiple angle/s.

The following input instructs PLUMED to calculate and print the angle between the vector
connecting atom 2 and atom 1 and the vector connecting atom 2 and atom 3.

```plumed
a1: ANGLE ATOMS=1,2,3
PRINT ARG=a1 FILE=colvar
```

In other words, the angle that is output by the input above is calculated as:

$$
\theta=\arccos\left(\frac{ r_{21}\cdot r_{23}}{
| r_{21}| | r_{23}|}\right)
$$

Here $r_{ij}$ is the vector connecting the $i$th and $j$th atoms, which by default is evaluated
in a way that takes periodic boundary conditions into account. If you wish to disregard the PBC you
can use the NOPBC flag as shown in the following input:

```plumed
a1: ANGLE ATOMS=1,2,3 NOPBC
PRINT ARG=a1 FILE=colvar
```

If the NOPBC flag is not included any sets of atoms that are broken by the periodic boundaries are made whole
using a procedure that is the same as that described in the documentation for [WHOLEMOLECULES](WHOLEMOLECULES.md).

We can also instruct PLUMED to calculate the angle between the vectors connecting
atoms 1 and atom 2 and atoms 3 and atom 4 by using the following input:

```plumed
a2: ANGLE ATOMS=1,2,3,4
PRINT ARG=a2 FILE=colvar

The angle in this input is calculated using:

$$
\theta=\arccos\left(\frac{ r_{21}\cdot r_{34}}{
| r_{21}| | r_{34}|}\right)
$$

Notice that angles defined in this way are non-periodic variables - their values must lie in between 0 and $\pi$.

You can specify multiple sets of three or four atoms to calculate vectors of angles as illustrated in the following
input which instructs PLUMED to calculate and output three angles:

```plumed
a3: ANGLE ATOMS1=1,2,3 ATOMS2=4,5,6 ATOMS3=7,8,9
PRINT ARG=a3 FILE=colvar
```

It is common to assume when using this feature that all the angles being computed are indistinguishable
so it makes sense to perform the same series of operations on every element of the output vector. The input
file below is more approrpriate if the angles are not indistinguishable:

```plumed
a4: ANGLE ATOMS=1,2,3
a5: ANGLE ATOMS=4,5,6
a6: ANGLE ATOMS=7,8,9
PRINT ARG=a4,a5,a6 FILE=colvar
```

*/
//+ENDPLUMEDOC

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
  static void calculateCV( const ColvarInput& cvin, ColvarOutput& cvout );
};

typedef ColvarShortcut<Angle> AngleShortcut;
PLUMED_REGISTER_ACTION(AngleShortcut,"ANGLE")
PLUMED_REGISTER_ACTION(Angle,"ANGLE_SCALAR")
typedef MultiColvarTemplate<Angle> AngleMulti;
PLUMED_REGISTER_ACTION(AngleMulti,"ANGLE_VECTOR")

void Angle::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords(keys);
  keys.setDisplayName("ANGLE");
  keys.add("atoms","ATOMS","the list of atoms involved in this collective variable (either 3 or 4 atoms)");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.setValueDescription("scalar/vector","the ANGLE involving these atoms");
  keys.reset_style("NUMERICAL_DERIVATIVES","hidden");
}

void Angle::parseAtomList( const int& num, std::vector<AtomNumber>& atoms, ActionAtomistic* aa ) {
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

unsigned Angle::getModeAndSetupValues( ActionWithValue* av ) {
  av->addValueWithDerivatives();
  av->setNotPeriodic();
  return 0;
}

Angle::Angle(const ActionOptions&ao):
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
void Angle::calculate() {

  if(pbc) {
    makeWhole();
  }
  ColvarOutput cvout = ColvarOutput::createColvarOutput(value,derivs,this);
  calculateCV( ColvarInput::createColvarInput( 0, getPositions(), this ), cvout );
  setValue( value[0] );
  for(unsigned i=0; i<getPositions().size(); ++i) {
    setAtomsDerivatives( i, cvout.getAtomDerivatives(0,i) );
  }
  setBoxDerivatives( cvout.virial[0] );
}

void Angle::calculateCV( const ColvarInput& cvin, ColvarOutput& cvout ) {
  Vector dij,dik;
  dij=delta(cvin.pos[2],cvin.pos[3]);
  dik=delta(cvin.pos[1],cvin.pos[0]);
  Vector ddij,ddik;
  PLMD::Angle a;
  cvout.values[0]=a.compute(dij,dik,ddij,ddik);
  cvout.derivs[0][0]=ddik;
  cvout.derivs[0][1]=-ddik;
  cvout.derivs[0][2]=-ddij;
  cvout.derivs[0][3]=ddij;
  ColvarInput::setBoxDerivativesNoPbc( cvin, cvout );
}

}
}



