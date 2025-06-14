/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
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
#include "ActionWithVirtualAtom.h"
#include "core/ActionRegister.h"
#include "tools/Vector.h"
#include "tools/Exception.h"

namespace PLMD {
namespace vatom {

//+PLUMEDOC VATOM FIXEDATOM
/*
Add a virtual atom in a fixed position.

This action creates [a virtual atom](specifying_atoms.md) at a fixed position.  The example input below illustrates
how this idea can be used to compute the angle between the vector connecting atoms 15 and 20 and the z axis and how this
quantity can then be restrained so that the angle stays close to zero.

```plumed
a: FIXEDATOM AT=0,0,0
b: FIXEDATOM AT=0,0,1
an: ANGLE ATOMS=a,b,15,20
RESTRAINT ARG=an AT=0.0 KAPPA=100.0
```

By default PLUMED assumes that any coordinates specified using the AT keyword specified are the cartesian coordinates of the fixed atom.
However, if you use the SCALED_COMPONENTS flag as shown below:

```plumed
a: FIXEDATOM AT=0.25,0.25,0.25 SCALED_COMPONENTS
DUMPATOMS ATOMS=a FILE=vatom.xyz
```

the coordinates specified using the AT keyword are interpretted as scaled coordinates. The positions output to the `vatom.xyz` file
in the input above is thus obtained by multiplying the input vector by the cell vectors on every step.  The position of the atom `a`
thus changes as the box size changes.

It is also possible to assign a predefined charge or mass to the atom by using the `SET_MASS` and `SET_CHARGE` keywords.

!!! caution ""

    This action, like [POSITION](POSITION.md) is not invariant for translation of the system so adding a force on it can cause trouble.

The problem is that the vector connecting any atom and a virtual atom created using the FIXEDATOM atoms command is not invariant to translation.
However, if, as has been done in the following example input, one first aligns atoms to a reference using [FIT_TO_TEMPLATE](FIT_TO_TEMPLATE.md),
then it is safe to add further fixed atoms without breaking translational invariance.

```plumed
#SETTINGS INPUTFILES=regtest/basic/rt63/align.pdb
FIT_TO_TEMPLATE STRIDE=1 REFERENCE=regtest/basic/rt63/align.pdb TYPE=SIMPLE
a: FIXEDATOM AT=10,20,30
d: DISTANCE ATOMS=a,20
PRINT ARG=d FILE=colvar
```

*/
//+ENDPLUMEDOC


class FixedAtom:
  public ActionWithVirtualAtom {
  Vector coord;
  std::vector<Tensor> deriv;
  double mass,charge;
  bool scaled_components;
public:
  explicit FixedAtom(const ActionOptions&ao);
  void calculate() override;
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(FixedAtom,"FIXEDATOM")

void FixedAtom::registerKeywords(Keywords& keys) {
  ActionWithVirtualAtom::registerKeywords(keys);
  keys.add("compulsory","AT","coordinates of the virtual atom");
  keys.add("compulsory","SET_MASS","1","mass of the virtual atom");
  keys.add("compulsory","SET_CHARGE","0","charge of the virtual atom");
  keys.addFlag("SCALED_COMPONENTS",false,"use scaled components");
  keys.reset_style("ATOMS","hidden");
}

FixedAtom::FixedAtom(const ActionOptions&ao):
  Action(ao),
  ActionWithVirtualAtom(ao) {
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if(atoms.size()!=0) {
    error("ATOMS should be empty");
  }

  parseFlag("SCALED_COMPONENTS",scaled_components);

  std::vector<double> at;
  parseVector("AT",at);
  if(at.size()!=3) {
    error("AT should be a list of three real numbers");
  }

  parse("SET_MASS",mass);
  parse("SET_CHARGE",charge);

  coord[0]=at[0];
  coord[1]=at[1];
  coord[2]=at[2];

  checkRead();
  log<<"  AT position "<<coord[0]<<" "<<coord[1]<<" "<<coord[2]<<"\n";
  if(scaled_components) {
    log<<"  position is in scaled components\n";
  }
}

void FixedAtom::calculate() {
  deriv.resize(getNumberOfAtoms());
  if(scaled_components) {
    setPosition(getPbc().scaledToReal(coord));
  } else {
    setPosition(coord);
  }
  setMass(mass);
  setCharge(charge);
  setAtomsDerivatives(deriv);
// Virial contribution
  if(!scaled_components) {
    setBoxDerivativesNoPbc();
  }
// notice that with scaled components there is no additional virial contribution
}

}
}
