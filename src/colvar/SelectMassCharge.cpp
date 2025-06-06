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
#include "core/ActionAtomistic.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MCOLVAR CHARGE
/*
Get the charges of one or multiple atoms

The following example shows how you can print the charge of atom one:

```plumed
q: CHARGE ATOMS=1
PRINT ARG=q FILE=colvar
```

If you want to output the charges of multiple atoms you would use an input similar to the one below:

```plumed
q: CHARGE ATOMS=1-10
PRINT ARG=q FILE=colvar
```

This input outputs a 10 dimensional vector that contains the charges of the first 10 atoms.

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR MASS
/*
Get the mass of one or multiple atoms

The following example shows how you can print the mass of atom one:

```plumed
m: MASS ATOMS=1
PRINT ARG=m FILE=colvar
```

If you want to output the masses of multiple atoms you would use an input similar to the one below:

```plumed
m: MASS ATOMS=1-10
PRINT ARG=m FILE=colvar
```

This input outputs a 10 dimensional vector that contains the masses of the first 10 atoms.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace colvar {

class SelectMassCharge :
  public ActionAtomistic,
  public ActionWithValue {
public:
  static void registerKeywords( Keywords& keys );
  explicit SelectMassCharge(const ActionOptions&);
// active methods:
  unsigned getNumberOfDerivatives() override {
    return 0;
  }
  void calculate() override;
  void apply() override {}
};

PLUMED_REGISTER_ACTION(SelectMassCharge,"MASS")
PLUMED_REGISTER_ACTION(SelectMassCharge,"CHARGE")

void SelectMassCharge::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  keys.add("atoms","ATOMS","the atom numbers that you would like to store the masses and charges of");
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.setValueDescription("scalar/vector","the " + keys.getDisplayName() + " of the atom/s");
}

SelectMassCharge::SelectMassCharge(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithValue(ao) {
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  log.printf("  getting %s of atoms : ", getName().c_str() );
  for(unsigned i=0; i<atoms.size(); ++i) {
    std::pair<std::size_t,std::size_t> p = getValueIndices( atoms[i] );
    if( getName()=="MASS" && !masv[p.first]->isConstant() ) {
      error("cannot deal with non-constant " + getName() + " values");
    }
    if( getName()=="CHARGE" && !chargev[p.first]->isConstant() ) {
      error("cannot deal with non-constant " + getName() + " values");
    }
    log.printf("%d ", atoms[i].serial() );
  }
  log.printf("\n");
  requestAtoms(atoms);
  std::vector<std::size_t> shape(1);
  if(atoms.size()==1) {
    shape.resize(0);
  } else {
    shape[0] = atoms.size();
  }
  addValue( shape );
  setNotPeriodic();
  getPntrToComponent(0)->setConstant();
}

// calculator
void SelectMassCharge::calculate() {
  Value* myval = getPntrToComponent(0);
  if( getName()=="CHARGES" ) {
    if( !chargesWereSet ) {
      error("cannot determine charges are charges were not set");
    }
    for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
      myval->set( i, getCharge(i) );
    }
  } else {
    for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
      myval->set( i, getMass(i) );
    }
  }
}

}
}



