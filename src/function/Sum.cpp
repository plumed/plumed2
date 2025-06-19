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
#include "core/ActionWithArguments.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"

//+PLUMEDOC COLVAR SUM
/*
Calculate the sum of the arguments

This action takes a single vector, a single matrix or a single grid in input. The output is a scalar
that contains the sum of all the elements in the input vector/matrix/grid. This action is
very useful if you want to do calculations like the one illustrated in this example:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
b: LESS_THAN ARG=d SWITCH={RATIONAL R_0=0.1}
s: SUM ARG=b PERIODIC=NO
PRINT ARG=s FILE=colvar
```

This example calculates and outputs the number of distances that are less than 0.1 nm.

You can do something similar by summing the elements of a matrix as shown below:

```plumed
c: CONTACT_MATRIX GROUP=1-100 SWITCH={RATIONAL R_0=0.1}
s: SUM ARG=c PERIODIC=NO
PRINT ARG=s FILE=colvar
```

If you want to sum all the elements in a grid you can. However, if the input is a function on a grid
it is more likely that you will want to compute the integral using [INTEGRATE_GRID](INTEGRATE_GRID.md).

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR MEAN
/*
Calculate the arithmetic mean of the elements in a vector

This action takes a single vector or a single matrix in input.  The output is a scalar
that contains the mean of all the elements in the input vector/matrix. This action is
useful if you want do to calculations like the one illustrated in this example:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
m: MEAN ARG=d PERIODIC=NO
PRINT ARG=m FILE=colvar
```

The output from the MEAN action here is the average over the four distances that are evaluated by the
DISTANCE action.  Notice that you can also do the calculation above using:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
s: SUM ARG=d PERIODIC=NO
m: CUSTOM ARG=s FUNC=x/4 PERIODIC=NO
PRINT ARG=m FILE=colvar
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace function {

class Sum :
  public ActionWithValue,
  public ActionWithArguments {
private:
  bool ismatrix;
public:
  static void registerKeywords(Keywords& keys);
  explicit Sum(const ActionOptions&);
  unsigned getNumberOfDerivatives() override ;
  void calculate() override ;
  void apply() override ;
};

PLUMED_REGISTER_ACTION(Sum,"SUM")
PLUMED_REGISTER_ACTION(Sum,"MEAN")

void Sum::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.addInputKeyword("compulsory","ARG","vector/matrix/grid","the vector/matrix/grid whose elements shuld be added together");
  keys.add("compulsory","PERIODIC","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO");
  keys.setValueDescription("scalar","the " + keys.getDisplayName() + " of the elements in the input value");
  keys.add("hidden","MASKED_INPUT_ALLOWED","turns on that you are allowed to use masked inputs ");
  keys.reset_style("NUMERICAL_DERIVATIVES","hidden");
}

Sum::Sum(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  ismatrix(false) {

  if( getNumberOfArguments()!=1 ) {
    error("should only be one argument in input for this action");
  }
  ismatrix = getPntrToArgument(0)->getRank()==2 && !getPntrToArgument(0)->hasDerivatives();

  addValue();
  std::vector<std::string> period;
  parseVector("PERIODIC",period);
  if( period.size()==1 ) {
    if( period[0]!="NO") {
      error("input to PERIODIC keyword does not make sense");
    }
    setNotPeriodic();
  } else if( period.size()==2 ) {
    setPeriodic( period[0], period[1] );
  } else {
    error("input to PERIODIC keyword does not make sense");
  }
}

unsigned Sum::getNumberOfDerivatives() {
  return 0;
}

void Sum::calculate() {
  double myval = 0;
  const Value* myarg = getPntrToArgument(0);
  if( ismatrix ) {
    unsigned ncols = myarg->getNumberOfColumns();
    for(unsigned i=0; i<myarg->getShape()[0]; ++i) {
      for(unsigned j=0; j<myarg->getRowLength(i); ++j) {
        myval += myarg->get(i*ncols + j, false);
      }
    }
  } else {
    unsigned nargs = myarg->getNumberOfStoredValues();
    for(unsigned i=0; i<nargs; ++i) {
      myval += myarg->get(i);
    }
  }
  if( getName()=="MEAN" ) {
    getPntrToComponent(0)->set( myval / myarg->getNumberOfValues() );
  } else {
    getPntrToComponent(0)->set( myval );
  }
}

void Sum::apply() {
  if( !getPntrToComponent(0)->forcesWereAdded() ) {
    return ;
  }

  Value* myarg = getPntrToArgument(0);
  double force = getPntrToComponent(0)->getForce(0);
  if( getName()=="MEAN" ) {
    force = force / myarg->getNumberOfValues();
  }
  if( ismatrix ) {
    unsigned ncols = myarg->getNumberOfColumns();
    for(unsigned i=0; i<myarg->getShape()[0]; ++i) {
      for(unsigned j=0; j<myarg->getRowLength(i); ++j) {
        myarg->addForce(i*ncols + j, force, false);
      }
    }
  } else {
    unsigned nargs = myarg->getNumberOfStoredValues();
    for(unsigned i=0; i<nargs; ++i) {
      myarg->addForce( i, force );
    }
  }
}

}


}
