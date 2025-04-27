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
#include "Sum.h"
#include "FunctionShortcut.h"
#include "FunctionOfScalar.h"
#include "FunctionOfVector.h"
#include "FunctionOfMatrix.h"
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

//+PLUMEDOC COLVAR SUM_VECTOR
/*
Calculate the sum of the elements in a vector

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR SUM_SCALAR
/*
Calculate the SUM of the set of input scalars

\par Examples

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

//+PLUMEDOC COLVAR MEAN_SCALAR
/*
Calculate the arithmetic mean of the set of input scalars

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR MEAN_VECTOR
/*
Calculate the arithmetic mean of the elements in a vector

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR SUM_MATRIX
/*
Sum all the elements in a matrix

\par Examples

*/
//+ENDPLUMEDOC


namespace PLMD {
namespace function {

typedef FunctionShortcut<Sum> SumShortcut;
PLUMED_REGISTER_ACTION(SumShortcut,"SUM")
PLUMED_REGISTER_ACTION(SumShortcut,"MEAN")
typedef FunctionOfScalar<Sum> ScalarSum;
PLUMED_REGISTER_ACTION(ScalarSum,"SUM_SCALAR")
PLUMED_REGISTER_ACTION(ScalarSum,"MEAN_SCALAR")
typedef FunctionOfVector<Sum> VectorSum;
PLUMED_REGISTER_ACTION(VectorSum,"SUM_VECTOR")
PLUMED_REGISTER_ACTION(VectorSum,"MEAN_VECTOR")
typedef FunctionOfMatrix<Sum> MatrixSum;
PLUMED_REGISTER_ACTION(MatrixSum,"SUM_MATRIX")

void Sum::registerKeywords( Keywords& keys ) {
  keys.use("PERIODIC");
  keys.setValueDescription("scalar","the sum");
  keys.add("hidden","MASKED_INPUT_ALLOWED","turns on that you are allowed to use masked inputs ");
}

void Sum::read( ActionWithArguments* action ) {
  if( action->getNumberOfArguments()!=1 ) {
    action->error("should only be one argument to sum actions");
  }
}

void Sum::setPrefactor( ActionWithArguments* action, const double pref ) {
  if(action->getName().find("MEAN")!=std::string::npos) {
    prefactor = pref / (action->getPntrToArgument(0))->getNumberOfValues();
  } else {
    prefactor = pref;
  }
}

bool Sum::zeroRank() const {
  return true;
}

void Sum::calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const {
  vals[0]=prefactor*args[0];
  derivatives(0,0)=prefactor;
}

}
}
