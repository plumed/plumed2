/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The plumed team
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
#include "core/ActionRegister.h"
#include "function/FunctionShortcut.h"
#include "function/FunctionOfScalar.h"
#include "function/FunctionOfVector.h"
#include "EvaluateGridFunction.h"
#include "ActionWithGrid.h"

//+PLUMEDOC GRIDCALC EVALUATE_FUNCTION_FROM_GRID
/*
Calculate the function stored on the input grid at an arbitrary point

This action can be used to evaluate a function that is stored on a grid at an arbitrary point in space via interpolation. For example, the 
following input illustrates how you can use this action to evaulate the square of the distance between atom 1 and atom 2

```plumed
d2: REFERENCE_GRID GRID_MIN=0 GRID_MAX=10 GRID_BIN=20 FUNC=d1*d1 VAR=d1 PERIODIC=NO
d1: DISTANCE ATOMS=1,2
ff: EVALUATE_FUNCTION_FROM_GRID GRID=d2 
PRINT ARG=d1,ff FILE=colvar
```

You can perform the same calculation using [CUSTOM](CUSTOM.md), which will do the calculation exactly and not perform any interpolation. However,
if for some reaon you only have the values of the function at a grid of point EVALUATE_FUNCTION_FROM_GRID may be a useful alternative as you can 
use [REFERENCE_GRID](REFERENCE_GRID.md) to read in the values of the function at your grid of point from a file as shown below:

```plumed
#SETTINGS INPUTFILES=regtest/gridtools/rt-weights-integral/kde.grid

ff: REFERENCE_GRID FILE=regtest/gridtools/rt-weights-integral/kde.grid VALUE=h
aa: ANGLE ATOMS=1,2,3
tt: TORSION ATOMS=1,2,3,4
val: EVALUATE_FUNCTION_FROM_GRID GRID=ff
PRINT ARG=aa,tt,val FILE=colvar
```

Notice that the ARG command is not used with the EVALUATE_FUNCTION_FROM_GRID command to specify the labels of the input arguments at which you 
would like the function to be evaluated.  Instead, yaou only specify the GRID that is being used as input for the interpolation.  The reason that it is not 
necessary to specify the arguments is that PLUMED stores labels for each coordinate axis of a grid.  These labels appear at the top of the 
columns that contain the grid points when you output grids using [DUMPGRID](DUMPGRID.md). Furthermore, in the first of the inputs above the labels of these 
input arguments are specified using the VAR option for [REFERENCE_GRID](REFERENCE_GRID.md). Meanwhile, for the second  the titles of the columsn are 
read in from the input file.  Consequently, PLUMED searches for a value named `d1` when evaluating `ff` in the first input.  For the second input it searches
for values labelled `aa` and `tt` and evaluates `val` from them.   

Notice that input the inputs above the input values are scalars.  The output from EVALUATE_FUNCTION_FROM_GRID is thus also a scalar.  If, however, one uses an 
input such as the one shown below:

```plumed
#SETTINGS INPUTFILES=regtest/gridtools/rt-weights-integral/kde.grid

ff: REFERENCE_GRID FILE=regtest/gridtools/rt-weights-integral/kde.grid VALUE=h
aa: ANGLE ATOMS1=1,2,3 ATOMS2=4,5,6 ATOMS3=7,8,9
tt: TORSION ATOMS1=1,2,3,4 ATOMS2=5,6,7,8 ATOMS3=9,10,11,12
val: EVALUATE_FUNCTION_FROM_GRID GRID=ff
PRINT ARG=aa,tt,val FILE=colvar
```

The value, `val`, that is output by the EVALUATE_FUNCTION_FROM_GRID action will be a vector.  In other words, in the same way you can use scalars, vectors, matrices and even 
functions on grids in the input to the [CUSTOM](CUSTOM.md) command you can also use scalars, vectors, matrices and functions on grid in the input to the EVALUATE_FUNCTION_FROM_GRID
command.

*/
//+ENDPLUMEDOC

//+PLUMEDOC GRIDCALC EVALUATE_FUNCTION_FROM_GRID_SCALAR
/*
Calculate the function stored on the input grid at an arbitrary point

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC GRIDCALC EVALUATE_FUNCTION_FROM_GRID_VECTOR
/*
Calculate the function stored on the input grid for each of the points in the input vector/s

\par Examples

*/
//+ENDPLUMEDOC


namespace PLMD {
namespace gridtools {

class EvaluateFunctionOnGrid : public ActionShortcut {
public:
  static void registerKeywords(Keywords&);
  explicit EvaluateFunctionOnGrid(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(EvaluateFunctionOnGrid,"EVALUATE_FUNCTION_FROM_GRID")
typedef function::FunctionOfScalar<EvaluateGridFunction> ScalarEvalGrid;
PLUMED_REGISTER_ACTION(ScalarEvalGrid,"EVALUATE_FUNCTION_FROM_GRID_SCALAR")
typedef function::FunctionOfVector<EvaluateGridFunction> VectorEvalGrid;
PLUMED_REGISTER_ACTION(VectorEvalGrid,"EVALUATE_FUNCTION_FROM_GRID_VECTOR")

void EvaluateFunctionOnGrid::registerKeywords(Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.addInputKeyword("compulsory","GRID","grid","the name of the grid that we are using to evaluate the function");
  keys.addInputKeyword("optional","ARG","scalar/vector","the arguments that you would like to use when evaluating the function.  If not specified these are determined from the names of the grid dimensions");
  keys.reserve("compulsory","PERIODIC","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO");
  EvaluateGridFunction ii;
  ii.registerKeywords( keys );
  keys.addActionNameSuffix("_SCALAR");
  keys.addActionNameSuffix("_VECTOR");
}

EvaluateFunctionOnGrid::EvaluateFunctionOnGrid(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  // Get the grid that we are evaluating here
  std::vector<std::string> gridn(1);
  parse("GRID",gridn[0]);
  std::vector<Value*> gridv;
  ActionWithArguments::interpretArgumentList( gridn, plumed.getActionSet(), this, gridv );
  // Read the input arguments from the input file
  std::vector<std::string> argn;
  parseVector("ARG",argn);
  // Now get the arguments
  ActionWithGrid* ag = dynamic_cast<ActionWithGrid*>( gridv[0]->getPntrToAction() );
  if( !ag ) {
    error("argument to GRID argument is not a grid");
  }
  if( argn.size()==0 ) {
    std::vector<std::string> argn2( ag->getGridCoordinateNames() );
    argn.resize( argn2.size() );
    for(unsigned i=0; i<argn2.size(); ++i) {
      argn[i]=argn2[i];
    }
  }
  if( argn.size()!=gridv[0]->getRank() ) {
    error("found wrong number of arguments in Evaluate function on grid");
  }
  // Now use this information to create a gridobject
  if( ag->getGridCoordinatesObject().getGridType()=="fibonacci" ) {
    error("cannot interpolate on fibonacci sphere");
  }
  // Now get the actual values we are using
  std::vector<Value*> vals;
  ActionWithArguments::interpretArgumentList( argn, plumed.getActionSet(), this, vals );
  if( vals.size()==0 ) {
    error("found no input arguments to function");
  }
  std::string allargs = gridn[0];
  for(unsigned i=0; i<argn.size(); ++i) {
    allargs += "," + argn[i];
  }
  function::FunctionShortcut<int>::createAction( this, vals, allargs );
}


}
}
