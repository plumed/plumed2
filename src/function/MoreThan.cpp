/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "MoreThan.h"
#include "FunctionShortcut.h"
#include "FunctionOfScalar.h"
#include "FunctionOfVector.h"
#include "FunctionOfMatrix.h"
#include "core/ActionRegister.h"

#include <cmath>

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION MORE_THAN
/*
Use a switching function to determine how many of the input variables are more than a certain cutoff.

This action takes one argument, $r$ and evaluates the following function:

$$
w(s) = 1 - s(r)
$$

In this equation $s(r)$ is one of the switching functions described in the documentation for the action [LESS_THAN](LESS_THAN.md).
The output value $w$ is thus a number between 0 and 1 that tells you if the input value is greater than some cutoff.  Furthermore,
the value of $w$ smoothly from zero to one as the input value $r$ crosses the threshold of interest so any function of this value 
is differentiable.

The following example, shows how we can apply the function above on the instantaneous value of the distance between atom 1 and 2. 
The MORE_THAN action here is used to determine whether the input distance is greater than 0.2 nm.

```plumed
d: DISTANCE ATOMS=1,2
b: MORE_THAN ARG=d SWITCH={RATIONAL R_0=0.2}
```

You can use all the switching function options described in the documentation for [LESS_THAN](LESS_THAN.md) here in place of RATIONAL.

## Non rank zero arguments

Instead of passing a single scalar in the input to the `MORE_THAN` action you can pass a single vector as shown here:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
b: MORE_THAN ARG=d SWITCH={RATIONAL R_0=0.2}
```

The input to the `MORE_THAN` action here is a vector with four elements. The output from the action `b` is similarly 
a vector with four elements. In calculating the elements of this vector PLUMED applies the function described in the previous 
section on each of the distances in turn. The first element of `b` thus tells you if the distance between atoms 1 and 2 is between
greater than 0.2 nm, the second element tells you if the distance between atoms 3 and 4 is greater than 0.2 nm and so on.

You can use the commands in the above example input to evaluate the number of distances that greater than a threshold as follows:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
b: MORE_THAN ARG=d SWITCH={RATIONAL R_0=0.2}
s: SUM ARG=b PERIODIC=NO
PRINT ARG=s FILE=colvar
```

The final scalar that is output here is evaluated using the following summation:

$$
s = \sum_i 1 - s(d_i)
$$

where the sum over $i$ here runs over the four distances in the above expression. This scalar tells you the number of distances that are 
more than 0.2 nm.

Notice that you can do something similar with a matrix as input as shown below:

```plumed
d: DISTANCE_MATRIX GROUPA=1-10 GROUPB=11-20
b: MORE_THAN ARG=d SWITCH={RATIONAL R_0=0.2}
s: SUM ARG=b PERIODIC=NO
PRINT ARG=s FILE=colvar
```

This input tells PLUMED to calculate the 100 distances between the atoms in the two input groups. The final value that is printed to the colvar file then 
tells you how many of these distances are greater than 0.2 nm.

*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION MORE_THAN_VECTOR
/*
Use a switching function to determine how many of elements in the input vector are more than a certain cutoff.

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR MORE_THAN_MATRIX
/*
Transform all the elements of a matrix using a switching function that is one when the input value is larger than a threshold

\par Examples

*/
//+ENDPLUMEDOC

typedef FunctionShortcut<MoreThan> MoreThanShortcut;
PLUMED_REGISTER_ACTION(MoreThanShortcut,"MORE_THAN")
typedef FunctionOfScalar<MoreThan> ScalarMoreThan;
PLUMED_REGISTER_ACTION(ScalarMoreThan,"MORE_THAN_SCALAR")
typedef FunctionOfVector<MoreThan> VectorMoreThan;
PLUMED_REGISTER_ACTION(VectorMoreThan,"MORE_THAN_VECTOR")
typedef FunctionOfMatrix<MoreThan> MatrixMoreThan;
PLUMED_REGISTER_ACTION(MatrixMoreThan,"MORE_THAN_MATRIX")

void MoreThan::registerKeywords(Keywords& keys) {
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.addFlag("SQUARED",false,"is the input quantity the square of the value that you would like to apply the switching function to");
  keys.setValueDescription("scalar/vector/matrix","a function that is one if the if the input is more than a threshold");
}

void MoreThan::read( ActionWithArguments* action ) {
  if( action->getNumberOfArguments()!=1 ) action->error("should only be one argument to more_than actions");
  if( action->getPntrToArgument(0)->isPeriodic() ) action->error("cannot use this function on periodic functions");


  std::string sw,errors;
  action->parse("SWITCH",sw);
  if(sw.length()>0) {
    switchingFunction.set(sw,errors);
    if( errors.length()!=0 ) action->error("problem reading SWITCH keyword : " + errors );
  } else {
    int nn=6; int mm=0; double d0=0.0; double r0=0.0; action->parse("R_0",r0);
    if(r0<=0.0) action->error("R_0 should be explicitly specified and positive");
    action->parse("D_0",d0); action->parse("NN",nn); action->parse("MM",mm);
    switchingFunction.set(nn,mm,r0,d0);
  }
  action->log<<"  using switching function with cutoff "<<switchingFunction.description()<<"\n";
  action->parseFlag("SQUARED",squared);
  if( squared ) action->log<<"  input quantity is square of quantity that switching function acts upon\n";
}

void MoreThan::calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const {
  plumed_dbg_assert( args.size()==1 );
  if( squared ) vals[0] = 1.0 - switchingFunction.calculateSqr( args[0], derivatives(0,0) );
  else vals[0] = 1.0 - switchingFunction.calculate( args[0], derivatives(0,0) );
  derivatives(0,0) = -args[0]*derivatives(0,0);
}

}
}


