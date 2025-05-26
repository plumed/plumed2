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
#include "FunctionShortcut.h"
#include "FunctionOfScalar.h"
#include "FunctionOfVector.h"
#include "FunctionOfMatrix.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION COMBINE
/*
Calculate a polynomial combination of a set of other variables.

This action takes in $N_{arg}$ arguments ($x_i$) and computes the following function of them:

\f[
C=\sum_{i=1}^{N_{arg}} c_i (x_i-a_i)^{p_i}
\f]

The $c_i$, $a_i$ and $p_i$ values in this expression are provided through the COEFFICIENTS, PARAMETERS and POWERS keywords respectively.
The following example illustrates how this action can be used to calculate and print the square of the distance between atom 1 and atom 2.

```plumed
d: DISTANCE ATOMS=1,2 COMPONENTS
c: COMBINE ARG=d.x,d.y,d.z COEFFICIENTS=1,1,1 PARAMETERS=0,0,0 POWERS=2,2,2 PERIODIC=NO
PRINT ARG=c FILE=colvar
```

Notice that if the COEFFICIENTS keyword is absent all the $c_i$ values are set equal to 1.
Furthermore, if the PARAMETERS keyword is absent all the $a_i$ values are set equal to 0.
We can thus make use of these defaults and rewrite the input above as:

```plumed
d: DISTANCE ATOMS=1,2 COMPONENTS
c: COMBINE ARG=d.x,d.y,d.z POWERS=2,2,2 PERIODIC=NO
PRINT ARG=c FILE=colvar
```

Notice that we cannot remove the POWERS keyword here as if it is absent all the $p_i$ values are set equal to 1.

## Periodic arguments

The COMBINE action is not able to predict the periodic domain for the function that is computed from the arguments
automatically.  The user is thus forced to specify it explicitly. Use PERIODIC=NO if the resulting variable is not periodic,
and PERIODIC=A,B where A and B are the two boundaries for the periodic domain if the resulting variable is periodic.
The following provides an example where the output from combine has a periodic domain.  In this input we are taking the
cube of a dihedral angle.  The dihedral angle has a periodic domain that runs from $-\pi$ to $\pi$.  The cube of this variable
thus has a periodic domain that runs from $-\pi^3$ to $\pi^3$ as indicated in the following input.

```plumed
t: TORSION ATOMS=1,3,5,7
#c: COMBINE ARG=t POWERS=3 PERIODIC=-31.0062766802998,31.0062766802998
c: COMBINE ARG=t POWERS=3 PERIODIC=-pi^3,pi^3
PRINT ARG=c FILE=colvar
```

## Vector arguments

The two examples in the previous section demonstrate how to use the COMBINE action with scalar arguments.  However, you can also
pass arguments to this action that have a rank greater than zero.  For example, in the example below COMBINE accepts three
vectors with 4 elements in input:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8 COMPONENTS
c: COMBINE ARG=d.x,d.y,d.z POWERS=2,2,2 PERIODIC=NO
PRINT ARG=c FILE=colvar
```

The output from the COMBINE action here is also a vector with four elements. The first element of this vector is the square of the
distance betwen atoms 1 and 2, the second is the square of the distance between atoms 3 and 4 and so on.

The COMBINE action can also take a mixture of scalars and vectors in input as long as the labels of vectors appear before the labels of scalars
in the input to the ARG keyword.  The following input illustrates an COMBINE action that takes vectors and scalars in input.

```plumed
p: CONSTANT VALUE=2
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
c: COMBINE ARG=d,p COEFFICIENTS=4,-1 POWERS=2.5,1 PERIODIC=NO
PRINT ARG=c FILE=colvar
```

The elements of the vector that is calculated by the COMBINE action here are given by:

$$
c_i = 4d_i^{2.5} - 2
$$

The $d_i$ in this expression are the distances that were calculated by the DISTANCE action. The 2 comes
from the scalar `p` that passed to the ARG action in input.  As a single scalar is passed the same number is
used when calculating all 4 elements of the output vector.  In other words, the scalar behaves as if it is a
vector with four components that all have the same value.

Lastly, notice that if you pass a single vector in input to COMBINE as in the following example, the output is still a vector:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
c: COMBINE ARG=d COEFFICIENTS=4 POWERS=2.5 PARAMETERS=0.5 PERIODIC=NO
PRINT ARG=c FILE=colvar
```

To calculate the 4-dimensional output by the COMBINE action here we subtract 0.5 from each of the input distances, raise
the result from this subtraction to the power 2.5 and then multiply the result by 4.

If you want to calculate a linear combination of the elements of a vector using the formula at the top of the page you should use
the [CUSTOM](CUSTOM.md) action to transform all the components of the input vector.  You can then add all the results from these
transformations together using [SUM](SUM.md).

## Matrix arguments

You can also use matrices in the input to the COMBINE action as illustrated below:

```plumed
d: DISTANCE_MATRIX GROUPA=1,2 GROUPB=3,4,5 COMPONENTS
c: COMBINE ARG=d.x,d.y,d.z POWERS=2,2,2 PERIODIC=NO
PRINT ARG=c FILE=colvar
```

The input to the combine action here consists of three $2\times3$ matrices. The output is thus a $2\times3$ matrix that contains the squares of the
distances between the atoms in GROUPA and the atoms in GROUPB.  Notice that all the input matrices must have the same size as the elements of the final
matrix are calculated by applying the formula in the first section of this input to each set of elements to the input matrices in turn.

The input to this action can be a combination of matrices and scalars as long as the labels of the matrices appear before the labels of the scalars in the input for the ARG keyword.  
If your input arguments are an $N\times M$ matrix and a scalar the scalar is treated as if
it is a $N\times M$ matrix which has all its elements equal to the input scalar. You __cannot__ use a mixture of vectors and matrices in the input to this action.

Furthermore, if you pass a single matrix to COMBINE the output is a matrix.  To calculate a linear combination of all the elements in a matrix using the formula at the top of the page you should use
the [CUSTOM](CUSTOM.md) action to transform all the components of the input vector.  You can then add all the results from these transformations together using [SUM](SUM.md).

*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION COMBINE_SCALAR
/*
Calculate a polynomial combination of a set of other variables.

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION COMBINE_VECTOR
/*
Add together the elements of a set of vectors elementwise

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR COMBINE_MATRIX
/*
Calculate the sum of a number of matrices

\par Examples

*/
//+ENDPLUMEDOC

class Combine  {
public:
  std::vector<double> coefficients;
  std::vector<double> parameters;
  std::vector<double> powers;
  std::vector<bool> periodic;
  std::vector<double> max_minus_min;
  std::vector<double> inv_max_minus_min;
  static void registerKeywords(Keywords& keys);
  static void read( Combine& func, ActionWithArguments* action, FunctionOptions& options );
  static void calc( const Combine& func, bool noderiv, const View<const double,helpers::dynamic_extent>& args, FunctionOutput& funcout );
};

typedef FunctionShortcut<Combine> CombineShortcut;
PLUMED_REGISTER_ACTION(CombineShortcut,"COMBINE")
typedef FunctionOfScalar<Combine> ScalarCombine;
PLUMED_REGISTER_ACTION(ScalarCombine,"COMBINE_SCALAR")
typedef FunctionOfVector<Combine> VectorCombine;
PLUMED_REGISTER_ACTION(VectorCombine,"COMBINE_VECTOR")
typedef FunctionOfMatrix<Combine> MatrixCombine;
PLUMED_REGISTER_ACTION(MatrixCombine,"COMBINE_MATRIX")

void Combine::registerKeywords(Keywords& keys) {
  keys.use("PERIODIC");
  keys.add("compulsory","COEFFICIENTS","1.0","the coefficients of the arguments in your function");
  keys.add("compulsory","PARAMETERS","0.0","the parameters of the arguments in your function");
  keys.add("compulsory","POWERS","1.0","the powers to which you are raising each of the arguments in your function");
  keys.addFlag("NORMALIZE",false,"normalize all the coefficients so that in total they are equal to one");
  keys.setValueDescription("scalar/vector/matrix","a linear compbination");
}

void Combine::read( Combine& func, ActionWithArguments* action, FunctionOptions& options ) {
  unsigned nargs = action->getNumberOfArguments();
  ActionWithVector* av=dynamic_cast<ActionWithVector*>(action);
  if(av && av->getNumberOfMasks()>0) {
    nargs = nargs - av->getNumberOfMasks();
  }
  func.coefficients.resize( nargs );
  func.parameters.resize( nargs );
  func.powers.resize( nargs );
  action->parseVector("COEFFICIENTS",func.coefficients);
  if(func.coefficients.size()!=nargs) {
    action->error("Size of COEFFICIENTS array should be the same as number for arguments");
  }
  action->parseVector("PARAMETERS",func.parameters);
  if(func.parameters.size()!=nargs) {
    action->error("Size of PARAMETERS array should be the same as number for arguments");
  }
  action->parseVector("POWERS",func.powers);
  if(func.powers.size()!=nargs) {
    action->error("Size of POWERS array should be the same as number for arguments");
  }

  bool normalize;
  action->parseFlag("NORMALIZE",normalize);
  if(normalize) {
    double n=0.0;
    for(unsigned i=0; i<func.coefficients.size(); i++) {
      n+=func.coefficients[i];
    }
    for(unsigned i=0; i<func.coefficients.size(); i++) {
      func.coefficients[i]*=(1.0/n);
    }
  }

  action->log.printf("  with coefficients:");
  for(unsigned i=0; i<func.coefficients.size(); i++) {
    action->log.printf(" %f",func.coefficients[i]);
  }
  action->log.printf("\n  with parameters:");
  for(unsigned i=0; i<func.parameters.size(); i++) {
    action->log.printf(" %f",func.parameters[i]);
  }
  action->log.printf("\n  and powers:");
  for(unsigned i=0; i<func.powers.size(); i++) {
    action->log.printf(" %f",func.powers[i]);
  }
  action->log.printf("\n");
  // Store periodicity stuff
  func.periodic.resize( nargs, false );
  func.max_minus_min.resize( nargs, 0 );
  func.inv_max_minus_min.resize( nargs, 0 );
  for(unsigned i=0; i<nargs; ++i) {
    if( (action->getPntrToArgument(i))->isPeriodic() ) {
      func.periodic[i] = true;
      std::string min, max;
      (action->getPntrToArgument(i))->getDomain( min, max );
      double dmin, dmax;
      Tools::convert( min, dmin );
      Tools::convert( max, dmax );
      func.max_minus_min[i] = dmax - dmin;
      func.inv_max_minus_min[i] = 1 / func.max_minus_min[i];
    }
  }
}

void Combine::calc( const Combine& func, bool noderiv, const View<const double,helpers::dynamic_extent>& args, FunctionOutput& funcout ) {
  funcout.values[0]=0.0;
  for(unsigned i=0; i<func.coefficients.size(); ++i) {
    double cv = args[i] - func.parameters[i];
    if( func.periodic[i] ) {
      cv = func.max_minus_min[i]*Tools::pbc( cv*func.inv_max_minus_min[i] );
    }
    funcout.values[0] += func.coefficients[i]*pow( cv, func.powers[i] );
    if( !noderiv ) {
      funcout.derivs[0][i] = func.coefficients[i]*func.powers[i]*pow(cv,func.powers[i]-1.0);
    }
  }
}

}
}


