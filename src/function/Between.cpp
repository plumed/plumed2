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
#include "Between.h"
#include "FunctionShortcut.h"
#include "FunctionOfScalar.h"
#include "FunctionOfVector.h"
#include "FunctionOfMatrix.h"
#include "core/ActionRegister.h"

#include <cmath>

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION BETWEEN
/*
Use a switching function to determine how many of the input variables are within a certain range.

This action takes one argument, $s$ and evaluates the following function:

$$
w(s) = \int_a^b K\left( \frac{x - s}{\sigma} \right) \textrm{d}x
$$

In this equation $K$ is a symmetric function that must integrate to one and that is often
called a [kernel function](https://en.wikipedia.org/wiki/Kernel_(statistics)) and $\sigma$ is a smearing parameter.
The above function can be used to evaluate whether $s$ is between $a$ and $b$. The advantage of using the function above is that
the resulting quantity has continuous derivatives. It can thus be used when calculating a collective variable upon which a simulation
bias will be applied.

The following example, shows how we can apply the function above on the instantaneous value of the distance between atom 1 and 2.
The BETWEEN action here is used to determine whether the input distance is between 0.1 and 0.2 nm.

```plumed
d: DISTANCE ATOMS=1,2
b: BETWEEN ARG=d SWITCH={GAUSSIAN LOWER=0.1 UPPER=0.2 SMEAR=0.5}
```

## The Kernel function

The $\sigma$ values in the expressions above is calculated from the parameters $a$, $b$ and $s$ that are provided using the
`LOWER`, `UPPER` and `SMEAR` parameters respectively using the following function:

$$
\sigma = s(b-a)
$$

Also note that the Kernel function $K$ that is used in the first example is a Gaussian. The actual integral that is evaluated is thus:

$$
w(s) = \frac{1}{\sqrt{2\pi}\sigma} \int_a^b \exp\left( -\frac{(x-s)^2}{2\sigma^2}\right) \textrm{d}x
$$

The Gaussian kernel in this expression can be replaced by a triangular Kernel by changing the input to:

```plumed
d: DISTANCE ATOMS=1,2
b: BETWEEN ARG=d SWITCH={TRIANGULAR LOWER=0.1 UPPER=0.2 SMEAR=0.5}
```

With this input the integral is then evaluated using:

$$
w(s) = \frac{1}{2\sigma} \int_a^b 1 - H\left( \frac{x-s}{\sigma} \right) \textrm{d}x
$$

where:

$$
H(x) = \begin{cases}
x & \textrm{if} \quad x<1 \\
0 & \textrm{otherwise}
\end{cases}
$$

## Non rank zero arguments

Instead of passing a single scalar in the input to the `BETWEEN` action you can pass a single vector as shown here:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
b: BETWEEN ARG=d SWITCH={GAUSSIAN LOWER=0.1 UPPER=0.2 SMEAR=0.5}
```

The input to the `BETWEEN` action here is a vector with four elements. The output from the action `b` is similarly
a vector with four elements. In calculating the elements of this vector PLUMED applies the function described in the previous
section on each of the distances in turn. The first element of `b` thus tells you if the distance between atoms 1 and 2 is between
0.1 and 0.2 nm, the second element tells you if the distance between atoms 3 and 4 is between 0.1 and 0.2 nm and so on.

You can use the commands in the above example input to evaluate the number of distances that are within the range of interest as follows:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
b: BETWEEN ARG=d SWITCH={GAUSSIAN LOWER=0.1 UPPER=0.2 SMEAR=0.5}
s: SUM ARG=b PERIODIC=NO
PRINT ARG=s FILE=colvar
```

The final scalar that is output here is evaluated using the following summation:

$$
s = \sum_i \int_a^b \left( \frac{x - d_i}{\sigma} \right) \textrm{d}x
$$

where the sum over $i$ here runs over the four distances in the above expression. This scalar tells you the number of distances that are
between 0.1 and 0.2.

Notice that you can do something similar with a matrix as input as shown below:

```plumed
d: DISTANCE_MATRIX GROUPA=1-10 GROUPB=11-20
b: BETWEEN ARG=d SWITCH={GAUSSIAN LOWER=0.1 UPPER=0.2 SMEAR=0.5}
s: SUM ARG=b PERIODIC=NO
PRINT ARG=s FILE=colvar
```

This input tells PLUMED to calculate the 100 distances between the atoms in the two input groups. The final value that is printed to the colvar file then
tells you how many of these distances are between 0.1 and 0.2 nm.

*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION BETWEEN_VECTOR
/*
Use a switching function to determine how many of the input components are within a certain range

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC COLVAR BETWEEN_MATRIX
/*
Transform all the elements of a matrix using a switching function that is oen when the input value is within a particular range

\par Examples

*/
//+ENDPLUMEDOC

typedef FunctionShortcut<Between> BetweenShortcut;
PLUMED_REGISTER_ACTION(BetweenShortcut,"BETWEEN")
typedef FunctionOfScalar<Between> ScalarBetween;
PLUMED_REGISTER_ACTION(ScalarBetween,"BETWEEN_SCALAR")
typedef FunctionOfVector<Between> VectorBetween;
PLUMED_REGISTER_ACTION(VectorBetween,"BETWEEN_VECTOR")
typedef FunctionOfMatrix<Between> MatrixBetween;
PLUMED_REGISTER_ACTION(MatrixBetween,"BETWEEN_MATRIX")

void Between::registerKeywords(Keywords& keys) {
  keys.add("compulsory","LOWER","the lower boundary for this particular bin");
  keys.add("compulsory","UPPER","the upper boundary for this particular bin");
  keys.add("compulsory","SMEAR","0.5","the ammount to smear the Gaussian for each value in the distribution");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous function defined above. "
           "The following provides information on the \\ref histogrambead that are available. "
           "When this keyword is present you no longer need the LOWER, UPPER, SMEAR and KERNEL keywords.");
  keys.setValueDescription("scalar/vector/matrix","a function that is one if the input falls within a particular range and zero otherwise");
}

void Between::read( ActionWithArguments* action ) {
  if( action->getNumberOfArguments()!=1 ) {
    ActionWithVector* av = dynamic_cast<ActionWithVector*>( action );
    if( !av || (av && action->getNumberOfArguments()-av->getNumberOfMasks()!=1) ) {
      action->error("should only be one argument to less_than actions");
    }
  }
  if( action->getNumberOfArguments()!=1 ) {
    action->error("should only be one argument to between actions");
  }

  std::string str_min, str_max, tstr_min, tstr_max;
  bool isPeriodic = action->getPntrToArgument(0)->isPeriodic();
  if( isPeriodic ) {
    action->getPntrToArgument(0)->getDomain( str_min, str_max );
  }

  std::string hinput;
  action->parse("SWITCH",hinput);
  if(hinput.length()==0) {
    std::string low, up, sme;
    action->parse("LOWER",low);
    action->parse("UPPER",up);
    action->parse("SMEAR",sme);
    hinput = "GAUSSIAN LOWER=" + low + " UPPER=" + up + " SMEAR=" + sme;
  }
  std::string errors;
  hist.set( hinput, errors );
  if( errors.size()!=0 ) {
    action->error( errors );
  }
  action->log.printf("  %s \n", hist.description().c_str() );

  if( !isPeriodic ) {
    hist.isNotPeriodic();
  } else {
    double min;
    Tools::convert( str_min, min );
    double max;
    Tools::convert( str_max, max );
    hist.isPeriodic( min, max );
  }
}

void Between::calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const {
  plumed_dbg_assert( args.size()==1 );
  vals[0] = hist.calculate( args[0], derivatives(0,0) );
}

}
}


