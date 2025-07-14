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
#ifdef __PLUMED_HAS_OPENACC
#define __PLUMED_USE_OPENACC 1
#endif //__PLUMED_HAS_OPENACC
#include "FunctionSetup.h"
#include "FunctionShortcut.h"
#include "FunctionOfScalar.h"
#include "FunctionOfVector.h"
#include "FunctionOfMatrix.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"

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
b: MORE_THAN ARG=d R_0=0.2 D_0=0.0 NN=6 MM=12
```

In the above input the rational switching function that is described in the documentation for [LESS_THAN](LESS_THAN.md) is used to
transform the input distance. However, we would recommend using the following syntax rather than the one above:

```plumed
d: DISTANCE ATOMS=1,2
b: MORE_THAN ARG=d SWITCH={RATIONAL R_0=0.2}
```

as this syntax allows you to use all the switching function options described in the documentation for [LESS_THAN](LESS_THAN.md) here in place of RATIONAL.

## Transforming the square

If you have computed the square of the distance you can use the flag SQUARED to indicate that the input
quantity is the square of the distance as indicated below:

```plumed
d: DISTANCE COMPONENTS ATOMS=1,2
dm: CUSTOM ARG=d.x,d.y,d.z FUNC=x*x+y*y+z*z PERIODIC=NO
s: MORE_THAN ARG=dm SQUARED SWITCH={RATIONAL R_0=0.2}
```

This option can be useful for improving performance by removing the expensive square root operations.

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

## The MASK keyword

Consider the following input:

```plumed
# Fixed virtual atom which serves as the probe volume's center (pos. in nm)
center: FIXEDATOM AT=2.5,2.5,2.5
# Vector in which element i is one if atom i is in sphere of interest and zero otherwise
sphere: INSPHERE ATOMS=1-400 CENTER=center RADIUS={GAUSSIAN D_0=0.5 R_0=0.01 D_MAX=0.52}
# Compute the coordination numbers
adj: CONTACT_MATRIX GROUP=1-400 SWITCH={RATIONAL R_0=0.3 D_MAX=1.0} MASK=sphere
ones: ONES SIZE=400
coord: MATRIX_VECTOR_PRODUCT ARG=adj,ones
# Determine if coordination numbers are more than 4 or not
lt: MORE_THAN ARG=coord MASK=sphere SWITCH={RATIONAL R_0=4.0}
# And calculate how many atoms in the region of interst have a coordination number that is less than four
ltsphere: CUSTOM ARG=lt,sphere FUNC=x*y PERIODIC=NO
cv: SUM ARG=ltsphere PERIODIC=NO
PRINT ARG=cv FILE=colvar
```

This input calculates the total number of atoms that are within a spherical region that is centered on the point $(2.5,2.5,2.5)$ and have a
coordination number that is more than 4.  Notice that to reduce the computational expense we use the MASK keyword in the input to
[CONTACT_MATRIX](CONTACT_MATRIX.md) so PLUMED knows to not bother calculating the coordination numbers of atoms that are not within the spherical region of
interest. Further notice that we have also used this MASK keyword in the input to MORE_THAN to prevent PLUMED from transforming the coordination numbers
that have not been calculated with the switching function to get a further speed up.

Using the MASK keyword in this way is necessary if the input argument is a vector. If the input argument is a matrix then you should never need to use the
MASK keyword. PLUMED will use the sparsity pattern for the input matrix to reduce the number of transformations that are performed within MORE_THAN.

*/
//+ENDPLUMEDOC

class MoreThan {
public:
  bool squared;
#ifdef __PLUMED_USE_OPENACC
  SwitchingFunctionAccelerable switchingFunction;
#else
  SwitchingFunction switchingFunction;
#endif //__PLUMED_USE_OPENACC
  static void registerKeywords( Keywords& keys );
  static void read( MoreThan& func,
                    ActionWithArguments* action,
                    FunctionOptions& options );
  static void calc( const MoreThan& func,
                    bool noderiv,
                    const View<const double>& args,
                    FunctionOutput& funcout );
#ifdef __PLUMED_USE_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1],squared)
    switchingFunction.toACCDevice();
  }
  void removeFromACCDevice() const {
    switchingFunction.removeFromACCDevice();
#pragma acc exit data delete(squared,this[0:1])
  }
#endif //__PLUMED_USE_OPENACC
};

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

void MoreThan::read( MoreThan& func, ActionWithArguments* action, FunctionOptions& options ) {
  options.derivativeZeroIfValueIsZero = true;
  if( action->getNumberOfArguments()!=1 ) {
    ActionWithVector* av = dynamic_cast<ActionWithVector*>( action );
    if( !av || (av && action->getNumberOfArguments()-av->getNumberOfMasks()!=1) ) {
      action->error("should only be one argument to less_than actions");
    }
  }
  if( action->getPntrToArgument(0)->isPeriodic() ) {
    action->error("cannot use this function on periodic functions");
  }


  std::string errors;
  std::string sfinput;
  action->parse("SWITCH", sfinput);
  if(sfinput.length()>0) {
    func.switchingFunction.set(sfinput, errors);
    if( errors.length()!=0 ) {
      action->error("problem reading SWITCH keyword : " + errors );
    }
  } else {
    int nn=6;
    int mm=0;
    double d0=0.0;
    double r0=0.0;
    action->parse("R_0",r0);
    if(r0<=0.0) {
      action->error("R_0 should be explicitly specified and positive");
    }
    action->parse("D_0",d0);
    action->parse("NN",nn);
    action->parse("MM",mm);
    func.switchingFunction.set(nn,mm,r0,d0);
  }
  action->log<<"  using switching function with cutoff "<<func.switchingFunction.description()<<"\n";
  action->parseFlag("SQUARED",func.squared);
  if( func.squared ) {
    action->log<<"  input quantity is square of quantity that switching function acts upon\n";
  }
}

void MoreThan::calc( const MoreThan& func, bool noderiv, const View<const double>& args, FunctionOutput& funcout ) {
  // the presence of NDEBUG seems to be ignored by nvcc...
  // plumed_dbg_assert( args.size()==1 );
  double d;
  if( func.squared ) {
    funcout.values[0] = 1.0 - func.switchingFunction.calculateSqr( args[0], d );
  } else {
    funcout.values[0] = 1.0 - func.switchingFunction.calculate( args[0], d );
  }
  if( !noderiv ) {
    funcout.derivs[0][0] = -args[0]*d;
  }
}

}
}


