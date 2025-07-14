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
#include "tools/HistogramBead.h"
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
b: BETWEEN ARG=d LOWER=0.1 UPPER=0.2 SMEAR=0.5
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

You can achieve the same result as the first example input by using the following input:

```plumed
d: DISTANCE ATOMS=1,2
b: BETWEEN ARG=d SWITCH={GAUSSIAN LOWER=0.1 UPPER=0.2 SMEAR=0.5}
```

The advantage of this syntax, however, is that you can replace the Gaussian kernel in the expression can be replaced by a
triangular Kernel by changing the input to:

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
# Determine if coordination numbers are less than 4 or not
lt: BETWEEN ARG=coord MASK=sphere SWITCH={GAUSSIAN LOWER=4.0 UPPER=6 SMEAR=0.5}
# And calculate how many atoms in the region of interst have a coordination number that is less than four
ltsphere: CUSTOM ARG=lt,sphere FUNC=x*y PERIODIC=NO
cv: SUM ARG=ltsphere PERIODIC=NO
PRINT ARG=cv FILE=colvar
```

This input calculates the total number of atoms that are within a spherical region that is centered on the point $(2.5,2.5,2.5)$ and have a
coordination number that is between 4 and 6.  Notice that to reduce the computational expense we use the MASK keyword in the input to
[CONTACT_MATRIX](CONTACT_MATRIX.md) so PLUMED knows to not bother calculating the coordination numbers of atoms that are not within the spherical region of
interest. Further notice that we have also used this MASK keyword in the input to BETWEEN to prevent PLUMED from transforming the coordination numbers
that have not been calculated with the switching function to get a further speed up.

Using the MASK keyword in this way is necessary if the input argument is a vector. If the input argument is a matrix then you should never need to use the
MASK keyword. PLUMED will use the sparsity pattern for the input matrix to reduce the number of transformations that are performed within BETWEEN.

*/
//+ENDPLUMEDOC

struct Between {
  HistogramBead hist{HistogramBead::KernelType::gaussian,0.0,1.0,0.5};
  static void registerKeywords( Keywords& keys );
  static void read( Between& func, ActionWithArguments* action,
                    FunctionOptions& options );
  static void calc( const Between& func,
                    bool noderiv,
                    const View<const double> args,
                    FunctionOutput& funcout );
#ifdef __PLUMED_USE_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1])
    hist.toACCDevice();
  }
  void removeFromACCDevice() const {
    hist.removeFromACCDevice();
#pragma acc exit data delete(order, this[0:1])
  }
#endif // __PLUMED_USE_OPENACC
};

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

void Between::read( Between& func, ActionWithArguments* action, FunctionOptions& options ) {
  options.derivativeZeroIfValueIsZero = true;
  if( action->getNumberOfArguments()!=1 ) {
    ActionWithVector* av = dynamic_cast<ActionWithVector*>( action );
    if( !av || (av && action->getNumberOfArguments()-av->getNumberOfMasks()!=1) ) {
      action->error("should only be one argument to less_than actions");
    }
  }

  std::string str_min, str_max;
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
  func.hist.set( hinput, errors );
  if( errors.size()!=0 ) {
    action->error( errors );
  }
  action->log.printf("  %s \n", func.hist.description().c_str() );

  if( !isPeriodic ) {
    func.hist.isNotPeriodic();
  } else {
    double min;
    double max;
    Tools::convert( str_min, min );
    Tools::convert( str_max, max );
    func.hist.isPeriodic( min, max );
  }
}

void Between::calc( const Between& func,
                    bool noderiv,
                    const View<const double> args,
                    FunctionOutput& funcout ) {
  // the presence of NDEBUG seems to be ignored by nvcc...
  // plumed_dbg_assert( args.size()==1 );
  double deriv;
  funcout.values[0] = func.hist.calculate( args[0], deriv );
  if( !noderiv ) {
    funcout.derivs[0][0] = deriv;
  }
}

}
}


