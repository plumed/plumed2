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
#include "FunctionSetup.h"
#include "FunctionShortcut.h"
#include "FunctionOfScalar.h"
#include "FunctionOfVector.h"
#include "core/ActionRegister.h"
#include <array>
#include <cmath>

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION BESSEL
/*
Calculate the value of a Bessel function.

This action can be used to evaluate [Bessel functions](https://en.wikipedia.org/wiki/Bessel_function).
This was action is needed because the modified Bessel function of the first kind of order 0 is needed in
order to calculate the [von-Misses distribution](https://en.wikipedia.org/wiki/Von_Mises_distribution)
that is used in the implementation of [KERNEL](KERNEL.md).  You can thus only use this function to evaluate the
modified Bessel function of the first kind of order 0.  You can evaluate the value of this Bessel function
at $x=1$ by using the following input:

```plumed
c: CONSTANT VALUE=1
b: BESSEL ARG=c ORDER=0
```

Notice that you can also use a vector in the input for this action as illustrated below:

```plumed
c: CONSTANT VALUES=1,1,5,2
b: BESSEL ARG=c ORDER=0
```

The value output by BESSEL in this case is also a vector with three components.  These components give the values
of the Bessel function at 1.0, 1.5 and 2.0 respectively.

## The MASK keyword

As this action only takes one argument in input the MASK keyword action works in the way that is described in the documentation for
[LESS_THAN](LESS_THAN.md).  This keyword is only used if the input argument is a vector.  You pass to it a vector that contains ones
and zeros and that has the same length as the vector that is passed in the input to ARG. If the $i$th element of the vector that was passed
in the input to MASK is one then the corresponding element of the output vector is evaluated.  If it is equal to zero then that element is not
evaluated.  Consequently in the example input shown below only the 1st and 4th element of the output vector are evaluated.

```plumed
m: CONSTANT VALUES=1,0,0,1
c: CONSTANT VALUES=1,1,5,2
b: BESSEL ARG=c MASK=m ORDER=0
```

*/
//+ENDPLUMEDOC

class Bessel {
  static constexpr std::array<double,30> A{-4.41534164647933937950E-18,
           3.33079451882223809783E-17,
           -2.43127984654795469359E-16,
           1.71539128555513303061E-15,
           -1.16853328779934516808E-14,
           7.67618549860493561688E-14,
           -4.85644678311192946090E-13,
           2.95505266312963983461E-12,
           -1.72682629144155570723E-11,
           9.67580903537323691224E-11,
           -5.18979560163526290666E-10,
           2.65982372468238665035E-9,
           -1.30002500998624804212E-8,
           6.04699502254191894932E-8,
           -2.67079385394061173391E-7,
           1.11738753912010371815E-6,
           -4.41673835845875056359E-6,
           1.64484480707288970893E-5,
           -5.75419501008210370398E-5,
           1.88502885095841655729E-4,
           -5.76375574538582365885E-4,
           1.63947561694133579842E-3,
           -4.32430999505057594430E-3,
           1.05464603945949983183E-2,
           -2.37374148058994688156E-2,
           4.93052842396707084878E-2,
           -9.49010970480476444210E-2,
           1.71620901522208775349E-1,
           -3.04682672343198398683E-1,
           6.76795274409476084995E-1};
  static constexpr std::array<double,25>B {-7.23318048787475395456E-18,
           -4.83050448594418207126E-18,
           4.46562142029675999901E-17,
           3.46122286769746109310E-17,
           -2.82762398051658348494E-16,
           -3.42548561967721913462E-16,
           1.77256013305652638360E-15,
           3.81168066935262242075E-15,
           -9.55484669882830764870E-15,
           -4.15056934728722208663E-14,
           1.54008621752140982691E-14,
           3.85277838274214270114E-13,
           7.18012445138366623367E-13,
           -1.79417853150680611778E-12,
           -1.32158118404477131188E-11,
           -3.14991652796324136454E-11,
           1.18891471078464383424E-11,
           4.94060238822496958910E-10,
           3.39623202570838634515E-9,
           2.26666899049817806459E-8,
           2.04891858946906374183E-7,
           2.89137052083475648297E-6,
           6.88975834691682398426E-5,
           3.36911647825569408990E-3,
           8.04490411014108831608E-1};
public:
  unsigned order{0};
  static void registerKeywords( Keywords& keys );
  explicit Bessel();
  static void read( Bessel& func, ActionWithArguments* action, FunctionOptions& options );
  template <size_t N>
  static double chbevl(double x, const std::array<double,N>& array);
  static void calc( const Bessel& func, bool noderiv, View<const double> args, FunctionOutput& funcout );
};

typedef FunctionShortcut<Bessel> BesselShortcut;
PLUMED_REGISTER_ACTION(BesselShortcut,"BESSEL")
typedef FunctionOfScalar<Bessel> ScalarBessel;
PLUMED_REGISTER_ACTION(ScalarBessel,"BESSEL_SCALAR")
typedef FunctionOfVector<Bessel> VectorBessel;
PLUMED_REGISTER_ACTION(VectorBessel,"BESSEL_VECTOR")

void Bessel::registerKeywords(Keywords& keys) {
  keys.add("compulsory","ORDER","0","the order of Bessel function to use.  Can only be zero at the moment.");
  keys.setValueDescription("scalar/vector","the value of the bessel function");
}

Bessel::Bessel() {
}

void Bessel::read( Bessel& func, ActionWithArguments* action, FunctionOptions& options ) {
  if( action->getNumberOfArguments()!=1 ) {
    ActionWithVector* av = dynamic_cast<ActionWithVector*>( action );
    if( !av || (av && action->getNumberOfArguments()-av->getNumberOfMasks()!=1) ) {
      action->error("should only be one argument to less_than actions");
    }
  }
  if( action->getPntrToArgument(0)->isPeriodic() ) {
    action->error("cannot use this function on periodic functions");
  }
  action->parse("ORDER",func.order);
  action->log.printf("  computing %dth order bessel function \n", func.order );
  if( func.order!=0 ) {
    action->error("only zero order bessel function is implemented");
  }
}

template <size_t N>
double Bessel::chbevl(const double x, const std::array<double,N>& array) {
  double b0, b1, b2;

  b0 = array[0];
  b1 = 0.0;
  b2 = b1 ;

  for(unsigned index = 1 ; index < N ; index++) {
    b2 = b1;
    b1 = b0;
    b0 = x * b1 - b2 + array[index];
  }
  return (0.5 * (b0 - b2));
}


void Bessel::calc( const Bessel& func,
                   bool noderiv,
                   const View<const double> args,
                   FunctionOutput& funcout ) {
// NDEBUG seems to be ignored by nvcc...
  plumed_dbg_assert( args.size()==1 );
  if( !noderiv ) {
    plumed_merror("derivatives have not been implemented for this function");
  }
  //I trust the parser into avoiding geting here with order!=0
  double x = fabs(args[0]);
  if (x <= 8.0) {
    double y = (x / 2.0) - 2.0;
    funcout.values[0] =  chbevl(y, Bessel::A) ;
    return;
  }
  funcout.values[0] = chbevl(32.0 / x - 2.0, Bessel::B) / sqrt(x) ;
  ///TODO: higher order bessel functions

}

}
}


