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
#include "function/FunctionShortcut.h"
#include "function/FunctionOfScalar.h"
#include "function/FunctionOfVector.h"
#include "core/ActionRegister.h"
#include "function/FunctionSetup.h"

#include <cmath>

namespace PLMD {
namespace refdist {

//+PLUMEDOC FUNCTION DIFFERENCE
/*
Calculate the differences between two scalars or vectors

This action can be used to calculate the difference between two values.  For example, if we want to calculate the difference
between the distances between atoms 1 and 2 and 3 and 4 we can use an input like the one shown below:

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
diff: DIFFERENCE ARG=d1,d2
PRINT ARG=diff FILE=colvar
```

At first sight this action appears pointless as you can achieve the same result with the following input:

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
diff: CUSTOM ARG=d1,d2 FUNC=x-y PERIODIC=NO
PRINT ARG=diff FILE=colvar
```

However, this second example will not give an equivalent result to the following input:

```plumed
t1: TORSION ATOMS=1,2,3,4
t2: TORSION ATOMS=5,6,7,8
diff: DIFFERENCE ARG=t1,t2
PRINT ARG=diff FILE=colvar
```

as the [CUSTOM](CUSTOM.md) command cannot caclulate differences between periodic input variables.  In this case you have to
use DIFFERENCE to ensure that the periodic boundary conditions are applied when the final difference is calculated.

## Difference from reference value/s

The DIFFERENCE command is frequently used to determine the difference between the instantaneous value of a quantity and
a reference value as illustrated in the example below:

```plumed
ref_psi: CONSTANT VALUES=2.25029540
t: TORSION ATOMS=1,2,3,4
diff: DIFFERENCE ARG=t,ref_psi
PRINT ARG=diff FILE=colvar
```

You can also use this action to calculate the difference between the instaneous values of a vector of quantities and a
reference vector as illustrated in the example below:

```plumed
ref: CONSTANT VALUES=1.387980,1.126120,1.269380,1.321120,1.212420
d: DISTANCES ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8 ATOMS5=9,10
diff: DIFFERENCE ARG=d,ref
PRINT ARG=diff FILE=colvar
```

The output in this case is a five dimensional vector.  Each element of this vector contains the difference between the
corresponding elements of the two input vectors.

Notice that you __cannot__ use a mixture of scalars and vectors in the input for this action. Only two values can be input
to this action and these two input values __must__ have the same rank and the same number of elements.

*/
//+ENDPLUMEDOC

class Difference {
public:
  bool periodic;
  double max_minus_min;
  double inv_max_minus_min;
  static void registerKeywords(Keywords& keys);
  static void read( Difference& func, ActionWithArguments* action, function::FunctionOptions& options );
  static void calc( const Difference& func, bool noderiv, const View<const double,helpers::dynamic_extent>& args, function::FunctionOutput& funcout );
};


typedef function::FunctionShortcut<Difference> DifferenceShortcut;
PLUMED_REGISTER_ACTION(DifferenceShortcut,"DIFFERENCE")
typedef function::FunctionOfScalar<Difference> ScalarDifference;
PLUMED_REGISTER_ACTION(ScalarDifference,"DIFFERENCE_SCALAR")
typedef function::FunctionOfVector<Difference> VectorDifference;
PLUMED_REGISTER_ACTION(VectorDifference,"DIFFERENCE_VECTOR")

void Difference::registerKeywords(Keywords& keys) {
  keys.setValueDescription("scalar/vector","a function that measures the difference");
}

void Difference::read( Difference& func, ActionWithArguments* action, function::FunctionOptions& options ) {
  if( action->getNumberOfArguments()!=2 ) {
    action->error("should be two arguments to this action");
  }
  if( action->getPntrToArgument(0)->getRank()==action->getPntrToArgument(1)->getRank() ) {
    std::vector<std::size_t> shape( action->getPntrToArgument(0)->getShape() );
    for(unsigned i=0; i<shape.size(); ++i) {
      if( shape[i]!=action->getPntrToArgument(1)->getShape()[i] ) {
        action->error("shapes of input actions do not match");
      }
    }
  }

  func.periodic=false;
  std::string min0, max0;
  if( action->getPntrToArgument(0)->isPeriodic() ) {
    func.periodic=true;
    action->getPntrToArgument(0)->getDomain( min0, max0 );
    if( !action->getPntrToArgument(1)->isConstant() && !action->getPntrToArgument(1)->isPeriodic() ) {
      action->error("period for input variables " + action->getPntrToArgument(0)->getName() + " and " + action->getPntrToArgument(1)->getName() + " should be the same 0");
    }
    if( !action->getPntrToArgument(1)->isConstant() ) {
      std::string min1, max1;
      action->getPntrToArgument(1)->getDomain( min1, max1 );
      if( min0!=min0 || max0!=max1 ) {
        action->error("domain for input variables should be the same");
      }
    } else {
      action->getPntrToArgument(1)->setDomain( min0, max0 );
    }
  } else if( action->getPntrToArgument(1)->isPeriodic() ) {
    func.periodic=true;
    action->getPntrToArgument(1)->getDomain( min0, max0 );
    if( !action->getPntrToArgument(0)->isConstant() ) {
      action->error("period for input variables " + action->getPntrToArgument(0)->getName() + " and " + action->getPntrToArgument(1)->getName() + " should be the same 1");
    } else {
      action->getPntrToArgument(0)->setDomain( min0, max0 );
    }
  }
  if( func.periodic ) {
    double min, max;
    Tools::convert( min0, min );
    Tools::convert( max0, max );
    func.max_minus_min=max-min;
    func.inv_max_minus_min=1.0/func.max_minus_min;
  }
}

void Difference::calc( const Difference& func, bool noderiv, const View<const double,helpers::dynamic_extent>& args, function::FunctionOutput& funcout ) {
  plumed_dbg_assert( args.size()==2 );
  if( func.periodic ) {
    funcout.values[0] = func.max_minus_min*Tools::pbc( func.inv_max_minus_min*(args[0] - args[1]) );
  } else {
    funcout.values[0] = args[0] - args[1];
  }
  if( !noderiv ) {
    funcout.derivs[0][0] = 1.0;
    funcout.derivs[0][1]=-1;
  }
}

}
}


