/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "core/ActionRegister.h"
#include "FunctionSetup.h"

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION SORT
/*
This function can be used to sort colvars according to their magnitudes.

This action sorts the input arguments according to their magnitudes. It will output the same number of values as it has arguments.
The lowest argument input will be output as a value labelled _label_.1, the second lowest input argument will be output as a value labelled _label_.2 and so on.
Thus, for example, the following input can be used to print the distance of the closest and of the farthest atoms to atom 1, chosen among atoms from 2 to 5

```plumed
d12:  DISTANCE ATOMS=1,2
d13:  DISTANCE ATOMS=1,3
d14:  DISTANCE ATOMS=1,4
d15:  DISTANCE ATOMS=1,5
sort: SORT ARG=d12,d13,d14,d15
PRINT ARG=sort.1,sort.4
```

Notice that you can also achieve the same result using the following input:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=1,3 ATOMS3=1,4 ATOMS4=1,5
sort: SORT ARG=d
PRINT ARG=sort.1,sort.4
```

In this second case the four distances are passed to the SORT action as a vector.  The SORT action then outputs 4 components - the same number of components as there
are elements in the input vector - that contain the elements of the input vector in order of increasing magnitude.

These examples are representative the only two ways you can use this action.  In input it can accept either a list of scalars or a single vector.
It does not accept matrices or a list of vectors in input.

*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION SORT_SCALAR
/*
Sort the input scalars in a vector according to their magnitudes

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION SORT_VECTOR
/*
Sort the elements in a vector according to their magnitudes

\par Examples

*/
//+ENDPLUMEDOC

class Sort {
public:
  static void registerKeywords(Keywords& keys);
  static void read( Sort& func, ActionWithArguments* action, FunctionOptions& options );
  static void calc( const Sort& func, bool noderiv, const View<const double,helpers::dynamic_extent>& args, FunctionOutput& funcout );
  Sort& operator=(const Sort& m) {
    return *this;
  }
};

typedef FunctionShortcut<Sort> SortShortcut;
PLUMED_REGISTER_ACTION(SortShortcut,"SORT")
typedef FunctionOfScalar<Sort> ScalarSort;
PLUMED_REGISTER_ACTION(ScalarSort,"SORT_SCALAR")
typedef FunctionOfVector<Sort> VectorSort;
PLUMED_REGISTER_ACTION(VectorSort,"SORT_VECTOR")

void Sort::registerKeywords(Keywords& keys) {
  keys.setValueDescription("vector","sorted");
  keys.setComponentsIntroduction("The names of the components in this action will be customized in accordance with the contents of the input file. "
                                 "The largest value is called label.1th, the second largest label.2th, the third label.3th and so on");
}


void Sort::read( Sort& func, ActionWithArguments* action, FunctionOptions& options ) {
  if( action->getNumberOfArguments()==1 ) action->warning("if there is only one argument to sort is this function really needed?");

  options.derivativeZeroIfValueIsZero=true;
  for(unsigned i=0; i<action->getNumberOfArguments(); ++i) {
    if((action->getPntrToArgument(i))->isPeriodic()) {
      action->error("Cannot sort periodic values (check argument "+ (action->getPntrToArgument(i))->getName() +")");
    }
    if(!(action->getPntrToArgument(i))->isDerivativeZeroWhenValueIsZero() ) {
      options.derivativeZeroIfValueIsZero=false;
    }
  }
}

void Sort::calc( const Sort& func, bool noderiv, const View<const double,helpers::dynamic_extent>& args, FunctionOutput& funcout ) {
  if( !noderiv ) {
    for(unsigned i=0; i<args.size(); ++i) {
      for(unsigned j=0; j<args.size(); ++j) {
        funcout.derivs[i][j] = 0;
      }
    }
  }
  std::vector<std::pair<double,int> > data(args.size());
  for(unsigned i=0; i<args.size(); ++i) {
    data[i].first=args[i];
// In this manner I remember from which argument the component depends:
    data[i].second=i;
  }
// STL sort sorts based on first element (value) then second (index)
  std::sort(data.begin(),data.end());
  for(int i=0; i<funcout.values.size(); ++i) {
    funcout.values[i] = data[i].first;
    if( !noderiv ) {
      funcout.derivs[i][data[i].second] = 1;
    }
  }
}

}
}


