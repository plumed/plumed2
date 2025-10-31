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

## Using multiple vectors in input

If you input multiple vectors with the same numbers of elements to this action, as shown below, the output will be a set of vectors.

```plumed
d1: DISTANCE ATOMS1=1,2 ATOMS2=3,4
d2: DISTANCE ATOMS1=5,6 ATOMS2=7,8
d3: DISTANCE ATOMS1=9,10 ATOMS2=11,12
sort: SORT ARG=d1,d2,d3
PRINT ARG=sort.1,sort.2,sort.3 FILE=colvar
```

The elements of the three output vector here are determined by doing an elementwise comparison of the elements in the input vectors.  In the above
input the first element of `sort.1` is equal to the distance between atoms 1 and 2 if this distance is smaller than the distances between atoms 5 and 6 and the distance between atoms 9 and 10.
By the same token the second element of `sort.1` is equal to the distance between atoms 3 and 4 if this is smaller than the distance between atoms 7 and 8 and the distance between atoms 11 and 12.
The elements of `sort.2` are the second largest of these distances, while the elements of `sort.3` are three largest distances

Notice that you can also use a combination of scalars and vectors in the input to this action as shown below:

```plumed
c: CONSTANT VALUE=0.5
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4
sort: SORT ARG=d,c
PRINT ARG=sort.1,sort.2 FILE=colvar
```

For the input above the SORT action outputs two vectors with two elements.  The elements of `sort.1` are equal to the distances between the pairs
of atoms that are specified in the DISTANCE command as long as those distances are less than 0.5 nm.  If either of the two input distances is more
than 0.5 nm then the corresponding value in the vector `h` is set equal to 0.5 nm.  Similarly, the elements of `sort.2` are equal to the distances between
of the atoms that are  specified in the DISTANCE command as long as those distances are more than 0.5 nm.

## The MASK keyword

You should only use the MASK keyword for this action if you are using multiple arguments and at least one vector in the input to this action.
The following input illustrates how this keyword works:

```plumed
m: CONSTANT VALUES=1,0
d1: DISTANCE ATOMS1=1,2 ATOMS2=3,4
d2: DISTANCE ATOMS1=5,6 ATOMS2=7,8
d3: DISTANCE ATOMS1=9,10 ATOMS2=11,12
sort: SORT ARG=d1,d2,d3 MASK=m
```

By using the MASK keyword here you ensure that only the first elements of the output vectors for `sort` is calculated. These elements of the output vectors
are evaluated as the corresponding element of the vector that is passed in the input to the MASK keyword is one.  By the same logic the second element
of the output vectors are not computed because the corresponding element in the mask vector is 0.

The input above illustrates the key idea for the MASK keyword and is probably not very useful. To see an example where using a MASK for this type of action
is useful you should look at the expanded shortcuts in the documentation for [PARABETARMSD](PARABETARMSD.md).

*/
//+ENDPLUMEDOC

class Sort {
public:
  static void registerKeywords(Keywords& keys);
  static void read( Sort& func, ActionWithArguments* action, FunctionOptions& options );
  static void calc( const Sort& func, bool noderiv, View<const double> args, FunctionOutput& funcout );
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
  if( action->getNumberOfArguments()==1 ) {
    action->warning("if there is only one argument to sort is this function really needed?");
  }

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

void Sort::calc( const Sort& func, bool noderiv, const View<const double> args, FunctionOutput& funcout ) {
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
  for(unsigned i=0; i<funcout.values.size(); ++i) {
    funcout.values[i] = data[i].first;
    if( !noderiv ) {
      funcout.derivs[i][data[i].second] = 1;
    }
  }
}

}
}


