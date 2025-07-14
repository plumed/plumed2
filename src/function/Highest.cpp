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
#include "core/ActionRegister.h"
#include "FunctionShortcut.h"
#include "FunctionOfScalar.h"
#include "FunctionOfVector.h"
#include "FunctionWithSingleArgument.h"
#include "FunctionSetup.h"

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION HIGHEST
/*
This function can be used to find the highest colvar by magnitude in a set.

This action allows you to find the highest of the input arguments.  As a first example of how it might be used consider the following input:

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
h1: HIGHEST ARG=d1,d2
PRINT ARG=h1 FILE=colvar
```

The value, `h1`, that is output to the file `colvar` here will be equal to `d1` if `d1>d2` and will be equal to `d2` if `d2>d1`.  In other words, if
all the arguments input to a HIGHEST action are scalars then the output value will be a scalar that is equal to the largest of the input arguments.
Notice that you can also use this command with more than two arguments as illustrated below:

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
d3: DISTANCE ATOMS=5,6
d4: DISTANCE ATOMS=7,8
h1: HIGHEST ARG=d1,d2,d3,d4
PRINT ARG=h1 FILE=colvar
```

## Using a single vector as input

Instead of inputting multiple scalars you can input a single vector to this action instead as is illustrated below:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4
h1: HIGHEST ARG=d
PRINT ARG=h1 FILE=colvar
```

The output from this action is a single scalar once again.  This single scalar is equal to the largest element of the input vector.

## Using multiple vectors in input

If you input multiple vectors with the same numbers of elements to this action, as shown below, the output will be a vector.


```plumed
d1: DISTANCE ATOMS1=1,2 ATOMS2=3,4
d2: DISTANCE ATOMS1=5,6 ATOMS2=7,8
d3: DISTANCE ATOMS1=9,10 ATOMS2=11,12
h2: HIGHEST ARG=d1,d2,d3
PRINT ARG=h2 FILE=colvar
```

The elements of the output vector here are determined by doing an elementwise comparison of the elements in the input vectors.  In the above
input the first element of `h2` is equal to the distance between atoms 1 and 2 if this distance is larger than the distances between atoms 5 and 6 and the distance between atoms 9 and 10.
By the same token the second element of `h2` is equal to the distance between atoms 3 and 4 if this is larger than the distance between atoms 7 and 8 and the distance between atoms 11 and 12.
In other words, if the elements of the $j$th input vector are given by $v_i^{(j)}$ then the elements of the output vector, $h_i$ are given by:

$$
h_i = \max_j v_i^{(j)}
$$

Notice that you can also use a combination of scalars and vectors in the input to this action as shown below:

```plumed
c: CONSTANT VALUE=0.05
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4
h: HIGHEST ARG=d,c
PRINT ARG=h FILE=colvar
```

For the input above the HIGHEST action outputs a vector with two elements.  The elements of this vector are equal to the distances between the pairs
of atoms that are specified in the DISTANCE command as long as those distances are greater than 0.05 nm.  If either of the two input distances is less
than 0.05 nm then the corresponding value in the vector `h` is set equal to 0.05 nm.

## The MASK keyword

You should only use the MASK keyword for this action if you are using multiple arguments and at least one vector in the input to this action.
The following input illustrates how this keyword works:

```plumed
m: CONSTANT VALUES=1,0
d1: DISTANCE ATOMS1=1,2 ATOMS2=3,4
d2: DISTANCE ATOMS1=5,6 ATOMS2=7,8
d3: DISTANCE ATOMS1=9,10 ATOMS2=11,12
h2: HIGHEST ARG=d1,d2,d3 MASK=m
```

By using the MASK keyword here you ensure that only the first element of the output vector for `h2` is calculated. This element of the output vector
is evaluated as the corresponding element of the vector that is passed in the input to the MASK keyword is one.  By the same logic the second element
of the output vector is not computed because the corresponding element in the mask vector is 0.

The input above illustrates the key idea for the MASK keyword and is probably not very useful. To see an example where using a MASK for this type of action
is useful you should look at the expanded shortcuts in the documentation for [PARABETARMSD](PARABETARMSD.md).

*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION LOWEST
/*
This function can be used to find the lowest colvar by magnitude in a set.

This action allows you to find the lowest of the input arguments.  As a first example of how it might be used consider the following input:

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
l1: LOWEST ARG=d1,d2
PRINT ARG=l1 FILE=colvar
```

The value, `l1`, that is output to the file `colvar` here will be equal to `d1` if `d1<d2` and will be equal to `d2` if `d2<d1`.  In other words, if
all the arguments input to a LOWEST action are scalars then the output value will be a scalar that is equal to the smallest of the input arguments.
Notice that you can also use this command with more than two arguments as illustrated below:

```plumed
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
d3: DISTANCE ATOMS=5,6
d4: DISTANCE ATOMS=7,8
l1: LOWEST ARG=d1,d2,d3,d4
PRINT ARG=l1 FILE=colvar
```

## Using a single vector as input

Instead of inputting multiple scalars you can input a single vector to this action instead as is illustrated below:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4
l1: LOWEST ARG=d
PRINT ARG=l1 FILE=colvar
```

The output from this action is a single scalar once again.  This single scalar is equal to the smallest element of the input vector.

## Using multiple vectors in input

If you input multiple vectors with the same numbers of elements to this action, as shown below, the output will be a vector.


```plumed
d1: DISTANCE ATOMS1=1,2 ATOMS2=3,4
d2: DISTANCE ATOMS1=5,6 ATOMS2=7,8
d3: DISTANCE ATOMS1=9,10 ATOMS2=11,12
l2: LOWEST ARG=d1,d2,d3
PRINT ARG=l2 FILE=colvar
```

The elements of the output vector here are determined by doing an elementwise comparison of the elements in the input vectors.  In the above
input the first element of `h2` is equal to the distance between atoms 1 and 2 if this distance is smaller than the distances between atoms 5 and 6 and the distance between atoms 9 and 10.
By the same token the second element of `h2` is equal to the distance between atoms 3 and 4 if this is smaller than the distance between atoms 7 and 8 and the distance between atoms 11 and 12.
In other words, if the elements of the $j$th input vector are given by $v_i^{(j)}$ then the elements of the output vector, $h_i$ are given by:

$$
h_i = \min_j v_i^{(j)}
$$

Notice that you can also use a combination of scalars and vectors in the input to this action as shown below:

```plumed
c: CONSTANT VALUE=0.5
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4
h: LOWEST ARG=d,c
PRINT ARG=h FILE=colvar
```

For the input above the LOWEST action outputs a vector with two elements.  The elements of this vector are equal to the distances between the pairs
of atoms that are specified in the DISTANCE command as long as those distances are less than 0.5 nm.  If either of the two input distances is more
than 0.5 nm then the corresponding value in the vector `h` is set equal to 0.5 nm.

## The MASK keyword

You should only use the MASK keyword for this action if you are using multiple arguments and at least one vector in the input to this action.
The following input illustrates how this keyword works:

```plumed
m: CONSTANT VALUES=1,0
d1: DISTANCE ATOMS1=1,2 ATOMS2=3,4
d2: DISTANCE ATOMS1=5,6 ATOMS2=7,8
d3: DISTANCE ATOMS1=9,10 ATOMS2=11,12
h2: LOWEST ARG=d1,d2,d3 MASK=m
```

By using the MASK keyword here you ensure that only the first element of the output vector for `h2` is calculated. This element of the output vector
is evaluated as the corresponding element of the vector that is passed in the input to the MASK keyword is one.  By the same logic the second element
of the output vector is not computed because the corresponding element in the mask vector is 0.

The input above illustrates the key idea for the MASK keyword and is probably not very useful. To see an example where using a MASK for this type of action
is useful you should look at the expanded shortcuts in the documentation for [PARABETARMSD](PARABETARMSD.md).

*/
//+ENDPLUMEDOC

class Highest {
public:
  bool min;
  static void registerKeywords( Keywords& keys );
  static void read( Highest& func, ActionWithArguments* action, FunctionOptions& options );
  static void calc( const Highest& func, bool noderiv, const View<const double>& args, FunctionOutput& funcout );
};

typedef FunctionShortcut<Highest> HighestShortcut;
PLUMED_REGISTER_ACTION(HighestShortcut,"HIGHEST")
PLUMED_REGISTER_ACTION(HighestShortcut,"LOWEST")
typedef FunctionWithSingleArgument<Highest> OneargHighest;
PLUMED_REGISTER_ACTION(OneargHighest,"HIGHEST_ONEARG")
PLUMED_REGISTER_ACTION(OneargHighest,"LOWEST_ONEARG")
typedef FunctionOfScalar<Highest> ScalarHighest;
PLUMED_REGISTER_ACTION(ScalarHighest,"HIGHEST_SCALAR")
PLUMED_REGISTER_ACTION(ScalarHighest,"LOWEST_SCALAR")
typedef FunctionOfVector<Highest> VectorHighest;
PLUMED_REGISTER_ACTION(VectorHighest,"HIGHEST_VECTOR")
PLUMED_REGISTER_ACTION(VectorHighest,"LOWEST_VECTOR")

void Highest::registerKeywords( Keywords& keys ) {
  if( keys.getDisplayName().find("LOWEST") ) {
    keys.setValueDescription("scalar","the lowest of the input values");
  } else {
    keys.setValueDescription("scalar","the highest of the input values");
  }
  keys.add("hidden","MASKED_INPUT_ALLOWED","turns on that you are allowed to use masked inputs");
}

void Highest::read( Highest& func, ActionWithArguments* action, FunctionOptions& options ) {
  func.min=action->getName().find("LOWEST")!=std::string::npos;
  if( !func.min ) {
    plumed_assert( action->getName().find("HIGHEST")!=std::string::npos );
  }
  for(unsigned i=0; i<action->getNumberOfArguments(); ++i) {
    if( action->getPntrToArgument(i)->isPeriodic() ) {
      action->error("Cannot sort periodic values (check argument "+ action->getPntrToArgument(i)->getName() +")");
    }
  }
}

void Highest::calc( const Highest& func, bool noderiv, const View<const double>& args, FunctionOutput& funcout ) {
  if( !noderiv ) {
    for(unsigned i=0; i<args.size(); ++i) {
      funcout.derivs[0][i] = 0;
    }
  }
  if( func.min ) {
    funcout.values[0] = *std::min_element(args.begin(), args.end());
    if( !noderiv ) {
      funcout.derivs[0][std::min_element(args.begin(), args.end()) - args.begin()] = 1;
    }
  } else {
    funcout.values[0] = *std::max_element(args.begin(), args.end());
    if( !noderiv ) {
      funcout.derivs[0][std::max_element(args.begin(), args.end()) - args.begin()] = 1;
    }
  }
}

}
}


