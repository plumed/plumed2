/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"
#include "tools/Random.h"

namespace PLMD {
namespace generic {

//+PLUMEDOC DIMRED CREATE_MASK
/*
Create a mask vector to use for landmark selection

This action should be used in conjuction with the [SELECT_WITH_MASK](SELECT_WITH_MASK.md) action. As is explained in the documentation
for [SELECT_WITH_MASK](SELECT_WITH_MASK.md), [SELECT_WITH_MASK](SELECT_WITH_MASK.md)
allows you to output a scalar vector or matrix that contains a subset of the elements in the input vector or matrix.  The following example
shows how this action works in practice.

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6
m: CONSTANT VALUES=1,0,1
v: SELECT_WITH_MASK ARG=d MASK=m
```

In the input above, The value, `m`, that is passed to the keyword MASK here is a vector with the same length as `d`.
Elements of `d` that whose corresponding elements in `m` are zero are copied to the output value `v`.
When elements of `m` are non-zero the corresponding elements in `d` are not transferred to the output
value - they are masked.  Consequently, although the input, `d`, to the select SELECT_WITH_MASK action is a vector
with 3 elements the output, `v`, is a scalar.

In the vast majority of cases you can use the [CONSTANT](CONSTANT.md) action to create the values that are passed to the
SELECT_WITH_MASK action using the `MASK` keyword as the size of the input vector is known when you write the input file.
The one exception to this is if you are using a [COLLECT](COLLECT.md) action like the one shown below:

```plumed
d: DISTANCE ATOMS=1,2
c: COLLECT ARG=d STRIDE=1
```

The vector `c` here will have as many frames as there are frames in the trajectory that is being generated or analysed.
Normally, with an input like this a user will want to analyse all the data in the whole trajectory at once without
specifying the trajectories length to PLUMED.  action should was desgined to create inputs to the MASK keyword for
[SELECT_WITH_MASK](SELECT_WITH_MASK.md) action in this particular case.  Basically it creates a vector of ones and
zeros that has the same length as the input vector.

There are three modes for creating these vectors and ones and zeros.  This first one creates a mask vector in which
every element is equal to zero.  The output from the SELECT_WITH_MASK command, `v`, is thus identical to the output from
the collect command, `d`.

```plumed
d: DISTANCE ATOMS=1,2
c: COLLECT ARG=d STRIDE=1
m: CREATE_MASK TYPE=nomask ARG=c
v: SELECT_WITH_MASK ARG=c MASK=m
```

The second mode, which is illustrated below, sets N randomly-chosen elements of the mask equal to zero.

```plumed
d: DISTANCE ATOMS=1,2
c: COLLECT ARG=d STRIDE=1
m: CREATE_MASK TYPE=random NZEROS=20 SEED=23623 ARG=c
v: SELECT_WITH_MASK ARG=c MASK=m
```

The vector `v` that is output by the SELECT_WITH_MASK command contains 20 randomly selected elements from the input vector `c`.

The third and final mode, which is illustrated below, selects a set of N evenly spaced points from the input vector.

```plumed
d: DISTANCE ATOMS=1,2
c: COLLECT ARG=d STRIDE=1
m: CREATE_MASK TYPE=stride NZEROS=20 ARG=c
v: SELECT_WITH_MASK ARG=c MASK=m
```

The vector `v` that is output by the SELECT_WITH_MASK command contains 20 elements from the input vector `c`.  If `c` has
200 elements then `v` will contain every 20th element of the input vector `c`.

For more examples that demonstrate how this action is used look at the actions in the [landmarks](module_landmarks.md) module.


*/
//+ENDPLUMEDOC

class CreateMask :
  public ActionWithValue,
  public ActionWithArguments {
private:
  Random r;
  unsigned nzeros;
  enum {nomask,stride,random} type;
public:
  static void registerKeywords( Keywords& keys );
  CreateMask( const ActionOptions& );
  unsigned getNumberOfDerivatives() override {
    return 0;
  }
  void prepare() override ;
  void calculate() override ;
  void apply() override {}
};

PLUMED_REGISTER_ACTION(CreateMask,"CREATE_MASK")

void CreateMask::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.addInputKeyword("compulsory","ARG","vector","the label of the vector that you would like to construct a mask for");
  keys.add("compulsory","TYPE","the way the zeros are supposed to be set");
  keys.add("compulsory","NZEROS","the number of zeros that you want to put in the mask");
  keys.add("optional","SEED","the seed to use for the random number generator");
  keys.setValueDescription("vector","a vector of zeros and ones that is used that can be used to mask some of the elements in a time series");
}


CreateMask::CreateMask( const ActionOptions& ao ) :
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  nzeros(0) {
  if( getNumberOfArguments()!=1 ) {
    error("should only be one argument to this action");
  }
  if( getPntrToArgument(0)->getRank()!=1 ) {
    error("argument should be a vector");
  }
  std::string stype;
  parse("TYPE",stype);
  if( stype!="nomask" ) {
    parse("NZEROS",nzeros);
  }

  if( stype=="nomask" ) {
    type=nomask;
    log.printf("  setting all points in output mask to zero \n");
  } else if( stype=="stride" ) {
    type=stride;
    log.printf("  setting every %d equally spaced points in output mask to zero \n", nzeros );
  } else if( stype=="random" ) {
    unsigned seed=230623;
    parse("SEED",seed);
    r.setSeed(-seed);
    type=random;
    log.printf("  choosing %d points to set to non-zero in mask in accordance with input weights \n", nzeros );
  } else {
    error( stype + " is not a valid way input for TYPE");
  }
  std::vector<std::size_t> shape(1);
  shape[0] = getPntrToArgument(0)->getShape()[0];
  addValue( shape );
  setNotPeriodic();
  for(unsigned i=0; i<shape[0]; ++i) {
    getPntrToComponent(0)->set( i, 1.0 );
  }
  // We run calculate here to make a selection so that we can parse the input.
  // This initial selection is never used
  if( shape[0]>0 ) {
    calculate();
  }
}

void CreateMask::prepare() {
  Value* out=getPntrToComponent(0);
  Value* arg=getPntrToArgument(0);
  if( out->getShape()[0]!=arg->getShape()[0] ) {
    std::vector<std::size_t> shape(1);
    shape[0] = arg->getShape()[0];
    out->setShape( shape );
  }
  if( type!=nomask ) {
    for(unsigned i=nzeros; i<out->getShape()[0]; ++i) {
      out->set( i, 1 );
    }
  }
}

void CreateMask::calculate() {
  Value* out=getPntrToComponent(0);
  Value* arg=getPntrToArgument(0);
  unsigned ns = arg->getShape()[0];
  for(unsigned i=0; i<ns; ++i) {
    out->set( i, 1.0 );
  }

  if( type==stride ) {
    std::size_t ss = int( std::floor( ns / nzeros ) );
    for(unsigned i=0; i<nzeros; ++i) {
      out->set( i*ss, 0.0 );
    }
  } else if( type==random ) {
    for(unsigned i=0; i<nzeros; ++i ) {
      double totweights = 0;
      for(unsigned j=0; j<ns; ++j) {
        if( out->get(j)>0 ) {
          totweights += arg->get(j);
        }
      }
      double rr = r.U01()*totweights;
      double accum=0;
      for(unsigned j=0; j<ns; ++j) {
        if( out->get(j)>0 ) {
          accum += arg->get(j);
        }
        if( accum<rr ) {
          continue;
        }
        out->set( j, 0 );
        break;
      }
    }
  } else if( type==nomask ) {
    for(unsigned i=0; i<ns; ++i) {
      out->set( i, 0.0 );
    }
  } else {
    error("invalid mask creation type");
  }
}

}
}
