/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2023 The plumed team
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
#include "core/ActionPilot.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

//+PLUMEDOC ANALYSIS COLLECT
/*
Collect data from the trajectory for later analysis

The way this command can be used is illustrated by the following example:

```plumed
d: DISTANCE ATOMS=1,2
c: COLLECT ARG=d STRIDE=1
DUMPVECTOR ARG=c FILE=timeseries
```

The COLLECT command in the input above stores the time series of distances over the whole
trajectory. The STRIDE keyword controls how frequently data is stored and the [DUMPVECTOR](DUMPVECTOR.md)
command then outputs the full time series of `d` values at the end of the calculation.

The above example is not particularly useful as you can achieve the same result by using simply
a DISTANCE command and the [PRINT](PRINT.md) command.  The COLLECT command is useful if you want
use the whole trajectory to perform a analysis such as  dimensionality reduction (see [dimred](module_dimred.md)).
The shortcut command [COLLECT_FRAMES](COLLECT_FRAMES.md) uses this action heavily and allows one to easily
deal with the fact that COLLECT can only collect the time series for one PLUMED Value at a time.

## Collecting a part of the trajectory

You can use the CLEAR keyword to collect a subset of the trajectory as illustrated below:

```plumed
d: DISTANCE ATOMS=1,2
c: COLLECT ARG=d STRIDE=1 CLEAR=1000
DUMPVECTOR ARG=c FILE=timeseries STRIDE=1000
```

This outputs files contain from this input contain 1000-frame blocks of the trajectory for the distance between atom 1 and 2.
This type of input might proove useful if you wish to perform separate analysis on different parts of the trajectory for
later comparison.

An alternative way to collect part of the trajectry is to use the UPDATE_FROM and UPDATE_UNTIL arguments as illustrated in the
example shown below:

```plumed
d: DISTANCE ATOMS=1,2
c: COLLECT ARG=d STRIDE=1 UPDATE_FROM=100 UPDATE_UNTIL=200
DUMPVECTOR ARG=c FILE=timeseries
```

In this input the COLLECT action only starts collecting the data on step 100 so the first 100 values of the trajectory are discarded.
The COLLECT action then stops collecting data on step 200. These 200 stored distances are then output at the end of the simulation
even if the trajectory contains more than 200 frames.

## Collecting vectors

You can use the collect command even if the input value has a rank that is greater than zero.  For example, the following
input collects vectors of distances from the trajectory:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
c: COLLECT ARG=d TYPE=vector STRIDE=1 CLEAR=100
DUMPVECTOR ARG=c FILE=timeseries STRIDE=100
```

Notice that if the input to the collect action is a value with rank>0 you __must__ use the TYPE keyword to specify whether
the output Value is a vector or matrix.  In the above input we are storing a vector so the DUMPVECTOR command outputs
a list of 400 distances - 4 distances for each frame.  The assumption in the input above is that the four distances that have
been computed by the DISTANCES command are indistinguishable.

By using the following input we can ensure the the four distances are treated as distinguishable when we do any analysis:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8
c: COLLECT ARG=d TYPE=matrix STRIDE=1 CLEAR=100
DUMPVECTOR ARG=c FILE=timeseries STRIDE=100
```

The value `c` is now a $100\times 4$ matrix.  Furthermore, when we use DUMPVECTOR to output the output file will contain
four columns of data.

## Collecting matrices

You can also use this action to collect a matrix as illustrated by the following example:

```plumed
d: DISTANCE_MATRIX GROUPA=1,2 GROUPB=3,4,5
c: COLLECT ARG=d TYPE=vector STRIDE=1 CLEAR=100
DUMPVECTOR ARG=c FILE=timeseries STRIDE=100
```

The value `c` here is a vector with 600 elements as the input matrix is converted to a vector. These vectors are then stored in
one contiguous object.

If by contrast we use `TYPE=matrix` as shown below:

```plumed
d: DISTANCE_MATRIX GROUPA=1,2 GROUPB=3,4,5
c: COLLECT ARG=d TYPE=matrix STRIDE=1 CLEAR=100
DUMPVECTOR ARG=c FILE=timeseries STRIDE=100
```

A $100 \times 6$ matrix is stored.  Each row of this matrix contains a vectorized version of the input matrix.
There is currently no way to store the collected data in a way that recognises that each of the input PLUMED Value
was a matrix. You also cannot use this action to store functions evaluated on a grid.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace generic {

class Collect :
  public ActionWithValue,
  public ActionWithArguments,
  public ActionPilot {
private:
  bool usefirstconf;
  unsigned clearstride;
public:
  static void registerKeywords( Keywords& keys );
  Collect( const ActionOptions& );
  unsigned getNumberOfDerivatives() override;
  bool calculateOnUpdate() override {
    return false;
  }
  bool calculateConstantValues( const bool& have_atoms ) override {
    return false;
  }
  void calculate() override {}
  void apply() override {}
  void update() override ;
};

PLUMED_REGISTER_ACTION(Collect,"COLLECT")

void Collect::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.addInputKeyword("compulsory","ARG","scalar/vector/matrix","the label of the value whose time series is being stored for later analysis");
  keys.add("compulsory","STRIDE","1","the frequency with which the data should be collected and added to the quantity being averaged");
  keys.add("compulsory","CLEAR","0","the frequency with which to clear all the accumulated data.  The default value "
           "of 0 implies that all the data will be used and that the grid will never be cleared");
  keys.add("compulsory","TYPE","auto","required if you are collecting an object with rank>0. Should be vector/matrix and determines how data is stored.  If rank==0 then data has to be stored as a vector");
  keys.setValueDescription("vector/matrix","the time series for the input quantity");
}

Collect::Collect( const ActionOptions& ao ):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  ActionPilot(ao),
  usefirstconf(false) {
  if( getNumberOfArguments()!=1 ) {
    error("there should only be one argument to this action");
  }
  if( getPntrToArgument(0)->getRank()>0 && getPntrToArgument(0)->hasDerivatives() ) {
    error("input to the collect argument cannot be a grid");
  }

  std::string type;
  parse("TYPE",type);
  if( getPntrToArgument(0)->getNumberOfValues()==1 && (type=="auto" || type=="vector") ) {
    type="vector";
  } else if( getPntrToArgument(0)->getNumberOfValues()==1 && type=="matrix" ) {
    error("invalid type specified. Cannot construct a matrix by collecting scalars");
  } else if(  getPntrToArgument(0)->getNumberOfValues()!=1 && type=="auto" ) {
    error("missing TYPE keyword.  TYPE should specify whether data is to be stored as a vector or a matrix");
  } else if( type!="vector" && type!="matrix" ) {
    error("invalid TYPE specified. Should be matrix/scalar found " + type);
  }

  if( type=="vector" ) {
    log.printf("  adding %ld elements to stored vector each time we collect\n", getPntrToArgument(0)->getNumberOfValues() );
  } else {
    log.printf("  constructing matrix with rows of length %ld from input data\n", getPntrToArgument(0)->getNumberOfValues() );
  }

  parse("CLEAR",clearstride);
  std::size_t nvals=0;
  if( clearstride==getStride() ) {
    nvals=1;
    usefirstconf=(getStride()==0);
  } else if( clearstride>0 ) {
    if( clearstride%getStride()!=0 ) {
      error("CLEAR parameter must be a multiple of STRIDE");
    }
    log.printf("  clearing collected data every %u steps \n",clearstride);
    nvals=(clearstride/getStride());
  }

  std::vector<std::size_t> shape(1);
  shape[0]=nvals;
  if( type=="matrix" ) {
    shape.resize(2);
    shape[1] = getPntrToArgument(0)->getNumberOfValues();
  }
  if( type=="vector" ) {
    shape[0] = nvals*getPntrToArgument(0)->getNumberOfValues();
  }
  addValue( shape );
  if( shape.size()==2 ) {
    getPntrToComponent(0)->reshapeMatrixStore( shape[1] );
  }
  if( getPntrToArgument(0)->isPeriodic() ) {
    std::string min, max;
    getPntrToArgument(0)->getDomain( min, max );
    setPeriodic( min, max );
  } else {
    setNotPeriodic();
  }
}

unsigned Collect::getNumberOfDerivatives() {
  return 0;
}

void Collect::update() {
  if( getStep()==0 || (!onStep() && !usefirstconf) ) {
    return ;
  }
  usefirstconf=false;

  Value* myin=getPntrToArgument(0);
  Value* myout=getPntrToComponent(0);
  unsigned nargs=myin->getNumberOfValues();
  if( clearstride==getStride() ) {
    for(unsigned i=0; i<nargs; ++i) {
      myout->set( i, myin->get(i) );
    }
  } else if( clearstride>0 ) {
    unsigned step = getStep() - clearstride*std::floor( getStep() / clearstride );
    if( getStep()%clearstride==0 ) {
      step = step + clearstride;
    }
    unsigned base = (step/getStride()-1)*nargs;
    for(unsigned i=0; i<nargs; ++i) {
      myout->set( base+i, myin->get(i) );
    }
  } else {
    for(unsigned i=0; i<nargs; ++i) {
      myout->push_back( myin->get(i) );
    }
    if( myout->getRank()==2 ) {
      myout->reshapeMatrixStore( nargs );
    }
  }
}

}
}
