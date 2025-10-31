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

//+PLUMEDOC GRIDCALC ACCUMULATE
/*
Sum the elements of this value over the course of the trajectory

This action is used to sum the outputs from another action over the course of the trajectory.  This is useful
if you want to calculate the average value that a CV took over the course of a simulation.  As an example, the following
input can be used to calculate the average distance between atom 1 and atom 2.

```plumed
c: CONSTANT VALUE=1
d: DISTANCE ATOMS=1,2
# This adds together the value of the distance on every step
s: ACCUMULATE ARG=d STRIDE=1
# This adds one every time we add a new distance to the value s
n: ACCUMULATE ARG=c STRIDE=1
# This is thus the average distance
a: CUSTOM ARG=s,n FUNC=x/y PERIODIC=NO
# This prints out the average over the whole trajectory (STRIDE=0 means print at end only)
PRINT ARG=a FILE=average.dat STRIDE=0
```

You can use this action for block averaging by using the `CLEAR` keyword as shown below:

```plumed
c: CONSTANT VALUE=1
d: DISTANCE ATOMS=1,2
# This adds together the value of the distance on every step
s: ACCUMULATE ARG=d STRIDE=1 CLEAR=1000
# This adds one every time we add a new distance to the value s
n: ACCUMULATE ARG=c STRIDE=1 CLEAR=1000
# This is thus the average distance
a: CUSTOM ARG=s,n FUNC=x/y PERIODIC=NO
# This prints out the average over the whole trajectory (STRIDE=0 means print at end only)
PRINT ARG=a FILE=average.dat STRIDE=1000
```

The instructions `CLEAR=1000` in the above input tells PLUMED to set the values `s` and `n` back to
zero after 1000 new steps have been performed. The PRINT action will thus print a block average that
is taken from the first 1000 steps of the trajectory, a second block average from the second 1000 steps
of the trajectory and so on.  Notice that you can achieve a similar effect using UPDATE_FROM and UPDATE_UNTIL as
shown below:

```plumed
c: CONSTANT VALUE=1
d: DISTANCE ATOMS=1,2
s1: ACCUMULATE ARG=d STRIDE=1 UPDATE_UNTIL=1000
n1: ACCUMULATE ARG=c STRIDE=1 UPDATE_UNTIL=1000
a1: CUSTOM ARG=s1,n1 FUNC=x/y PERIODIC=NO
s2: ACCUMULATE ARG=d STRIDE=1 UPDATE_FROM=1000
n2: ACCUMULATE ARG=c STRIDE=1 UPDATE_FROM=1000
a2: CUSTOM ARG=s2,n2 FUNC=x/y PERIODIC=NO
diff: CUSTOM ARG=a1,a2 FUNC=x-y PERIODIC=NO
PRINT ARG=a1,a2,diff FILE=colver STRIDE=2000
```

This output calculates the average distance for the first 1000 frames of the trajectory and for the second 1000 frames of the trajectory.
These two averages are then output to a file called colvar as well as the difference between them.

## Estimating histograms

We can also use this action to construct histograms. The example below shows how you can estimate the
distribution of distances between atoms 1 and 2 that are sampled over the course of the trajectory.

```plumed
c: CONSTANT VALUE=1
d: DISTANCE ATOMS=1,2
# Construct the instantaneous histogram from the instantaneous value of the distance
kde: KDE ARG=d BANDWIDTH=0.05 GRID_MIN=0 GRID_MAX=5 GRID_BIN=250
# Now add together all the instantaneous histograms
hist: ACCUMULATE ARG=kde STRIDE=1
# And normalise the histogram
n: ACCUMULATE ARG=c STRIDE=1
a: CUSTOM ARG=hist,n FUNC=x/y PERIODIC=NO
# And print out the final histogram
DUMPGRID ARG=a FILE=histo.grid
```

At first glance the fact that we use a [KDE](KDE.md) action to construct an instaneous histogram from a single
distance may appear odd.  The reason for doing this, however, is to introduce a clear distinction between
the syntaxes that are used for spatial and temporal averaging.  To see what I mean consider the following input:


```plumed
c: CONSTANT VALUE=5
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8 ATOMS5=9,10
# Construct the instantaneous histogram from the instantaneous value of the distance
kde: KDE ARG=d BANDWIDTH=0.05 GRID_MIN=0 GRID_MAX=5 GRID_BIN=250
# Now add together all the instantaneous histograms
hist: ACCUMULATE ARG=kde STRIDE=1
# And normalise the histogram
n: ACCUMULATE ARG=c STRIDE=1
a: CUSTOM ARG=hist,n FUNC=x/y PERIODIC=NO
# And print out the final histogram
DUMPGRID ARG=a FILE=histo.grid
```

This input computes 5 distances. Kernels correpsonding to all five of these distances are added to the instaneous
histogram that is constructed using the [KDE](KDE.md) action.  When we call the accumulate action here we are thus
not simply adding a single kernel to the accumulated grid when we add the the elements from `kde`.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace generic {

class Accumulate :
  public ActionWithValue,
  public ActionWithArguments,
  public ActionPilot {
private:
  bool clearnextstep;
  unsigned clearstride;
public:
  static void registerKeywords( Keywords& keys );
  Accumulate( const ActionOptions& );
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

PLUMED_REGISTER_ACTION(Accumulate,"ACCUMULATE")

void Accumulate::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  keys.use("UPDATE_FROM");
  keys.use("UPDATE_UNTIL");
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.addInputKeyword("compulsory","ARG","scalar/grid","the label of the argument that is being added to on each timestep");
  keys.add("compulsory","STRIDE","1","the frequency with which the data should be collected and added to the quantity being averaged");
  keys.add("compulsory","CLEAR","0","the frequency with which to clear all the accumulated data.  The default value "
           "of 0 implies that all the data will be used and that the grid will never be cleared");
  keys.setValueDescription("scalar/grid","a sum calculated from the time series of the input quantity");
}

Accumulate::Accumulate( const ActionOptions& ao ):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  ActionPilot(ao),
  clearnextstep(true) {
  if( getNumberOfArguments()!=1 ) {
    error("there should only be one argument to this action");
  }
  if( !getPntrToArgument(0)->hasDerivatives() && getPntrToArgument(0)->getRank()!=0 ) {
    error("input to the accumulate action should be a scalar or a grid");
  }

  parse("CLEAR",clearstride);
  if( clearstride>0 ) {
    if( clearstride%getStride()!=0 ) {
      error("CLEAR parameter must be a multiple of STRIDE");
    }
    log.printf("  clearing average every %u steps \n",clearstride);
  }
  std::vector<std::size_t> shape( getPntrToArgument(0)->getShape() );
  addValueWithDerivatives( shape );
  setNotPeriodic();
  if( getPntrToArgument(0)->isPeriodic() ) {
    error("you cannot accumulate a periodic quantity");
  }
}

unsigned Accumulate::getNumberOfDerivatives() {
  if( getPntrToArgument(0)->getRank()>0 ) {
    return getPntrToArgument(0)->getNumberOfGridDerivatives();
  }
  return getPntrToArgument(0)->getNumberOfDerivatives();
}

void Accumulate::update() {
  if( clearnextstep ) {
    if( getPntrToComponent(0)->getNumberOfValues()!=getPntrToArgument(0)->getNumberOfValues() ) {
      getPntrToComponent(0)->setShape( getPntrToArgument(0)->getShape() );
    }
    clearnextstep=false;
    getPntrToComponent(0)->set(0,0.0);
    getPntrToComponent(0)->clearDerivatives(true);
  }
  if( getStep()==0 ) {
    return;
  }

  Value* myarg=getPntrToArgument(0);
  Value* myout = getPntrToComponent(0);
  if( getPntrToArgument(0)->getRank()>0 ) {
    unsigned nvals = myarg->getNumberOfValues(), nder = myarg->getNumberOfGridDerivatives();
    for(unsigned i=0; i<nvals; ++i) {
      myout->set( i, myout->get(i) + myarg->get(i) );
      for(unsigned j=0; j<nder; ++j) {
        myout->addGridDerivatives( i, j, myarg->getGridDerivative( i, j ) );
      }
    }
  } else {
    getPntrToComponent(0)->add( getPntrToArgument(0)->get() );
  }

  // Clear if required
  if( clearstride>0 && getStep()%clearstride==0 ) {
    clearnextstep=true;
  }
}

}
}
