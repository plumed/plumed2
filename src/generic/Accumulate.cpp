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
#include "core/ActionWithGrid.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "tools/GridCoordinatesObject.h"
#include "tools/Checkpoint.h"

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
kde: KDE ARG=d BANDWIDTH=0.05 GRID_SPACING=0.01
# Now add together all the instantaneous histograms
hist: ACCUMULATE ARG=kde STRIDE=1
# And normalise the histogram
n: ACCUMULATE ARG=c STRIDE=1
a: CUSTOM ARG=hist,n FUNC=x/y PERIODIC=NO
# And print out the final histogram
DUMPGRID ARG=a FILE=histo.grid
```

At first glance the fact that we use a [KDE](KDE.md) action to construct an instantaneous histogram from a single
distance may appear odd.  The reason for doing this, however, is to introduce a clear distinction between
the syntaxes that are used for spatial and temporal averaging.  To see what I mean consider the following input:


```plumed
c: CONSTANT VALUE=5
d: DISTANCE ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6 ATOMS4=7,8 ATOMS5=9,10
# Construct the instantaneous histogram from the instantaneous value of the distance
kde: KDE ARG=d BANDWIDTH=0.05 GRID_SPACING=0.01
# Now add together all the instantaneous histograms
hist: ACCUMULATE ARG=kde STRIDE=1
# And normalise the histogram
n: ACCUMULATE ARG=c STRIDE=1
a: CUSTOM ARG=hist,n FUNC=x/y PERIODIC=NO
# And print out the final histogram
DUMPGRID ARG=a FILE=histo.grid
```

This input computes 5 distances. Kernels corresponding to all five of these distances are added to the instantaneous
histogram that is constructed using the [KDE](KDE.md) action.  When we call the accumulate action here we are thus
not simply adding a single kernel to the accumulated grid when we add the elements from `kde`.

## Restart

With [CHECKPOINT](CHECKPOINT.md) the accumulated data (scalar or grid) is stored in the checkpoint, so that sums, averages and histograms continue exactly after a restart.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace generic {

class Accumulate :
  public ActionWithValue,
  public ActionWithArguments,
  public ActionPilot,
  public GridOutput {
private:
  bool clearnextstep;
  unsigned clearstride;
/// The grid of the accumulated data. It is the grid of the argument when the data are accumulated; when the grid of
/// the argument is extended (see KDE) the data are placed in the new grid at the next accumulation.
  GridCoordinatesObject gridcoords;
  bool gridset=false;
/// the grid of the accumulated data becomes the one of the argument
  void copyGrid();
  void extendGrid();
public:
  std::vector<std::string> getGridCoordinateNames() const override {
    return ActionWithGrid::getInputActionWithGrid( getPntrToArgument(0)->getPntrToAction() )->getGridCoordinateNames();
  }
  const GridCoordinatesObject& getGridCoordinatesObject() const override {
    return gridcoords;
  }
  static void registerKeywords( Keywords& keys );
  bool hasCheckpointState() const override {
    return true;
  }
  void writeCheckpoint(CheckpointWriter& w) override;
  void readCheckpoint(CheckpointReader& r) override;
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
  keys.addFlag("DO_AFTER_FINAL_STEP",false,"do a final accumulation step at end the trajectory even if the trajectory's length is not a multiple of the stride");
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

  if( getStride()==0 ) {
    error("STRIDE should be a positive number of steps");
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
  if( shape.size()>0 ) {
    if( !dynamic_cast<const GridOutput*>( getPntrToArgument(0)->getPntrToAction() ) ) {
      error("the input grid should be calculated by an action that outputs a grid (e.g. KDE)");
    }
    // the type, dimension and periodicity of the grid (the bounds are set at the first accumulation)
    gridcoords=ActionWithGrid::getInputActionWithGrid( getPntrToArgument(0)->getPntrToAction() )->getGridCoordinatesObject();
  }
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

// The accumulated data is the state that survives from one step to the next (this
// is the value of type "average"): the data of the output value, with the
// derivatives for grids, and whether it has to be cleared at the next update
void Accumulate::writeCheckpoint(CheckpointWriter& w) {
  w.write("clearnextstep",(long long) clearnextstep);
  if( clearnextstep ) {
    return;
  }
  // the grid of the accumulated data, which can differ from the one of the argument at the restart (see extendGrid)
  if( gridset && gridcoords.getGridType()=="flat" ) {
    w.write("grid_set",1LL);
    gridcoords.writeCheckpoint( w, "grid_" );
  }
  const Value* myout = getPntrToComponent(0);
  std::vector<long long> shape;
  for(const auto& d : myout->getShape()) {
    shape.push_back(d);
  }
  w.write("shape",shape);
  std::vector<double> data;
  const std::size_t nvals = myout->getNumberOfValues();
  for(std::size_t i=0; i<nvals; ++i) {
    data.push_back( myout->get(i) );
  }
  w.write("data",data);
  if( myout->getRank()>0 ) {
    std::vector<double> der;
    const unsigned nder = myout->getNumberOfGridDerivatives();
    for(std::size_t i=0; i<nvals; ++i) {
      for(unsigned j=0; j<nder; ++j) {
        der.push_back( myout->getGridDerivative(i,j) );
      }
    }
    w.write("derivatives",der);
  }
}

void Accumulate::readCheckpoint(CheckpointReader& r) {
  long long cl;
  r.read("clearnextstep",cl);
  clearnextstep = cl!=0;
  if( clearnextstep ) {
    return;
  }
  if( r.has("grid_set") ) {
    gridcoords.readCheckpoint( r, "grid_" );
    gridset=true;
  }
  Value* myout = getPntrToComponent(0);
  std::vector<long long> shape;
  r.read("shape",shape);
  std::vector<std::size_t> ss(shape.begin(),shape.end());
  if( ss!=getPntrToArgument(0)->getShape() ) {
    error("the checkpoint was written with a different shape of the argument");
  }
  if( myout->getShape()!=ss ) {
    myout->setShape(ss);
  }
  std::vector<double> data;
  r.read("data",data);
  const std::size_t nvals = myout->getNumberOfValues();
  if( data.size()!=nvals ) {
    error("the checkpoint was written with a different number of values");
  }
  // the data of a grid holds the values and the derivatives: clear first, then restore
  myout->clearDerivatives(true);
  for(std::size_t i=0; i<nvals; ++i) {
    myout->set( i, data[i] );
  }
  if( myout->getRank()>0 ) {
    std::vector<double> der;
    r.read("derivatives",der);
    const unsigned nder = myout->getNumberOfGridDerivatives();
    if( der.size()!=nvals*nder ) {
      error("the checkpoint was written with a different number of derivatives");
    }
    for(std::size_t i=0; i<nvals; ++i) {
      for(unsigned j=0; j<nder; ++j) {
        myout->addGridDerivatives( i, j, der[i*nder+j] );
      }
    }
  }
}

// The grid of the argument has been extended (e.g. KDE, when the kernels do not fit in it): the accumulated data
// are placed in the new grid, which has the same spacing and contains the old one.
void Accumulate::extendGrid() {
  Value* myarg=getPntrToArgument(0);
  Value* myout=getPntrToComponent(0);
  const std::size_t nold=myout->getNumberOfValues();
  const unsigned nder=myout->getNumberOfGridDerivatives();
  std::vector<double> old( nold*(1+nder) );
  bool zero=true;
  for(std::size_t i=0; i<nold; ++i) {
    old[i*(1+nder)]=myout->get(i);
    for(unsigned j=0; j<nder; ++j) {
      old[i*(1+nder)+1+j]=myout->getGridDerivative(i,j);
    }
  }
  for(const auto v : old) {
    if( v!=0.0 ) {
      zero=false;
      break;
    }
  }
  const std::vector<std::size_t> oldshape( myout->getShape() );
  myout->setShape( myarg->getShape() );
  myout->clearDerivatives(true);
  const GridCoordinatesObject oldgrid( gridcoords );
  copyGrid();
  if( zero ) {
    return;
  }
  const GridCoordinatesObject& gc=gridcoords;
  const std::size_t ndim=myout->getRank();
  if( gc.getGridType()!="flat" || oldgrid.getNumberOfPoints()!=nold ) {
    error("the grid of the argument has changed: the accumulated data cannot be placed in it");
  }
  // the grid of the argument has been extended on the same lattice: the shift of the old points is exact
  if( !gc.onSameLattice( oldgrid ) ) {
    error("the grid of the argument has changed spacing or lattice: the accumulated data cannot be placed in it");
  }
  const std::vector<long long> offset( gc.getShiftFrom( oldgrid ) );
  const std::vector<unsigned> oldnbin( oldshape.begin(), oldshape.end() );
  const std::vector<std::size_t>& newshape=myout->getShape();
  std::vector<unsigned> oldind( ndim ), newind( ndim );
  for(std::size_t i=0; i<nold; ++i) {
    gc.convertIndexToIndices( i, oldnbin, oldind );
    bool inside=true;
    for(unsigned j=0; j<ndim; ++j) {
      const long long k=static_cast<long long>(oldind[j])+offset[j];
      if( k<0 || k>=static_cast<long long>(newshape[j]) ) {
        inside=false;
        break;
      }
      newind[j]=k;
    }
    const double* const point=old.data()+i*(1+nder);
    if( !inside ) {
      for(unsigned j=0; j<=nder; ++j) {
        if( point[j]!=0.0 ) {
          error("the grid of the argument does not contain the accumulated data");
        }
      }
      continue;
    }
    const unsigned ind=gc.getIndex( newind );
    myout->set( ind, point[0] );
    for(unsigned j=0; j<nder; ++j) {
      myout->setGridDerivatives( ind, j, point[1+j] );
    }
  }
}

void Accumulate::copyGrid() {
  gridcoords=ActionWithGrid::getInputActionWithGrid( getPntrToArgument(0)->getPntrToAction() )->getGridCoordinatesObject();
  gridset=true;
}

void Accumulate::update() {
  Value* myarg=getPntrToArgument(0);
  Value* myout = getPntrToComponent(0);
  if( clearnextstep ) {
    if( myout->getNumberOfValues()!=myarg->getNumberOfValues() ) {
      myout->setShape( myarg->getShape() );
    }
    clearnextstep=false;
    myout->set(0,0.0);
    myout->clearDerivatives(true);
    if( myarg->getRank()>0 ) {
      copyGrid();
    }
  }
  if( getStep()==0 ) {
    return;
  }

  if( myarg->getRank()>0 ) {
    if( myarg->getNumberOfValues()!=myout->getNumberOfValues() ) {
      extendGrid();
    }
    const std::size_t nvals = myarg->getNumberOfValues();
    const unsigned nder = myarg->getNumberOfGridDerivatives();
    for(std::size_t i=0; i<nvals; ++i) {
      myout->set( i, myout->get(i) + myarg->get(i) );
      for(unsigned j=0; j<nder; ++j) {
        myout->addGridDerivatives( i, j, myarg->getGridDerivative( i, j ) );
      }
    }
    if( !gridset ) {
      copyGrid();
    }
  } else {
    myout->add( myarg->get() );
  }

  // Clear if required
  if( clearstride>0 && getStep()%clearstride==0 ) {
    clearnextstep=true;
  }
}

}
}
