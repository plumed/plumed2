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
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ParallelTaskManager.h"
#include "ActionWithGrid.h"

//+PLUMEDOC GRIDANALYSIS MARGINAL
/*
Compute a marginal from the input multi-dimensional function

The marginals of the function $f(x,y)$ are defined as:

$$
g(y) = \int f(x,y) \textrm{d}x \qquad \textrm{and} \qquad h(x) = \int f(x,y) \textrm{d}y
$$

This action takes a function of two or more coordinates that is stored on a grid in input.  The marginal
of the input function is then computed via numerical integration.  To see how this works in practice
consider the following input:

```plumed
d1: DISTANCE ATOMS1=1,3 ATOMS2=1,4 ATOMS3=1,5 ATOMS4=1,6 COMPONENTS

k: KDE ARG=d1.x,d1.y BANDWIDTH=0.02,0.02 KERNEL=gaussian GRID_BIN=50,50 GRID_MIN=0,0 GRID_MAX=1,1
DUMPGRID ARG=k FILE=func_xy.grid STRIDE=1
m1: MARGINAL ARG=k DIR=1
DUMPGRID ARG=m1 FILE=func_x.grid STRIDE=1
m2: MARGINAL ARG=k DIR=2
DUMPGRID ARG=m1 FILE=func_y.grid STRIDE=1
```

This input uses the [KDE](KDE.md) action to compute a distribution from  the x and y components of the distances
that were computed using the [DISTANCE](DISTANCE.md) command.  We then output this two dimensional histogram and
use the MARGINAL command to compute the distributions for the x and y components of the vector separately.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class MarginalInput {
public:
/// This is the direction we are calculating the marginal along
  int dir;
/// This is the function whose marginal we are calculating
  Value* function=NULL;
/// This is the action that calculates the function we need
  ActionWithGrid* gridact=NULL;
  static void registerKeywords( Keywords& keys );
  static void read( MarginalInput& func, ActionWithArguments* action );
  MarginalInput& operator=( const MarginalInput& m ) {
    dir = m.dir;
    function = m.function;
    gridact = m.gridact;
    return *this;
  }
/// Get the vector containing the minimum value of the grid in each dimension
  std::vector<std::string> getMin() const ;
/// Get the vector containing the maximum value of the grid in each dimension
  std::vector<std::string> getMax() const ;
/// Get the periodicity of the grid
  std::vector<bool> getPbc() const ;
/// Get the number of grid points in each direction
  std::vector<std::size_t> getNbin() const ;
/// This gets the grid object
  const GridCoordinatesObject & getGridObject() const ;
};

void MarginalInput::registerKeywords( Keywords& keys ) {
  keys.add("compulsory","DIR","the direction along which the marginal is being computed");
}

std::vector<std::string> MarginalInput::getMin() const {
  std::vector<std::string> min(1, (gridact->getGridCoordinatesObject()).getMin()[dir] );
  return min;
}

std::vector<std::string> MarginalInput::getMax() const {
  std::vector<std::string> max(1, (gridact->getGridCoordinatesObject()).getMax()[dir] );
  return max;
}

std::vector<bool> MarginalInput::getPbc() const {
  std::vector<bool> pbc(1, (gridact->getGridCoordinatesObject()).isPeriodic(dir) );
  return pbc;
}

std::vector<std::size_t> MarginalInput::getNbin() const {
  std::vector<std::size_t> nbin(1, (gridact->getGridCoordinatesObject()).getNbin(false)[dir] );
  return nbin;
}

/// This gets the grid object
const GridCoordinatesObject & MarginalInput::getGridObject() const {
  return (gridact->getGridCoordinatesObject());
}

void MarginalInput::read( MarginalInput& func, ActionWithArguments* action ) {
  func.function = action->getPntrToArgument(0);
  if( func.function->getRank()==0 || !func.function->hasDerivatives() ) {
    action->error("should have one grid as input to this action");
  }
  if( func.function->getRank()==1 ) {
    action->error("cannot calculate marginal for one dimensional grid");
  }
  // Get the input grid
  func.gridact = ActionWithGrid::getInputActionWithGrid( func.function->getPntrToAction() );
  plumed_assert( func.gridact );
  if( func.gridact->getGridCoordinatesObject().getGridType()!="flat" ) {
    action->error("cannot calculate marginals on fibonacci sphere");
  }
  action->parse("DIR",func.dir);
  action->log.printf("  calculating marginal along %dth axis of input grid\n", func.dir );
  func.dir = func.dir - 1;
  if( func.dir<0 || static_cast<unsigned>(func.dir)>=(action->getPntrToArgument(0))->getRank() ) {
    action->error("not enough directions in input grid to compute desired marginal");
  }
}

class Marginal : public ActionWithGrid {
public:
  using input_type = MarginalInput;
  using PTM = ParallelTaskManager<Marginal>;
private:
  bool firststep;
/// The parallel task manager
  PTM taskmanager;
  GridCoordinatesObject output_grid;
public:
  static void registerKeywords( Keywords& keys );
  explicit Marginal(const ActionOptions&ao);
  unsigned getNumberOfDerivatives() override ;
  const GridCoordinatesObject& getGridCoordinatesObject() const override ;
  std::vector<std::string> getGridCoordinateNames() const override ;
  void calculate() override ;
  // Calculate the value of the function at a grid point
  static void performTask( std::size_t task_index,
                           const MarginalInput& actiondata,
                           ParallelActionsInput& input,
                           ParallelActionsOutput& output );
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
  static int getNumberOfValuesPerTask( std::size_t task_index,
                                       const MarginalInput& actiondata );
  static void getForceIndices( std::size_t task_index,
                               std::size_t colno,
                               std::size_t ntotal_force,
                               const MarginalInput& actiondata,
                               const ParallelActionsInput& input,
                               ForceIndexHolder force_indices );
};

PLUMED_REGISTER_ACTION(Marginal,"MARGINAL")

void Marginal::registerKeywords( Keywords& keys ) {
  ActionWithGrid::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","grid","the label for function on the grid that you would like to interpolate");
  keys.setValueDescription("grid","the function evaluated onto the interpolated grid");
  MarginalInput::registerKeywords( keys );
  PTM::registerKeywords( keys );
}

Marginal::Marginal(const ActionOptions&ao):
  Action(ao),
  ActionWithGrid(ao),
  firststep(true),
  taskmanager(this) {
  if( getNumberOfArguments()!=1 ) {
    error("should only be one argument to this action");
  }
  if( getPntrToArgument(0)->getRank()==0 || !getPntrToArgument(0)->hasDerivatives() ) {
    error("input to this action should be a grid");
  }
  if( getPntrToArgument(0)->getRank()==1 ) {
    error("should have more than one axis in grid to compute a marginal");
  }

  // Create the input grid
  MarginalInput::read( taskmanager.getActionInput(), this );
  // Need this for creation of tasks
  output_grid.setup( "flat", taskmanager.getActionInput().getPbc(), 0, 0.0 );

  // Now add a value
  std::vector<std::size_t> shape( 1, getPntrToArgument(0)->getShape()[taskmanager.getActionInput().dir] );
  addValueWithDerivatives( shape );

  if( getPntrToArgument(0)->isPeriodic() ) {
    std::string min, max;
    getPntrToArgument(0)->getDomain( min, max );
    setPeriodic( min, max );
  } else {
    setNotPeriodic();
  }
}

unsigned Marginal::getNumberOfDerivatives() {
  return 1;
}

const GridCoordinatesObject& Marginal::getGridCoordinatesObject() const {
  return output_grid;
}

std::vector<std::string> Marginal::getGridCoordinateNames() const {
  ActionWithGrid* ag = ActionWithGrid::getInputActionWithGrid( getPntrToArgument(0)->getPntrToAction() );
  plumed_assert( ag );
  std::vector<std::string> gnames(1);
  gnames[0] = ag->getGridCoordinateNames()[taskmanager.getActionInput().dir];
  return gnames;
}

void Marginal::calculate() {
  if( firststep ) {
    std::vector<double> gspacing;
    output_grid.setBounds( taskmanager.getActionInput().getMin(), taskmanager.getActionInput().getMax(), taskmanager.getActionInput().getNbin(), gspacing );
    getPntrToComponent(0)->setShape( output_grid.getNbin(true) );
    taskmanager.setupParallelTaskManager( taskmanager.getActionInput().function->getNumberOfValues()/getPntrToComponent(0)->getNumberOfValues(), 0 );
    firststep=false;
  }
  taskmanager.runAllTasks();
}

void Marginal::performTask( std::size_t task_index,
                            const MarginalInput& actiondata,
                            ParallelActionsInput& input,
                            ParallelActionsOutput& output ) {
  const GridCoordinatesObject & gridobject = actiondata.getGridObject();
  std::vector<unsigned> indices( gridobject.getDimension() );
  double volelement = 1;
  for(unsigned i=0; i<indices.size(); ++i ) {
    if( static_cast<int>(i)==actiondata.dir ) {
      continue;
    }
    volelement *= gridobject.getGridSpacing()[i];
  }
  std::size_t k = 0;
  for(unsigned i=0; i<actiondata.function->getNumberOfValues(); ++i) {
    gridobject.getIndices( i, indices );
    if( indices[actiondata.dir]!=task_index ) {
      continue;
    }
    output.values[0] += actiondata.function->get(i);
    output.values[1] += actiondata.function->getGridDerivative(i,0);
    output.derivatives[k] = volelement;
    k++;
  }
  output.values[0] *= volelement;
  output.values[1] *= volelement;
}

void Marginal::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

int Marginal::getNumberOfValuesPerTask( std::size_t task_index,
                                        const MarginalInput& actiondata ) {
  return 1;
}

void Marginal::getForceIndices( std::size_t task_index,
                                std::size_t colno,
                                std::size_t ntotal_force,
                                const MarginalInput& actiondata,
                                const ParallelActionsInput& input,
                                ForceIndexHolder force_indices ) {
  std::size_t k = 0;
  const GridCoordinatesObject & gridobject = actiondata.getGridObject();
  std::vector<unsigned> indices( gridobject.getDimension() );
  for(unsigned i=0; i<actiondata.function->getNumberOfValues(); ++i) {
    gridobject.getIndices( i, indices );
    if( indices[actiondata.dir]!=task_index ) {
      continue;
    }
    force_indices.indices[0][k] = i;
    k++;
  }
  force_indices.threadsafe_derivatives_end[0] = k;
  force_indices.tot_indices[0] = k;
}

}
}
