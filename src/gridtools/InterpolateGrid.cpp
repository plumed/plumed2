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
#include "EvaluateGridFunction.h"
#include "ActionWithGrid.h"

//+PLUMEDOC GRIDANALYSIS INTERPOLATE_GRID
/*
Interpolate a smooth function stored on a grid onto a grid with a smaller grid spacing.

This action takes a function evaluated on a grid as input and can be used to interpolate the values of that
function on to a finer grained grid.  By default the interpolation within this algorithm is done using splines.

## Examples

The input below can be used to post process a trajectory.  It calculates a [HISTOGRAM](HISTOGRAM.md) as a function the
distance between atoms 1 and 2 using kernel density estimation.  During the calculation the values of the kernels
are evaluated at 100 points on a uniform grid between 0.0 and 3.0.  Prior to outputting this function at the end of the
simulation this function is interpolated onto a finer grid of 200 points between 0.0 and 3.0.

```plumed
x: DISTANCE ATOMS=1,2
hA1: HISTOGRAM ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1
ii: INTERPOLATE_GRID ARG=hA1 GRID_BIN=200
DUMPGRID ARG=ii FILE=histo.dat
```

When specifying the parameters of the interpolating grid you can use GRID_BIN to specify the number of grid points
or GRID_SPACING to specify the spacing between grid points as shown below:

```plumed
x: DISTANCE ATOMS=1,2
hA1: HISTOGRAM ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1
ii: INTERPOLATE_GRID ARG=hA1 GRID_SPACING=0.05 INTERPOLATION_TYPE=floor
DUMPGRID ARG=ii FILE=histo.dat
```

In the above input spline interpolation is not used. Instead, if you evaluating the function at value $x$, which
is $x_n < x < x_{n+1}$ where $x_n$ and $x_{n+1}$ are two points where the values for the function on the input grid
are $f(x_n)$ and $f(x_{n+1})$ then $f(x)=f(x_n)$. If you want $f(x)=f(x_{n+1})$ you use `INTERPOLATION_TYPE=ceiling`.
Alternatively, you can use linear interpolation as has been done in the following input.

```plumed
x: DISTANCE ATOMS=1,2
hA1: HISTOGRAM ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1
ii: INTERPOLATE_GRID ARG=hA1 GRID_MIN=0.015 GRID_MAX=2.985 GRID_BIN=99 INTERPOLATION_TYPE=linear
DUMPGRID ARG=ii FILE=histo.dat
```

Notice that the way we have specified `GRID_MIN`, `GRID_MAX` and `GRID_BIN` in the above input ensures that function is evaluated
at the $n-1$ mid points that lie between the points on the input grid.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class InterpolateGrid : public ActionWithGrid {
public:
  using input_type = EvaluateGridFunction;
  using PTM = ParallelTaskManager<InterpolateGrid>;
private:
  bool firststep;
/// The parallel task manager
  PTM taskmanager;
  std::vector<std::size_t> nbin;
  std::vector<std::string> gmin, gmax;
  std::vector<double> gspacing;
  GridCoordinatesObject output_grid;
public:
  static void registerKeywords( Keywords& keys );
  explicit InterpolateGrid(const ActionOptions&ao);
  unsigned getNumberOfDerivatives() override ;
  const GridCoordinatesObject& getGridCoordinatesObject() const override ;
  std::vector<std::string> getGridCoordinateNames() const override ;
  void calculate() override ;
  void getInputData( std::vector<double>& inputdata ) const override ;
  // Calculate the value of the function at a grid point
  static void performTask( std::size_t task_index,
                           const EvaluateGridFunction& actiondata,
                           ParallelActionsInput& input,
                           ParallelActionsOutput& output );
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
  static int getNumberOfValuesPerTask( std::size_t task_index,
                                       const EvaluateGridFunction& actiondata );
  static void getForceIndices( std::size_t task_index,
                               std::size_t colno,
                               std::size_t ntotal_force,
                               const EvaluateGridFunction& actiondata,
                               const ParallelActionsInput& input,
                               ForceIndexHolder force_indices );
};

PLUMED_REGISTER_ACTION(InterpolateGrid,"INTERPOLATE_GRID")

void InterpolateGrid::registerKeywords( Keywords& keys ) {
  ActionWithGrid::registerKeywords( keys );
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.addInputKeyword("compulsory","ARG","grid","the label for function on the grid that you would like to interpolate");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
  keys.add("compulsory","GRID_MIN","auto","the lower bounds for the grid. By default this action uses the lower bound for the input grid");
  keys.add("compulsory","GRID_MAX","auto","the upper bounds for the grid. By default this action uses the upper bound for the input grid");
  EvaluateGridFunction::registerKeywords( keys );
  keys.remove("ZERO_OUTSIDE_GRID_RANGE");
  keys.setValueDescription("grid","the function evaluated onto the interpolated grid");
  PTM::registerKeywords( keys );
}

InterpolateGrid::InterpolateGrid(const ActionOptions&ao):
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

  parseVector("GRID_BIN",nbin);
  parseVector("GRID_SPACING",gspacing);
  unsigned dimension = getPntrToArgument(0)->getRank();
  gmin.resize( dimension );
  gmax.resize( dimension );
  parseVector("GRID_MIN",gmin);
  parseVector("GRID_MAX",gmax);
  if( nbin.size()!=dimension && gspacing.size()!=dimension ) {
    error("MIDPOINTS, GRID_BIN or GRID_SPACING must be set");
  }
  if( nbin.size()==dimension ) {
    log.printf("  number of bins in grid %ld", nbin[0]);
    for(unsigned i=1; i<nbin.size(); ++i) {
      log.printf(", %ld", nbin[i]);
    }
    log.printf("\n");
  } else if( gspacing.size()==dimension ) {
    log.printf("  spacing for bins in grid %f", gspacing[0]);
    for(unsigned i=1; i<gspacing.size(); ++i) {
      log.printf(", %f", gspacing[i]);
    }
    log.printf("\n");
  }
  // Create the input grid
  function::FunctionOptions foptions;
  EvaluateGridFunction::read( taskmanager.getActionInput(), this, foptions );
  // Need this for creation of tasks
  output_grid.setup( "flat", taskmanager.getActionInput().getPbc(), 0, 0.0 );

  // Now add a value
  std::vector<std::size_t> shape( getPntrToArgument(0)->getShape() );
  addValueWithDerivatives( shape );

  if( getPntrToArgument(0)->isPeriodic() ) {
    std::string min, max;
    getPntrToArgument(0)->getDomain( min, max );
    setPeriodic( min, max );
  } else {
    setNotPeriodic();
  }
  // Setup the task manager
  if( taskmanager.getActionInput().interpolation_type==EvaluateGridFunction::spline ) {
    taskmanager.setupParallelTaskManager( 0, 0 );
  } else if( taskmanager.getActionInput().interpolation_type==EvaluateGridFunction::linear ) {
    taskmanager.setupParallelTaskManager( 1+dimension, getPntrToArgument(0)->getNumberOfStoredValues() );
  } else if( taskmanager.getActionInput().interpolation_type==EvaluateGridFunction::floor || taskmanager.getActionInput().interpolation_type==EvaluateGridFunction::ceiling ) {
    taskmanager.setupParallelTaskManager( 1, getPntrToArgument(0)->getNumberOfStoredValues() );
  } else {
    error("interpolation type not defined");
  }
}

unsigned InterpolateGrid::getNumberOfDerivatives() {
  return getPntrToArgument(0)->getRank();
}

const GridCoordinatesObject& InterpolateGrid::getGridCoordinatesObject() const {
  return output_grid;
}

std::vector<std::string> InterpolateGrid::getGridCoordinateNames() const {
  ActionWithGrid* ag = ActionWithGrid::getInputActionWithGrid( getPntrToArgument(0)->getPntrToAction() );
  plumed_assert( ag );
  return ag->getGridCoordinateNames();
}

void InterpolateGrid::calculate() {
  if( firststep ) {
    ActionWithGrid* ag=ActionWithGrid::getInputActionWithGrid( getPntrToArgument(0)->getPntrToAction() );
    plumed_assert( ag );
    const GridCoordinatesObject& mygrid = ag->getGridCoordinatesObject();
    for(unsigned i=0; i<gmin.size(); ++i) {
      if( gmin[i]=="auto" ) {
        gmin[i] = mygrid.getMin()[i];
      }
      if( gmax[i]=="auto" ) {
        gmax[i] = mygrid.getMax()[i];
      }
    }
    output_grid.setBounds( gmin, gmax, nbin, gspacing );
    getPntrToComponent(0)->setShape( output_grid.getNbin(true) );
    if( taskmanager.getActionInput().interpolation_type==EvaluateGridFunction::linear ) {
      taskmanager.setupParallelTaskManager( 1+getPntrToComponent(0)->getRank(), getPntrToArgument(0)->getNumberOfStoredValues() );
    } else if( taskmanager.getActionInput().interpolation_type==EvaluateGridFunction::floor || taskmanager.getActionInput().interpolation_type==EvaluateGridFunction::ceiling ) {
      taskmanager.setupParallelTaskManager( 1, getPntrToArgument(0)->getNumberOfStoredValues() );
    } else if( taskmanager.getActionInput().interpolation_type!=EvaluateGridFunction::spline ) {
      error("interpolation type not defined");
    }
    firststep=false;
  }
  taskmanager.runAllTasks();
}

void InterpolateGrid::getInputData( std::vector<double>& inputdata ) const {
  std::size_t ndim = output_grid.getDimension();
  std::size_t nstored = getConstPntrToComponent(0)->getNumberOfStoredValues();
  std::vector<double> pos( ndim );
  if( inputdata.size()!=nstored*ndim ) {
    inputdata.resize( ndim*nstored );
  }

  for(unsigned i=0; i<nstored; ++i) {
    output_grid.getGridPointCoordinates( i, pos );
    for(unsigned j=0; j<ndim; ++j) {
      inputdata[ i*ndim + j ] = pos[j];
    }
  }
}

void InterpolateGrid::performTask( std::size_t task_index,
                                   const EvaluateGridFunction& actiondata,
                                   ParallelActionsInput& input,
                                   ParallelActionsOutput& output ) {
  auto funcout = function::FunctionOutput::create( 1,
                 output.values.data(),
                 input.ranks[0],
                 output.values.data()+1 );
  EvaluateGridFunction::calc( actiondata,
                              false,
                              View<const double>(input.inputdata + task_index*input.ranks[0],input.ranks[0]),
                              funcout );

  if( actiondata.interpolation_type==EvaluateGridFunction::linear ) {
    View<const double> args(input.inputdata + task_index*input.ranks[0],input.ranks[0]);
    const GridCoordinatesObject & gridobject = actiondata.getGridObject();
    unsigned dimension = gridobject.getDimension();
    std::vector<double> xfloor(dimension);
    std::vector<unsigned> indices(dimension), nindices(dimension), ind(dimension);
    gridobject.getIndices( args, indices );
    unsigned nn=gridobject.getIndex(indices);
    gridobject.getGridPointCoordinates( nn, nindices, xfloor );
    output.derivatives[0] = 0;
    for(unsigned i=0; i<dimension; ++i) {
      int x0=1;
      if(nindices[i]==indices[i]) {
        x0=0;
      }
      double ddx=gridobject.getGridSpacing()[i];
      double X = fabs((args[i]-xfloor[i])/ddx-(double)x0);
      output.derivatives[0] += (1-X);
      output.derivatives[1+i] = X;
    }
  } else if( actiondata.interpolation_type==EvaluateGridFunction::floor ) {
    output.derivatives[0] = 1;
  } else if( actiondata.interpolation_type==EvaluateGridFunction::ceiling ) {
    output.derivatives[0] = 1;
  }
}

void InterpolateGrid::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

int InterpolateGrid::getNumberOfValuesPerTask( std::size_t task_index,
    const EvaluateGridFunction& actiondata ) {
  return 1;
}

void InterpolateGrid::getForceIndices( std::size_t task_index,
                                       std::size_t colno,
                                       std::size_t ntotal_force,
                                       const EvaluateGridFunction& actiondata,
                                       const ParallelActionsInput& input,
                                       ForceIndexHolder force_indices ) {

  const GridCoordinatesObject & gridobject = actiondata.getGridObject();
  unsigned dimension = gridobject.getDimension();
  std::vector<unsigned> indices(dimension);
  gridobject.getIndices( View<const double>(input.inputdata + task_index*input.ranks[0],input.ranks[0]), indices );
  force_indices.threadsafe_derivatives_end[0] = 0;
  if( actiondata.interpolation_type==EvaluateGridFunction::linear ) {
    force_indices.tot_indices[0] = 1 + dimension;
    force_indices.indices[0][0] = gridobject.getIndex(indices);
    std::vector<unsigned> ind( dimension );
    for(unsigned i=0; i<dimension; ++i) {
      for(unsigned j=0; j<dimension; ++j) {
        ind[j] = indices[j];
      }
      if( gridobject.isPeriodic(i) && (ind[i]+1)==gridobject.getNbin(false)[i] ) {
        ind[i]=0;
      } else {
        ind[i] = ind[i] + 1;
      }
      force_indices.indices[0][1+i] = gridobject.getIndex(ind);
    }
  } else if( actiondata.interpolation_type==EvaluateGridFunction::floor ) {
    force_indices.indices[0][0] = gridobject.getIndex(indices);
    force_indices.tot_indices[0] = 1;
  } else if( actiondata.interpolation_type==EvaluateGridFunction::ceiling ) {
    for(unsigned i=0; i<dimension; ++i) {
      if( gridobject.isPeriodic(i) && (indices[i]+1)==gridobject.getNbin(false)[i] ) {
        indices[i]=0;
      } else {
        indices[i] = indices[i] + 1;
      }
    }
    force_indices.indices[0][0] = gridobject.getIndex(indices);
    force_indices.tot_indices[0] = 1;
  } else {
    plumed_merror("cannot calculate derivatives for this type of interpolation");
  }
}

}
}
