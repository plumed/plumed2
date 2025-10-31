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
#include "ContourFindingObject.h"
#include "gridtools/ActionWithGrid.h"
#include "gridtools/EvaluateGridFunction.h"
#include "core/ParallelTaskManager.h"

//+PLUMEDOC GRIDANALYSIS FIND_CONTOUR_SURFACE
/*
Find an isocontour by searching along either the x, y or z direction.

As discussed in the documentation for the [gridtools](module_gridtools.md), PLUMED contains a number of tools that allow you to calculate
a function on a grid.  The function on this grid might be a [HISTOGRAM](HISTOGRAM.md)  or it might be one of the phase fields that are
discussed [here](module_contour.md).  If this function has one or two input
arguments it is relatively straightforward to plot the function.  If by contrast the data has a three dimensions it can be
difficult to visualize.

This action provides one tool for visualizing these functions.  It can be used to search for a set of points on a contour
where the function takes a particular value.  In other words, for the function $f(x,y,z)$ this action would find a set
of points $\{x_c,y_c,z_c\}$ that have:

$$
f(x_c,y_c,z_c) - c = 0
$$

where $c$ is some constant value that is specified by the user.  The points on this contour are find by searching along lines
that run parallel to the $x$, $y$ or $z$ axis of the simulation cell.  The result is, therefore, a two dimensional
function evaluated on a grid that gives us the height of the interface as a function of two coordinates.

It is important to note that this action can only be used to detect contours in three dimensional functions.  In addition, this action will fail to
find the full set of contour  points if the contour does not have the same topology as an infinite plane.  If you are uncertain that the isocontours in your
function have the appropriate topology you should use [FIND_CONTOUR](FIND_CONTOUR.md) in place of this action.

## Examples

The input shown below was used to analyze the results from a simulation of an interface between solid and molten Lennard Jones.  The interface between
the solid and the liquid was set up in the plane perpendicular to the $z$ direction of the simulation cell. The input below calculates something
akin to a Willard-Chandler dividing surface (see [contour](module_contour.md)) between the solid phase and the liquid phase.  There are two of these interfaces within the
simulation box because of the periodic boundary conditions but we were able to determine that one of these two surfaces lies in a particular part of the
simulation box.  The input below detects the height profile of one of these two interfaces.  It does so by computing a phase field average from the values, $s_i$, of the
[FCCUBIC](FCCUBIC.md) symmetry functions for each of the atoms using the following expression.

$$
\rho'(x,y,z) = \frac{ \sum_{i=1}^N s_i K\left(\frac{x-x_i}{\lambda}, \frac{y-y_i}{\lambda}, \frac{z-z_i}{\lambda}\right) }{ \sum_{i=1}^N K\left(\frac{x-x_i}{\lambda}, \frac{y-y_i}{\lambda}, \frac{z-z_i}{\lambda}\right) }
$$

where $(x_i, y_i, z_i)$ is the position of atom $i$ relative to the position of atom 1, $K$ is a Gaussian kernel function and $\lambda=1.0$.

Notice that we use the fact that we know roughly where the interface is when specifying how this phase field is to be calculated and specify the region over the $z$-axis
in which the [KDE](KDE.md) is computed.  Once we have calculated the phase field we search for contour points on the lines that run parallel to the $z$-direction of the cell
box using the FIND_CONTOUR_SURFACE command.  The final result is a $14 \times 14$ grid of values for the height of the interface as a function of the $(x,y)$
position.  This grid is then output to a file called `contour2.dat`.

Notice that the commands below calculate the instantaneous position of the surface separating the solid and liquid and that as such the accumulated average is cleared
on every step.

```plumed
UNITS NATURAL

# This calculates the value of a set of symmetry functions for the atoms of interest
fcc: FCCUBIC ...
  SPECIES=1-96000 SWITCH={CUBIC D_0=1.2 D_MAX=1.5}
  ALPHA=27 PHI=0.0 THETA=-1.5708 PSI=-2.35619
...

# This determines the positions of the atoms of interest relative to the position of atom 1
dens2_dist: DISTANCES ORIGIN=1 ATOMS=fcc COMPONENTS
# This computes the numerator in the expression above for the phase field
dens2_numer: KDE ...
   VOLUMES=fcc_n ARG=dens2_dist.x,dens2_dist.y,dens2_dist.z
   GRID_BIN=14,14,50 GRID_MIN=auto,auto,6.0
   GRID_MAX=auto,auto,11.0 BANDWIDTH=1.0,1.0,1.0
...
# This computes the denominator
dens2_denom: KDE ...
   ARG=dens2_dist.x,dens2_dist.y,dens2_dist.z
   GRID_BIN=14,14,50 GRID_MIN=auto,auto,6.0
   GRID_MAX=auto,auto,11.0 BANDWIDTH=1.0,1.0,1.0
...
# This computes the final phase field
dens2: CUSTOM ARG=dens2_numer,dens2_denom FUNC=x/y PERIODIC=NO

# We can now find and print the location of the two dimensional contour surface
ss2: FIND_CONTOUR_SURFACE ARG=dens2 CONTOUR=0.42 SEARCHDIR=dens2_dist.z
DUMPGRID ARG=ss2 FILE=contour2.dat STRIDE=1
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace contour {

class FindContourSurfaceObject {
public:
  unsigned dir_n;
  std::vector<unsigned> gdirs;
  std::vector<double> direction;
  ContourFindingObject<gridtools::EvaluateGridFunction> cf;
  static void registerKeywords( Keywords& keys ) {
    ContourFindingObject<gridtools::EvaluateGridFunction>::registerKeywords( keys );
    keys.add("compulsory","SEARCHDIR","In which directions do you wish to search for the contour.");
  }
  static void read( FindContourSurfaceObject& func, ActionWithArguments* action, function::FunctionOptions& options ) {
    ContourFindingObject<gridtools::EvaluateGridFunction>::read( func.cf, action, options );
    std::string dir;
    action->parse("SEARCHDIR",dir);
    action->log.printf("  calculating location of contour on %d dimensional grid \n", (action->getPntrToArgument(0))->getRank()-1 );
    Value* gval=func.cf.function.function;
    gridtools::ActionWithGrid* ag = gridtools::ActionWithGrid::getInputActionWithGrid( gval->getPntrToAction() );
    if( !ag ) {
      action->error("input argument must be a grid");
    }
    std::vector<std::string> argn( ag->getGridCoordinateNames() );
    func.gdirs.resize( gval->getRank()-1 );
    unsigned n=0;
    for(unsigned i=0; i<gval->getRank(); ++i) {
      if( argn[i]==dir ) {
        func.dir_n=i;
      } else {
        if( n==func.gdirs.size() ) {
          action->error("could not find " + dir + " direction in input grid");
        }
        func.gdirs[n]=i;
        n++;
      }
    }
    if( n!=(gval->getRank()-1) ) {
      action->error("output of grid is not understood");
    }
  }
};

class FindContourSurface : public gridtools::ActionWithGrid {
public:
  using input_type = FindContourSurfaceObject;
  using PTM = ParallelTaskManager<FindContourSurface>;
private:
  bool firststep;
/// The parallel task manager
  PTM taskmanager;
  std::vector<std::string> gnames;
  gridtools::GridCoordinatesObject gridcoords;
/// Get the input grid object for the grid that we are finding the contour in
  const gridtools::GridCoordinatesObject& getInputGridObject() const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit FindContourSurface(const ActionOptions&ao);
  unsigned getNumberOfDerivatives() override ;
  std::vector<std::string> getGridCoordinateNames() const override ;
  const gridtools::GridCoordinatesObject& getGridCoordinatesObject() const override ;
  void calculate() override ;
  void getInputData( std::vector<double>& inputdata ) const override;
  static void performTask( std::size_t task_index,
                           const FindContourSurfaceObject& actiondata,
                           ParallelActionsInput& input,
                           ParallelActionsOutput& output );
};

PLUMED_REGISTER_ACTION(FindContourSurface,"FIND_CONTOUR_SURFACE")

void FindContourSurface::registerKeywords( Keywords& keys ) {
  ActionWithGrid::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","grid","the labels of the grid in which the contour will be found");
  FindContourSurfaceObject::registerKeywords( keys );
  keys.setValueDescription("grid","a grid containing the location of the points in the Willard-Chandler surface along the chosen direction");
  keys.addDOI("10.1088/1361-648X/aa893d");
  PTM::registerKeywords( keys );
}

FindContourSurface::FindContourSurface(const ActionOptions&ao):
  Action(ao),
  ActionWithGrid(ao),
  firststep(true),
  taskmanager(this) {
  if( getPntrToArgument(0)->getRank()<2 ) {
    error("cannot find dividing surface if input grid is one dimensional");
  }

  Value* gval=getPntrToArgument(0);
  gnames.resize( getPntrToArgument(0)->getRank()-1 );

  gridtools::ActionWithGrid* ag = ActionWithGrid::getInputActionWithGrid( gval->getPntrToAction() );
  if( !ag ) {
    error("input argument must be a grid");
  }
  const gridtools::GridCoordinatesObject& mygridobj = ag->getGridCoordinatesObject();
  if( mygridobj.getGridType()=="fibonacci") {
    error("cannot search for contours in fibonacci grids");
  }
  // Now add a value
  std::vector<std::size_t> shape( mygridobj.getDimension()-1 );
  addValueWithDerivatives( shape );
  setNotPeriodic();
  // Read in the action input
  function::FunctionOptions options;
  FindContourSurfaceObject::read( taskmanager.getActionInput(), this, options );
  // Prepare the grid stuff for this action
  std::vector<bool> ipbc( mygridobj.getDimension()-1 );
  for(unsigned i=0; i<gnames.size(); ++i) {
    ipbc[i] = mygridobj.isPeriodic(taskmanager.getActionInput().gdirs[i]);
    gnames[i] = ag->getGridCoordinateNames()[taskmanager.getActionInput().gdirs[i]];
  }
  gridcoords.setup( "flat", ipbc, 0, 0.0 );
  log.printf("  calculating dividing surface along which function equals %f \n", taskmanager.getActionInput().cf.contour);
  taskmanager.setupParallelTaskManager( 0, 0 );
}

const gridtools::GridCoordinatesObject& FindContourSurface::getInputGridObject() const {
  return taskmanager.getActionInput().cf.function.getGridObject();
}

void FindContourSurface::calculate() {
  if( firststep ) {
    std::vector<double> fspacing;
    std::vector<std::size_t> snbins( gridcoords.getDimension() );
    std::vector<std::string> smin( gridcoords.getDimension() ), smax( gridcoords.getDimension() );
    for(unsigned i=0; i<taskmanager.getActionInput().gdirs.size(); ++i) {
      smin[i]=getInputGridObject().getMin()[taskmanager.getActionInput().gdirs[i]];
      smax[i]=getInputGridObject().getMax()[taskmanager.getActionInput().gdirs[i]];
      snbins[i]=getInputGridObject().getNbin(false)[taskmanager.getActionInput().gdirs[i]];
    }
    gridcoords.setBounds( smin, smax, snbins, fspacing );
    getPntrToComponent(0)->setShape( gridcoords.getNbin(true) );

    std::vector<unsigned> find( gridcoords.getDimension() );
    std::vector<unsigned> ind( gridcoords.getDimension() );
    for(unsigned i=0; i<gridcoords.getNumberOfPoints(); ++i) {
      find.assign( find.size(), 0 );
      gridcoords.getIndices( i, ind );
      for(unsigned j=0; j<taskmanager.getActionInput().gdirs.size(); ++j) {
        find[taskmanager.getActionInput().gdirs[j]]=ind[j];
      }
    }

    // Set the direction in which to look for the contour
    taskmanager.getActionInput().direction.resize( getInputGridObject().getDimension(), 0 );
    taskmanager.getActionInput().direction[taskmanager.getActionInput().dir_n] = 0.999999999*getInputGridObject().getGridSpacing()[taskmanager.getActionInput().dir_n];
    firststep=false;
  }
  taskmanager.runAllTasks();
}

unsigned FindContourSurface::getNumberOfDerivatives() {
  return gridcoords.getDimension();
}

std::vector<std::string> FindContourSurface::getGridCoordinateNames() const {
  return gnames;
}

const gridtools::GridCoordinatesObject& FindContourSurface::getGridCoordinatesObject() const {
  return gridcoords;
}

void FindContourSurface::getInputData( std::vector<double>& inputdata ) const {
#ifndef DNDEBUG
  std::size_t rank = gridcoords.getDimension();
  std::size_t ndata = gridcoords.getNumberOfPoints();
  if( inputdata.size()!=ndata*rank ) {
    inputdata.resize( ndata*rank );
  }
  std::vector<double> point( rank );
  for(unsigned i=0; i<ndata; ++i) {
    gridcoords.getGridPointCoordinates( i, point );
    for(unsigned j=0; j<rank; ++j) {
      inputdata[i*rank + j] = point[j];
    }
  }
#endif
}

void FindContourSurface::performTask( std::size_t task_index,
                                      const FindContourSurfaceObject& actiondata,
                                      ParallelActionsInput& input,
                                      ParallelActionsOutput& output ) {
  unsigned nfound=0;
  const gridtools::GridCoordinatesObject& gridobj = actiondata.cf.function.getGridObject();
  std::size_t nbins = gridobj.getNbin(false)[actiondata.dir_n];
  std::size_t shiftn=task_index;
  std::vector<unsigned> ind( gridobj.getDimension() );
  std::vector<double> point( gridobj.getDimension() );
#ifndef DNDEBUG
  std::size_t rank = gridobj.getDimension()-1;
  View<const double> oind( input.inputdata+task_index*rank, rank );
#endif
  for(unsigned i=0; i<nbins; ++i) {
#ifndef DNDEBUG
    std::vector<double> base_ind( gridobj.getDimension() );
    gridobj.getGridPointCoordinates( shiftn, base_ind );
    for(unsigned j=0; j<actiondata.gdirs.size(); ++j) {
      plumed_dbg_assert( base_ind[actiondata.gdirs[j]]==oind[j] );
    }
#endif
    // Get the index of the current grid point
    gridobj.getIndices( shiftn, ind );
    // Exit if we are at the edge of the grid
    if( !gridobj.isPeriodic(actiondata.dir_n) && (ind[actiondata.dir_n]+1)==nbins ) {
      shiftn += gridobj.getStride()[actiondata.dir_n];
      continue;
    }

    // Now get the function value at two points
    double val1=actiondata.cf.function.function->get( shiftn ) - actiondata.cf.contour;
    double val2;
    if( (ind[actiondata.dir_n]+1)==nbins ) {
      val2 = actiondata.cf.function.function->get( task_index ) - actiondata.cf.contour;
    } else {
      val2 = actiondata.cf.function.function->get( shiftn + gridobj.getStride()[actiondata.dir_n] ) - actiondata.cf.contour;
    }

    // Check if the minimum is bracketed
    if( val1*val2<0 ) {
      gridobj.getGridPointCoordinates( shiftn, point );
      ContourFindingObject<gridtools::EvaluateGridFunction>::findContour( actiondata.cf, actiondata.direction, point );
      output.values[0] = point[actiondata.dir_n];
      nfound++;
      break;
    }

    // This moves us on to the next point
    shiftn += gridobj.getStride()[actiondata.dir_n];
  }
  if( nfound==0 ) {
    plumed_merror("failed to find required grid point");
  }
}

}
}
