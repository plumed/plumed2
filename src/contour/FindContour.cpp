/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "FindContour.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"

//+PLUMEDOC GRIDANALYSIS FIND_CONTOUR
/*
Find an isocontour in a smooth function.

As discussed in the documentation for the [gridtools](module_gridtools.md), PLUMED contains a number of tools that allow you to calculate
a function on a grid.  The function on this grid might be a [HISTOGRAM](HISTOGRAM.md)  or it might be one of the phase fields that are
discussed [here](module_contour.md).  If this function has one or two input
arguments it is relatively straightforward to plot the function.  If by contrast the data has a three dimensions it can be
difficult to visualize.

This action provides one tool for visualizing these functions.  It can be used to search for a set of points on a contour
where the function takes a particular values.  In other words, for the function $f(x,y)$ this action would find a set
of points $\{x_c,y_c\}$ that have:

$$
f(x_c,y_c) - c = 0
$$

where $c$ is some constant value that is specified by the user.  The points on this contour are detected using a variant
on the [marching squares](https://en.wikipedia.org/wiki/Marching_squares) or [marching cubes](https://en.wikipedia.org/wiki/Marching_cubes) algorithm.
As such, and unlike [FIND_CONTOUR_SURFACE](FIND_CONTOUR_SURFACE.md) or [FIND_SPHERICAL_CONTOUR](FIND_SPHERICAL_CONTOUR.md), the function input to this action is not required to have three arguments.
You can find a contour in any function with 2 or more arguments.  Furthermore, the topology of the contour will be determined by the algorithm and does not need to be specified by the user.

## Examples

The input shown below was used to analyze the results from a simulation of an interface between solid and molten Lennard Jones.  The interface between
the solid and the liquid was set up in the plane perpendicular to the $z$ direction of the simulation cell. The input below calculates something
akin to a Willard-Chandler dividing surface (see [contour](module_contour.md)) between the solid phase and the liquid phase.  There are two of these interfaces within the
simulation box because of the periodic boundary conditions but we were able to determine that one of these two surfaces lies in a particular part of the
simulation box.  The input below detects the height profile of one of these two interfaces.  It does so by computing a phase field average from the transformed values, $f(s_i)$, of the
[FCCUBIC](FCCUBIC.md) symmetry functions for each of the atoms using the following expression.

$$
\rho'(x,y,z) = \frac{ \sum_{i=1}^N f(s_i) K\left(\frac{x-x_i}{\lambda}, \frac{y-y_i}{\lambda}, \frac{z-z_i}{\lambda}\right) }{ \sum_{i=1}^N K\left(\frac{x-x_i}{\lambda}, \frac{y-y_i}{\lambda}, \frac{z-z_i}{\lambda}\right) }
$$

where $(x_i, y_i, z_i)$ is the position of atom $i$ relative to the position of atom 1, $f$ is a switching function, $K$ is a Gaussian kernel function and $\lambda=1.0$.

The virtual atom that is computed using [CENTER](CENTER.md) is located in the center of the region where the atoms are solid like and the positions of the atoms
$(x_i, y_i, z_i)$ are given relative to this position when computing the phase field in the expression above.  This phase field provides a continuous function that
gives a measure of the average degree of solidness at each point in the simulation cell.  The Willard-Chandler dividing surface is calculated by finding a a set of points
at which the value of this phase field is equal to 0.5.  This set of points is output to file called `mycontour.dat`.  A new contour
is found on every single step for each frame that is read in.

```plumed
UNITS NATURAL

# This calculates the value of a set of symmetry functions for the atoms of interest
fcc: FCCUBIC ...
  SPECIES=1-96000 SWITCH={CUBIC D_0=1.2 D_MAX=1.5}
  ALPHA=27 PHI=0.0 THETA=-1.5708 PSI=-2.35619
...

# Transform the symmetry functions with a switching function
tfcc: LESS_THAN ARG=fcc SWITCH={SMAP R_0=0.5 A=8 B=8}

# Now compute the center of the solid like region
center: CENTER ATOMS=1-96000 WEIGHTS=tfcc

# This determines the positions of the atoms of interest relative to the center of the solid region
dens_dist: DISTANCES ORIGIN=center ATOMS=1-96000 COMPONENTS
# This computes the numerator in the expression above for the phase field
dens_numer: KDE ...
   VOLUMES=tfcc ARG=dens_dist.x,dens_dist.y,dens_dist.z
   GRID_BIN=80,80,80 BANDWIDTH=1.0,1.0,1.0
...
# This computes the denominator
dens_denom: KDE ...
   ARG=dens_dist.x,dens_dist.y,dens_dist.z
   GRID_BIN=80,80,80 BANDWIDTH=1.0,1.0,1.0
...
# This computes the final phase field
dens: CUSTOM ARG=dens_numer,dens_denom FUNC=x/y PERIODIC=NO

# Find the isocontour
cont: FIND_CONTOUR ARG=dens CONTOUR=0.5
# Use the special method for outputting the contour to a file
DUMPCONTOUR ARG=cont FILE=surface.xyz
```

Notice that with this method you have to use the [DUMPCONTOUR](DUMPCONTOUR.md) action to output any contour you find to a file.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace contour {

PLUMED_REGISTER_ACTION(FindContour,"FIND_CONTOUR")

void FindContour::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords( keys );
  ActionWithValue::useCustomisableComponents(keys);
  keys.addInputKeyword("compulsory","ARG","grid","the labels of the grid in which the contour will be found");
  ContourFindingObject<gridtools::EvaluateGridFunction>::registerKeywords( keys );
  keys.add("compulsory","BUFFER","0","number of buffer grid points around location where grid was found on last step.  If this is zero the full grid is calculated on each step");
  keys.addDOI("10.1039/B805786A");
  keys.addDOI("10.1021/jp909219k");
  keys.addDOI("10.1088/1361-648X/aa893d");
  keys.addDOI("10.1063/1.5134461");
  PTM::registerKeywords( keys );
}

FindContour::FindContour(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  firststep(true),
  taskmanager(this) {
  parse("BUFFER",gbuffer);
  if( gbuffer>0 ) {
    log.printf("  after first step a subset of only %u grid points around where the countour was found will be checked\n",gbuffer);
  }

  gridtools::ActionWithGrid* ag=dynamic_cast<gridtools::ActionWithGrid*>( getPntrToArgument(0)->getPntrToAction() );
  std::vector<std::string> argn( ag->getGridCoordinateNames() );

  std::vector<std::size_t> shape(1);
  shape[0]=0;
  for(unsigned i=0; i<argn.size(); ++i ) {
    addComponent( argn[i], shape );
    componentIsNotPeriodic( argn[i] );
  }
  function::FunctionOptions options;
  ContourFindingObject<gridtools::EvaluateGridFunction>::read( taskmanager.getActionInput(), this, options );
  log.printf("  calculating dividing surface along which function equals %f \n", taskmanager.getActionInput().contour);
}

std::string FindContour::getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const {
  return "a vector of coordinates for the contour along the " + cname + " direction";
}

void FindContour::calculate() {
  if( firststep ) {
    std::vector<std::size_t> shape(1);
    shape[0] = getPntrToArgument(0)->getRank()*getPntrToArgument(0)->getNumberOfValues();
    for(unsigned i=0; i<getNumberOfComponents(); ++i) {
      getPntrToComponent(i)->setShape( shape );
    }
    active_cells.resize( shape[0] );
    taskmanager.setupParallelTaskManager( 0, 0 );
    firststep = false;
  }
  taskmanager.runAllTasks();
}

unsigned FindContour::getNumberOfDerivatives() {
  return 0;
}

void FindContour::getNumberOfTasks( unsigned& ntasks ) {
  ntasks = active_cells.size();

  Value* gval=getPntrToArgument(0);
  unsigned npoints = gval->getNumberOfValues();
  std::vector<unsigned> ind( gval->getRank() );
  std::vector<unsigned> ones( gval->getRank(), 1 );
  std::vector<std::size_t> nbin( getInputGridObject().getNbin( false ) );
  unsigned num_neighbours;
  std::vector<unsigned> neighbours;

  std::fill( active_cells.begin(), active_cells.end(), 0 );
  for(unsigned i=0; i<npoints; ++i) {
    // Get the index of the current grid point
    getInputGridObject().getIndices( i, ind );
    getInputGridObject().getNeighbors( ind, ones, num_neighbours, neighbours );
    // Get the value of a point on the grid
    double val1=gval->get( i ) - taskmanager.getActionInput().contour;
    bool edge=false;
    for(unsigned j=0; j<gval->getRank(); ++j) {
      // Make sure we don't search at the edge of the grid
      if( !getInputGridObject().isPeriodic(j) && (ind[j]+1)==nbin[j] ) {
        continue;
      } else if( (ind[j]+1)==nbin[j] ) {
        edge=true;
        ind[j]=0;
      } else {
        ind[j]+=1;
      }
      double val2=gval->get( getInputGridObject().getIndex(ind) ) - taskmanager.getActionInput().contour;
      if( val1*val2<0 ) {
        active_cells[gval->getRank()*i + j] = 1;
      }
      if( getInputGridObject().isPeriodic(j) && edge ) {
        edge=false;
        ind[j]=nbin[j]-1;
      } else {
        ind[j]-=1;
      }
    }
  }
}

int FindContour::checkTaskIsActive( const unsigned& taskno ) const {
  if( active_cells[taskno]>0 ) {
    return 1;
  }
  return -1;
}

void FindContour::getInputData( std::vector<double>& inputdata ) const {
  std::size_t rank = getPntrToArgument(0)->getRank();
  std::size_t ndata = getInputGridObject().getNumberOfPoints()*rank;
  if( inputdata.size()!=2*rank*ndata ) {
    inputdata.resize( 2*rank*ndata );
  }
  std::vector<double> point( getPntrToArgument(0)->getRank() );
  for(unsigned i=0; i<getInputGridObject().getNumberOfPoints(); ++i) {
    getInputGridObject().getGridPointCoordinates( i, point );
    for(unsigned j=0; j<rank; ++j) {
      for(unsigned k=0; k<rank; ++k) {
        inputdata[i*rank*2*rank + j*2*rank + k] = point[k];
        inputdata[i*rank*2*rank + j*2*rank + rank + k] = 0;
      }
      inputdata[i*rank*2*rank + j*2*rank + rank + j] = 0.999999999*getInputGridObject().getGridSpacing()[j];
    }
  }
}

void FindContour::performTask( std::size_t task_index,
                               const ContourFindingObject<gridtools::EvaluateGridFunction>& actiondata,
                               ParallelActionsInput& input,
                               ParallelActionsOutput& output ) {

  std::size_t rank = actiondata.function.getGridObject().getDimension();
  std::vector<double> direction( rank ), point( rank );
  for(unsigned i=0; i<rank; ++i) {
    point[i] = input.inputdata[ 2*rank*task_index + i];
    direction[i] = input.inputdata[ 2*rank*task_index + rank + i];
  }
  // Now find the contour
  ContourFindingObject<gridtools::EvaluateGridFunction>::findContour( actiondata, direction, point );
  for(unsigned i=0; i<rank; ++i) {
    output.values[i] = point[i];
  }
}

}
}
