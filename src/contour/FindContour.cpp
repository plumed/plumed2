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

As discussed in the part of the manual on \ref Analysis PLUMED contains a number of tools that allow you to calculate
a function on a grid.  The function on this grid might be a \ref HISTOGRAM as a function of a few collective variables
or it might be a phase field that has been calculated using \ref MULTICOLVARDENS.  If this function has one or two input
arguments it is relatively straightforward to plot the function.  If by contrast the data has a three or more dimensions
it can be difficult to visualize.

This action provides one tool for visualizing these functions.  It can be used to search for a set of points on a contour
where the function takes a particular values.  In other words, for the function \f$f(x,y)\f$ this action would find a set
of points \f$\{x_c,y_c\}\f$ that have:

\f[
f(x_c,y_c) - c = 0
\f]

where \f$c\f$ is some constant value that is specified by the user.  The points on this contour are detected using a variant
on the marching squares or marching cubes algorithm, which you can find information on here:

https://en.wikipedia.org/wiki/Marching_squares
https://en.wikipedia.org/wiki/Marching_cubes

As such, and unlike \ref FIND_CONTOUR_SURFACE or \ref FIND_SPHERICAL_CONTOUR, the function input to this action can have any dimension.
Furthermore, the topology of the contour will be determined by the algorithm and does not need to be specified by the user.

\par Examples

The input below allows you to calculate something akin to a Willard-Chandler dividing surface \cite wcsurface.
The simulation cell in this case contains a solid phase and a liquid phase.  The Willard-Chandler surface is the
surface that separates the parts of the box containing the solid from the parts containing the liquid.  To compute the position
of this surface  the \ref FCCUBIC symmetry function is calculated for each of the atoms in the system from on the geometry of the
atoms in the first coordination sphere of each of the atoms.  These quantities are then transformed using a switching function.
This procedure generates a single number for each atom in the system and this quantity has a value of one for atoms that are in
parts of the box that resemble the solid structure and zero for atoms that are in parts of the box that resemble the liquid.
The position of a virtual atom is then computed using \ref CENTER_OF_MULTICOLVAR and a phase field model is constructed using
\ref MULTICOLVARDENS.  These procedure ensures that we have a continuous function that gives a measure of the average degree of
solidness at each point in the simulation cell.  The Willard-Chandler dividing surface is calculated by finding a a set of points
at which the value of this phase field is equal to 0.5.  This set of points is output to file called mycontour.dat.  A new contour
is found on every single step for each frame that is read in.

\plumedfile
UNITS NATURAL
FCCUBIC ...
  SPECIES=1-96000 SWITCH={CUBIC D_0=1.2 D_MAX=1.5}
  ALPHA=27 PHI=0.0 THETA=-1.5708 PSI=-2.35619 LABEL=fcc
... FCCUBIC

tfcc: MTRANSFORM_MORE DATA=fcc LOWMEM SWITCH={SMAP R_0=0.5 A=8 B=8}
center: CENTER_OF_MULTICOLVAR DATA=tfcc

dens: MULTICOLVARDENS ...
  DATA=tfcc ORIGIN=center DIR=xyz
  NBINS=80,80,80 BANDWIDTH=1.0,1.0,1.0 STRIDE=1 CLEAR=1
...

FIND_CONTOUR GRID=dens CONTOUR=0.5 FILE=mycontour.xyz
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace contour {

PLUMED_REGISTER_ACTION(FindContour,"FIND_CONTOUR")

void FindContour::registerKeywords( Keywords& keys ) {
  ContourFindingBase::registerKeywords( keys );
  ActionWithValue::useCustomisableComponents(keys);
// We want a better way of doing this bit
  keys.add("compulsory","BUFFER","0","number of buffer grid points around location where grid was found on last step.  If this is zero the full grid is calculated on each step");
}

FindContour::FindContour(const ActionOptions&ao):
  Action(ao),
  ContourFindingBase(ao) {
  parse("BUFFER",gbuffer);
  if( gbuffer>0 ) {
    log.printf("  after first step a subset of only %u grid points around where the countour was found will be checked\n",gbuffer);
  }
  checkRead();

  gridtools::ActionWithGrid* ag=dynamic_cast<gridtools::ActionWithGrid*>( getPntrToArgument(0)->getPntrToAction() );
  std::vector<std::string> argn( ag->getGridCoordinateNames() );

  std::vector<unsigned> shape(1);
  shape[0]=0;
  for(unsigned i=0; i<argn.size(); ++i ) {
    addComponent( argn[i], shape );
    componentIsNotPeriodic( argn[i] );
    getPntrToComponent(i)->buildDataStore();
  }
  // Check for task reduction
  updateTaskListReductionStatus();
}

std::string FindContour::getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const {
  return "a vector of coordinates for the contour along the " + cname + " direction";
}

void FindContour::setupValuesOnFirstStep() {
  std::vector<unsigned> shape(1);
  shape[0] = getPntrToArgument(0)->getRank()*getPntrToArgument(0)->getNumberOfValues();
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    getPntrToComponent(i)->setShape( shape );
  }
  active_cells.resize( shape[0] );
}

unsigned FindContour::getNumberOfDerivatives() {
  return 0;
}

void FindContour::areAllTasksRequired( std::vector<ActionWithVector*>& task_reducing_actions ) {
  task_reducing_actions.push_back(this);
}

void FindContour::getNumberOfTasks( unsigned& ntasks ) {
  ntasks = active_cells.size();

  Value* gval=getPntrToArgument(0);
  unsigned npoints = gval->getNumberOfValues();
  std::vector<unsigned> ind( gval->getRank() );
  std::vector<unsigned> ones( gval->getRank(), 1 );
  std::vector<unsigned> nbin( getInputGridObject().getNbin( false ) );
  unsigned num_neighbours;
  std::vector<unsigned> neighbours;

  std::fill( active_cells.begin(), active_cells.end(), 0 );
  for(unsigned i=0; i<npoints; ++i) {
    // Get the index of the current grid point
    getInputGridObject().getIndices( i, ind );
    getInputGridObject().getNeighbors( ind, ones, num_neighbours, neighbours );
    // Get the value of a point on the grid
    double val1=gval->get( i ) - contour;
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
      double val2=gval->get( getInputGridObject().getIndex(ind) ) - contour;
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

int FindContour::checkTaskStatus( const unsigned& taskno, int& flag ) const {
  if( active_cells[taskno]>0 ) {
    return 1;
  }
  return 0;
}

void FindContour::performTask( const unsigned& current, MultiValue& myvals ) const {
  // Retrieve the initial grid point coordinates
  unsigned gpoint = std::floor( current / getPntrToArgument(0)->getRank() );
  std::vector<double> point( getPntrToArgument(0)->getRank() );
  getInputGridObject().getGridPointCoordinates( gpoint, point );

  // Retrieve the direction we are searching for the contour
  unsigned gdir = current%(getPntrToArgument(0)->getRank() );
  std::vector<double> direction( getPntrToArgument(0)->getRank(), 0 );
  direction[gdir] = 0.999999999*getInputGridObject().getGridSpacing()[gdir];

  // Now find the contour
  findContour( direction, point );
  // And transfer to the store data vessel
  for(unsigned i=0; i<getPntrToArgument(0)->getRank(); ++i) {
    myvals.setValue( getConstPntrToComponent(i)->getPositionInStream(), point[i] );
  }
}

}
}
