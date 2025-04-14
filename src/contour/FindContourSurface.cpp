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
#include "ContourFindingBase.h"
#include "core/PlumedMain.h"

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
dens2_numer: KDE VOLUMES=fcc_n ARG=dens2_dist.x,dens2_dist.y,dens2_dist.z GRID_BIN=14,14,50 GRID_MIN=auto,auto,6.0 GRID_MAX=auto,auto,11.0 BANDWIDTH=1.0,1.0,1.0
# This computes the denominator
dens2_denom: KDE ARG=dens2_dist.x,dens2_dist.y,dens2_dist.z GRID_BIN=14,14,50 GRID_MIN=auto,auto,6.0 GRID_MAX=auto,auto,11.0 BANDWIDTH=1.0,1.0,1.0
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

class FindContourSurface : public ContourFindingBase {
private:
  unsigned dir_n;
  std::vector<unsigned> ones;
  std::vector<unsigned> gdirs;
  std::vector<double> direction;
  std::vector<std::string> gnames;
  gridtools::GridCoordinatesObject gridcoords;
public:
  static void registerKeywords( Keywords& keys );
  explicit FindContourSurface(const ActionOptions&ao);
  void setupValuesOnFirstStep() override;
  unsigned getNumberOfDerivatives() override ;
  std::vector<std::string> getGridCoordinateNames() const override ;
  const gridtools::GridCoordinatesObject& getGridCoordinatesObject() const override ;
  void performTask( const unsigned& current, MultiValue& myvals ) const override;
  void gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                          const unsigned& bufstart, std::vector<double>& buffer ) const override ;
};

PLUMED_REGISTER_ACTION(FindContourSurface,"FIND_CONTOUR_SURFACE")

void FindContourSurface::registerKeywords( Keywords& keys ) {
  ContourFindingBase::registerKeywords( keys );
  keys.add("compulsory","SEARCHDIR","In which directions do you wish to search for the contour.");
  keys.setValueDescription("grid","a grid containing the location of the points in the Willard-Chandler surface along the chosen direction");
  keys.addDOI("10.1088/1361-648X/aa893d");
}

FindContourSurface::FindContourSurface(const ActionOptions&ao):
  Action(ao),
  ContourFindingBase(ao),
  ones(getPntrToArgument(0)->getRank(),1) {
  if( getPntrToArgument(0)->getRank()<2 ) {
    error("cannot find dividing surface if input grid is one dimensional");
  }

  std::string dir;
  parse("SEARCHDIR",dir);
  log.printf("  calculating location of contour on %d dimensional grid \n", getPntrToArgument(0)->getRank()-1 );
  checkRead();

  Value* gval=getPntrToArgument(0);
  unsigned n=0;
  gdirs.resize( gval->getRank()-1 );
  gnames.resize( getPntrToArgument(0)->getRank()-1 );

  gridtools::ActionWithGrid* ag=dynamic_cast<gridtools::ActionWithGrid*>( gval->getPntrToAction() );
  if( !ag ) {
    error("input argument must be a grid");
  }
  if( getInputGridObject().getGridType()=="fibonacci") {
    error("cannot search for contours in fibonacci grids");
  }
  std::vector<std::string> argn( ag->getGridCoordinateNames() );

  for(unsigned i=0; i<gval->getRank(); ++i) {
    if( argn[i]==dir ) {
      dir_n=i;
    } else {
      if( n==gdirs.size() ) {
        error("could not find " + dir + " direction in input grid");
      }
      gdirs[n]=i;
      gnames[n]=argn[i];
      n++;
    }
  }
  if( n!=(gval->getRank()-1) ) {
    error("output of grid is not understood");
  }

  std::vector<bool> ipbc( getInputGridObject().getDimension()-1 );
  for(unsigned i=0; i<gdirs.size(); ++i) {
    ipbc[i] = getInputGridObject().isPeriodic(gdirs[i]);
  }
  gridcoords.setup( "flat", ipbc, 0, 0.0 );

  // Now add a value
  std::vector<unsigned> shape( getInputGridObject().getDimension()-1 );
  addValueWithDerivatives( shape );
  setNotPeriodic();
  getPntrToComponent(0)->buildDataStore();
}

void FindContourSurface::setupValuesOnFirstStep() {
  std::vector<double> fspacing;
  std::vector<unsigned> snbins( gridcoords.getDimension() );
  std::vector<std::string> smin( gridcoords.getDimension() ), smax( gridcoords.getDimension() );
  for(unsigned i=0; i<gdirs.size(); ++i) {
    smin[i]=getInputGridObject().getMin()[gdirs[i]];
    smax[i]=getInputGridObject().getMax()[gdirs[i]];
    snbins[i]=getInputGridObject().getNbin(false)[gdirs[i]];
  }
  gridcoords.setBounds( smin, smax, snbins, fspacing );
  getPntrToComponent(0)->setShape( gridcoords.getNbin(true) );

  std::vector<unsigned> find( gridcoords.getDimension() );
  std::vector<unsigned> ind( gridcoords.getDimension() );
  for(unsigned i=0; i<gridcoords.getNumberOfPoints(); ++i) {
    find.assign( find.size(), 0 );
    gridcoords.getIndices( i, ind );
    for(unsigned j=0; j<gdirs.size(); ++j) {
      find[gdirs[j]]=ind[j];
    }
  }

  // Set the direction in which to look for the contour
  direction.resize( getInputGridObject().getDimension(), 0 );
  direction[dir_n] = 0.999999999*getInputGridObject().getGridSpacing()[dir_n];
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

void FindContourSurface::performTask( const unsigned& current, MultiValue& myvals ) const {
  std::vector<unsigned> neighbours;
  unsigned num_neighbours;
  unsigned nfound=0;
  double minv=0, minp;
  std::vector<unsigned> bins_n( getInputGridObject().getNbin(false) );
  unsigned shiftn=current;
  std::vector<unsigned> ind( getInputGridObject().getDimension() );
  std::vector<double> point( getInputGridObject().getDimension() );
#ifndef DNDEBUG
  std::vector<unsigned> oind( gridcoords.getDimension() );
  gridcoords.getIndices( current, oind );
#endif
  for(unsigned i=0; i<bins_n[dir_n]; ++i) {
#ifndef DNDEBUG
    std::vector<unsigned> base_ind( getInputGridObject().getDimension() );
    getInputGridObject().getIndices( shiftn, base_ind );
    for(unsigned j=0; j<gdirs.size(); ++j) {
      plumed_dbg_assert( base_ind[gdirs[j]]==oind[j] );
    }
#endif
    // Get the index of the current grid point
    getInputGridObject().getIndices( shiftn, ind );
    // Exit if we are at the edge of the grid
    if( !getInputGridObject().isPeriodic(dir_n) && (ind[dir_n]+1)==bins_n[dir_n] ) {
      shiftn += getInputGridObject().getStride()[dir_n];
      continue;
    }

    // Ensure points with inactive neighbours are ignored
    getInputGridObject().getNeighbors( ind, ones, num_neighbours, neighbours );

    // Now get the function value at two points
    double val1=getPntrToArgument(0)->get( shiftn ) - contour;
    double val2;
    if( (ind[dir_n]+1)==bins_n[dir_n] ) {
      val2 = getPntrToArgument(0)->get( current ) - contour;
    } else {
      val2=getPntrToArgument(0)->get( shiftn + getInputGridObject().getStride()[dir_n] ) - contour;
    }

    // Check if the minimum is bracketed
    if( val1*val2<0 ) {
      getInputGridObject().getGridPointCoordinates( shiftn, point );
      findContour( direction, point );
      minp=point[dir_n];
      nfound++;
      break;
    }

    // This moves us on to the next point
    shiftn += getInputGridObject().getStride()[dir_n];
  }
  if( nfound==0 ) {
    std::string num;
    Tools::convert( getStep(), num );
    error("On step " + num + " failed to find required grid point");
  }
  myvals.setValue( getConstPntrToComponent(0)->getPositionInStream(), minp );
}

void FindContourSurface::gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
    const unsigned& bufstart, std::vector<double>& buffer ) const {
  plumed_dbg_assert( valindex==0 );
  unsigned istart = bufstart + (1+gridcoords.getDimension())*code;
  unsigned valout = getConstPntrToComponent(0)->getPositionInStream();
  buffer[istart] += myvals.get( valout );
}

}
}
