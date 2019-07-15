/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#include "core/Atoms.h"

//+PLUMEDOC GRIDANALYSIS FIND_CONTOUR_SURFACE
/*
Find an isocontour by searching along either the x, y or z direction.

As discussed in the part of the manual on \ref Analysis PLUMED contains a number of tools that allow you to calculate
a function on a grid.  The function on this grid might be a \ref HISTOGRAM as a function of a few collective variables
or it might be a phase field that has been calculated using \ref MULTICOLVARDENS.  If this function has one or two input
arguments it is relatively straightforward to plot the function.  If by contrast the data has a three dimensions it can be
difficult to visualize.

This action provides one tool for visualizing these functions.  It can be used to search for a set of points on a contour
where the function takes a particular value.  In other words, for the function \f$f(x,y,z)\f$ this action would find a set
of points \f$\{x_c,y_c,z_c\}\f$ that have:

\f[
f(x_c,y_c,z_c) - c = 0
\f]

where \f$c\f$ is some constant value that is specified by the user.  The points on this contour are find by searching along lines
that run parallel to the \f$x\f$, \f$y\f$ or \f$z\f$ axis of the simulation cell.  The result is, therefore, a two dimensional
function evaluated on a grid that gives us the height of the interface as a function of two coordinates.

It is important to note that this action can only be used to detect contours in three dimensional functions.  In addition, this action will fail to
find the full set of contour  points if the contour does not have the same topology as an infinite plane.  If you are uncertain that the isocontours in your
function have the appropriate topology you should use \ref FIND_CONTOUR in place of \ref FIND_CONTOUR_SURFACE.


\par Examples

The input shown below was used to analyze the results from a simulation of an interface between solid and molten Lennard Jones.  The interface between
the solid and the liquid was set up in the plane perpendicular to the \f$z\f$ direction of the simulation cell.   The input below calculates something
akin to a Willard-Chandler dividing surface \cite wcsurface between the solid phase and the liquid phase.  There are two of these interfaces within the
simulation box because of the periodic boundary conditions but we were able to determine that one of these two surfaces lies in a particular part of the
simulation box.  The input below detects the height profile of one of these two interfaces.  It does so by computing a phase field average of the
\ref FCCUBIC symmetry function using the \ref MULTICOLVARDENS action.  Notice that we use the fact that we know roughly where the interface is when specifying
how this phase field is to be calculated and specify the region over the \f$z\f$-axis in which we are going to search for the phase field in the line defining
the \ref MULTICOLVARDENS.  Once we have calculated the phase field we search for contour points on the lines that run parallel to the \f$z\f$-direction of the cell
box using the FIND_CONTOUR_SURFACE command.  The final result is a \f$14 \times 14\f$ grid of values for the height of the interface as a function of the \f$(x,y)\f$
position.  This grid is then output to a file called contour2.dat.

Notice that the commands below calculate the instantaneous position of the surface separating the solid and liquid and that as such the accumulated average is cleared
on every step.

\plumedfile
UNITS NATURAL
FCCUBIC ...
  SPECIES=1-96000 SWITCH={CUBIC D_0=1.2 D_MAX=1.5}
  ALPHA=27 PHI=0.0 THETA=-1.5708 PSI=-2.35619 LABEL=fcc
... FCCUBIC

dens2: MULTICOLVARDENS DATA=fcc ORIGIN=1 DIR=xyz NBINS=14,14,50 ZREDUCED ZLOWER=6.0 ZUPPER=11.0 BANDWIDTH=1.0,1.0,1.0 CLEAR=1

ss2: FIND_CONTOUR_SURFACE GRID=dens2 CONTOUR=0.42 SEARCHDIR=z STRIDE=1 CLEAR=1
DUMPGRID GRID=ss2 FILE=contour2.dat FMT=%8.4f STRIDE=1
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class FindContourSurface : public ContourFindingBase {
private:
  bool firsttime;
  unsigned dir_n;
  unsigned gbuffer;
  std::vector<unsigned> ones;
  std::vector<unsigned> gdirs;
  std::vector<double> direction;
  std::vector<std::string> gnames;
  GridCoordinatesObject gridcoords;
public:
  static void registerKeywords( Keywords& keys );
  explicit FindContourSurface(const ActionOptions&ao);
  void finishOutputSetup();
  void getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
                             std::vector<std::string>& max, std::vector<unsigned>& nbin,
                             std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const ;
  void getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const ;
  void getGridPointAsCoordinate( const unsigned& ind, const bool& setlength, std::vector<double>& coords ) const ;
  unsigned getNumberOfDerivatives() const ;
  void performTask( const unsigned& current, MultiValue& myvals ) const ;
  void gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                          const unsigned& bufstart, std::vector<double>& buffer ) const ;
  void jobsAfterLoop();
};

PLUMED_REGISTER_ACTION(FindContourSurface,"FIND_CONTOUR_SURFACE")

void FindContourSurface::registerKeywords( Keywords& keys ) {
  ContourFindingBase::registerKeywords( keys );
  keys.add("compulsory","SEARCHDIR","In which directions do you wish to search for the contour.");
  keys.add("compulsory","BUFFER","0","number of buffer grid points around location where grid was found on last step.  If this is zero the full grid is calculated on each step");
}

FindContourSurface::FindContourSurface(const ActionOptions&ao):
  Action(ao),
  ContourFindingBase(ao),
  firsttime(true),
  ones(getPntrToArgument(0)->getRank(),1)
{
  if( getPntrToArgument(0)->getRank()<2 ) error("cannot find dividing surface if input grid is one dimensional");

  std::string dir; parse("SEARCHDIR",dir); parse("BUFFER",gbuffer);
  log.printf("  calculating location of contour on %d dimensional grid \n", getPntrToArgument(0)->getRank()-1 );
  if( gbuffer>0 ) log.printf("  after first step a subset of only %u grid points around where the countour was found will be checked\n",gbuffer);
  checkRead();

  unsigned n=0; gdirs.resize( getPntrToArgument(0)->getRank()-1 ); gnames.resize( getPntrToArgument(0)->getRank()-1 );
  Value* gval=getPntrToArgument(0); std::vector<unsigned> nbin( gval->getRank() ); std::string gtype;
  std::vector<double> spacing( gval->getRank() ); std::vector<bool> pbc( gval->getRank() );
  std::vector<std::string> argn( gval->getRank() ), min( gval->getRank() ), max( gval->getRank() );
  gval->getPntrToAction()->getInfoForGridHeader( gtype, argn, min, max, nbin, spacing, pbc, false );
  if( gtype=="fibonacci") error("cannot search for contours in fibonacci grids");
  for(unsigned i=0; i<getPntrToArgument(0)->getRank(); ++i) {
    if( argn[i]==dir ) {
      dir_n=i;
    } else {
      if( n==gdirs.size() ) error("could not find " + dir + " direction in input grid");
      gdirs[n]=i; gnames[n]=argn[i]; n++;
    }
  }
  if( n!=(getPntrToArgument(0)->getRank()-1) ) error("output of grid is not understood");

  std::vector<bool> ipbc( gridobject.getDimension()-1 );
  for(unsigned i=0; i<gdirs.size(); ++i) ipbc[i] = gridobject.isPeriodic(gdirs[i]);
  gridcoords.setup( "flat", ipbc, 0, 0.0 );

  // Now add a value
  std::vector<unsigned> shape( gridobject.getDimension()-1 ); 
  addValueWithDerivatives( shape ); setNotPeriodic();
  getPntrToOutput(0)->alwaysStoreValues();
}

void FindContourSurface::finishOutputSetup() {
  std::vector<double> fspacing; std::vector<unsigned> snbins( gridcoords.getDimension() );
  std::vector<std::string> smin( gridcoords.getDimension() ), smax( gridcoords.getDimension() );
  for(unsigned i=0; i<gdirs.size(); ++i) {
    smin[i]=gridobject.getMin()[gdirs[i]]; smax[i]=gridobject.getMax()[gdirs[i]];
    snbins[i]=gridobject.getNbin(false)[gdirs[i]];
  }
  gridcoords.setBounds( smin, smax, snbins, fspacing );
  getPntrToOutput(0)->setShape( gridcoords.getNbin(true) );

  std::vector<unsigned> find( gridcoords.getDimension() );
  std::vector<unsigned> ind( gridcoords.getDimension() );
  for(unsigned i=0; i<gridcoords.getNumberOfPoints(); ++i) {
    find.assign( find.size(), 0 ); gridcoords.getIndices( i, ind );
    for(unsigned j=0; j<gdirs.size(); ++j) find[gdirs[j]]=ind[j];
    // Current will be set equal to the start point for this grid index
    addTaskToList( gridcoords.getIndex(find) );
  }

  // Set the direction in which to look for the contour
  direction.resize( gridobject.getDimension(), 0 );
  direction[dir_n] = 0.999999999*gridobject.getGridSpacing()[dir_n];
}

void FindContourSurface::getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
    std::vector<std::string>& max, std::vector<unsigned>& nbin,
    std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const {
  bool isdists=dumpcube; double units=1.0; gtype="flat";
  for(unsigned i=0; i<getPntrToOutput(0)->getRank(); ++i) {
    if( gnames[i].find(".")==std::string::npos ) { isdists=false; break; }
    std::size_t dot = gnames[i].find("."); std::string name = gnames[i].substr(dot+1);
    if( name!="x" && name!="y" && name!="z" ) { isdists=false; break; }
  }
  if( isdists ) {
    if( plumed.getAtoms().usingNaturalUnits() ) units = 1.0/0.5292;
    else units = plumed.getAtoms().getUnits().getLength()/.05929;
  }
  std::vector<unsigned> nn( gridcoords.getNbin( false ) );
  for(unsigned i=0; i<getPntrToOutput(0)->getRank(); ++i) {
    argn[i] = gnames[i];
    double gmin, gmax;
    if( gridobject.getMin().size()>0 ) {
      Tools::convert( gridcoords.getMin()[i], gmin ); Tools::convert( gmin*units, min[i] );
      Tools::convert( gridcoords.getMax()[i], gmax ); Tools::convert( gmax*units, max[i] );
    }
    if( gridcoords.getGridSpacing().size()>0 ) spacing[i]=units*gridcoords.getGridSpacing()[i];
    nbin[i]=nn[i]; pbc[i]=gridcoords.isPeriodic(i);
  }
}

void FindContourSurface::getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const {
  gridcoords.getGridPointCoordinates( ind, indices, coords );
}

void FindContourSurface::getGridPointAsCoordinate( const unsigned& ind, const bool& setlength, std::vector<double>& coords ) const {
  if( coords.size()!=gridobject.getDimension() ) coords.resize( gridobject.getDimension() );
  std::vector<double> point( gridcoords.getDimension() ); gridcoords.getGridPointCoordinates( ind, point );
  for(unsigned i=0; i<gdirs.size(); ++i) coords[gdirs[i]]=point[i]; coords[dir_n] = getPntrToOutput(0)->get(ind);
}

unsigned FindContourSurface::getNumberOfDerivatives() const {
  return gridcoords.getDimension();
}

void FindContourSurface::jobsAfterLoop() {
  // And update the list of active grid points
  // if( gbuffer>0 ) {
  //   std::vector<double> dx( ingrid->getGridSpacing() );
  //   std::vector<double> point( ingrid->getDimension() );
  //   std::vector<double> lpoint( mygrid->getDimension() );
  //   std::vector<unsigned> neighbours; unsigned num_neighbours;
  //   std::vector<unsigned> ugrid_indices( ingrid->getDimension() );
  //   std::vector<bool> active( ingrid->getNumberOfPoints(), false );
  //   std::vector<unsigned> gbuffer_vec( ingrid->getDimension(), gbuffer );
  //   for(unsigned i=0; i<mygrid->getNumberOfPoints(); ++i) {
  //     // Retrieve the coordinates of this grid point
  //     mygrid->getGridPointCoordinates( i, lpoint );
  //     point[dir_n] = mygrid->getGridElement( i, 0 );
  //     // 0.5*dx added here to prevent problems with flooring of grid points
  //     for(unsigned j=0; j<gdirs.size(); ++j) point[gdirs[j]]=lpoint[j] + 0.5*dx[gdirs[j]];
  //     ingrid->getIndices( point, ugrid_indices );
  //     // Now activate buffer region
  //     ingrid->getNeighbors( ugrid_indices, gbuffer_vec, num_neighbours, neighbours );
  //     for(unsigned n=0; n<num_neighbours; ++n) active[ neighbours[n] ]=true;
  //   }
  //   ingrid->activateThesePoints( active );
  // }
}

void FindContourSurface::performTask( const unsigned& current, MultiValue& myvals ) const {
  std::vector<unsigned> neighbours; unsigned num_neighbours; unsigned nfound=0; double minv=0, minp;
  std::vector<unsigned> bins_n( gridobject.getNbin(false) ); unsigned shiftn=current;
  std::vector<unsigned> ind( gridobject.getDimension() ); std::vector<double> point( gridobject.getDimension() );
#ifndef DNDEBUG
  std::vector<unsigned> oind( gridcoords.getDimension() ); gridcoords.getIndices( current, oind );
#endif
  for(unsigned i=0; i<bins_n[dir_n]; ++i) {
#ifndef DNDEBUG
    std::vector<unsigned> base_ind( gridobject.getDimension() ); gridobject.getIndices( shiftn, base_ind );
    for(unsigned j=0; j<gdirs.size(); ++j) plumed_dbg_assert( base_ind[gdirs[j]]==oind[j] );
#endif
    // Ensure inactive grid points are ignored
    //if( ingrid->inactive( shiftn ) ) { shiftn += ingrid->getStride()[dir_n]; continue; }
    // Get the index of the current grid point
    gridobject.getIndices( shiftn, ind );
    // Exit if we are at the edge of the grid
    if( !gridobject.isPeriodic(dir_n) && (ind[dir_n]+1)==bins_n[dir_n] ) {
      shiftn += gridobject.getStride()[dir_n]; continue;
    }

    // Ensure points with inactive neighbours are ignored
    gridobject.getNeighbors( ind, ones, num_neighbours, neighbours );
    bool cycle=false;
    // for(unsigned j=0; j<num_neighbours; ++j) {
    //   if( ingrid->inactive( neighbours[j]) ) { cycle=true; break; }
    // }
    // if( cycle ) { shiftn += ingrid->getStride()[dir_n]; continue; }

    // Now get the function value at two points
    double val1=getFunctionValue( shiftn ) - contour; double val2;
    if( (ind[dir_n]+1)==bins_n[dir_n] ) val2 = getFunctionValue( current ) - contour;
    else val2=getFunctionValue( shiftn + gridobject.getStride()[dir_n] ) - contour;

    // Check if the minimum is bracketed
    if( val1*val2<0 ) {
      gridobject.getGridPointCoordinates( shiftn, point ); findContour( direction, point );
      minp=point[dir_n]; nfound++; break;
    }


    // This moves us on to the next point
    shiftn += gridobject.getStride()[dir_n];
  }
  if( nfound==0 ) {
    std::string num; Tools::convert( getStep(), num );
    error("On step " + num + " failed to find required grid point");
  }
  myvals.setValue( getPntrToOutput(0)->getPositionInStream(), minp );
}

void FindContourSurface::gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                                            const unsigned& bufstart, std::vector<double>& buffer ) const {
  plumed_dbg_assert( valindex==0 ); unsigned istart = bufstart + (1+getNumberOfDerivatives())*code;
  unsigned valout = getPntrToOutput(0)->getPositionInStream(); buffer[istart] += myvals.get( valout );
}

}
}
