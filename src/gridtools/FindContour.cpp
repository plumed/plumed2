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
#include "core/ActionRegister.h"
#include "vesselbase/StoreDataVessel.h"
#include "ContourFindingBase.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

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
namespace gridtools {

class FindContour : public ContourFindingBase {
private:
  bool firsttime;
  unsigned gbuffer;
/// Stuff for output
  OFile of;
  double lenunit;
  std::string fmt_xyz;
/// The data is stored in a grid
  vesselbase::StoreDataVessel* mydata;
public:
  static void registerKeywords( Keywords& keys );
  explicit FindContour(const ActionOptions&ao);
  bool checkAllActive() const override { return gbuffer==0; }
  void prepareForAveraging() override;
  bool isPeriodic() override { return false; }
  unsigned getNumberOfQuantities() const override { return 1 + ingrid->getDimension(); }
  void compute( const unsigned& current, MultiValue& myvals ) const override;
  void finishAveraging() override;
};

PLUMED_REGISTER_ACTION(FindContour,"FIND_CONTOUR")

void FindContour::registerKeywords( Keywords& keys ) {
  ContourFindingBase::registerKeywords( keys );
// We want a better way of doing this bit
  keys.add("compulsory","BUFFER","0","number of buffer grid points around location where grid was found on last step.  If this is zero the full grid is calculated on each step");
  keys.add("compulsory","FILE","file on which to output coordinates");
  keys.add("compulsory","UNITS","PLUMED","the units in which to print out the coordinates. PLUMED means internal PLUMED units");
  keys.add("optional", "PRECISION","The number of digits in trajectory file");
}

FindContour::FindContour(const ActionOptions&ao):
  Action(ao),
  ContourFindingBase(ao),
  firsttime(true)
{

  parse("BUFFER",gbuffer);
  if( gbuffer>0 ) log.printf("  after first step a subset of only %u grid points around where the countour was found will be checked\n",gbuffer);

  std::string file; parse("FILE",file);
  if( file.length()==0 ) error("name out output file was not specified");
  std::string type=Tools::extension(file);
  log<<"  file name "<<file<<"\n";
  if(type!="xyz") error("can only print xyz file type with contour finding");

  fmt_xyz="%f";
  std::string precision; parse("PRECISION",precision);
  if(precision.length()>0) {
    int p; Tools::convert(precision,p);
    log<<"  with precision "<<p<<"\n";
    std::string a,b;
    Tools::convert(p+5,a);
    Tools::convert(p,b);
    fmt_xyz="%"+a+"."+b+"f";
  }
  std::string unitname; parse("UNITS",unitname);
  if(unitname!="PLUMED") {
    Units myunit; myunit.setLength(unitname);
    lenunit=plumed.getAtoms().getUnits().getLength()/myunit.getLength();
  }
  else lenunit=1.0;
  of.link(*this); of.open(file);
  checkRead(); mydata=buildDataStashes( NULL );
}

void FindContour::prepareForAveraging() {
  // Create a task list if first time
  if( firsttime ) {
    for(unsigned i=0; i<ingrid->getDimension()*ingrid->getNumberOfPoints(); ++i) addTaskToList( i );
  }
  firsttime=false; deactivateAllTasks();

  // We now need to identify the grid points that we need to search through
  std::vector<unsigned> nbin( ingrid->getNbin() );
  std::vector<unsigned> ind( ingrid->getDimension() );
  std::vector<unsigned> ones( ingrid->getDimension(), 1 );
  unsigned num_neighbours; std::vector<unsigned> neighbours;
  for(unsigned i=0; i<ingrid->getNumberOfPoints(); ++i) {
    // Ensure inactive grid points are ignored
    if( ingrid->inactive(i) ) continue;

    // Get the index of the current grid point
    ingrid->getIndices( i, ind );
    ingrid->getNeighbors( ind, ones, num_neighbours, neighbours );
    bool cycle=false;
    for(unsigned j=0; j<num_neighbours; ++j) {
      if( ingrid->inactive( neighbours[j]) ) { cycle=true; break; }
    }
    if( cycle ) continue;

    // Get the value of a point on the grid
    double val1=getFunctionValue( i ) - contour;
    bool edge=false;
    for(unsigned j=0; j<ingrid->getDimension(); ++j) {
      // Make sure we don't search at the edge of the grid
      if( !ingrid->isPeriodic(j) && (ind[j]+1)==nbin[j] ) continue;
      else if( (ind[j]+1)==nbin[j] ) { edge=true; ind[j]=0; }
      else ind[j]+=1;
      double val2=getFunctionValue( ind ) - contour;
      if( val1*val2<0 ) taskFlags[ ingrid->getDimension()*i + j ] = 1;
      if( ingrid->isPeriodic(j) && edge ) { edge=false; ind[j]=nbin[j]-1; }
      else ind[j]-=1;
    }
  }
  lockContributors();
}

void FindContour::compute( const unsigned& current, MultiValue& myvals ) const {
  // Retrieve the initial grid point coordinates
  unsigned gpoint = std::floor( current / ingrid->getDimension() );
  std::vector<double> point( ingrid->getDimension() );
  ingrid->getGridPointCoordinates( gpoint, point );

  // Retrieve the direction we are searching for the contour
  unsigned gdir = current%(ingrid->getDimension() );
  std::vector<double> direction( ingrid->getDimension(), 0 );
  direction[gdir] = 0.999999999*ingrid->getGridSpacing()[gdir];

  // Now find the contour
  findContour( direction, point );
  // And transfer to the store data vessel
  for(unsigned i=0; i<ingrid->getDimension(); ++i) myvals.setValue( 1+i, point[i] );
}

void FindContour::finishAveraging() {
  // And update the list of active grid points
  if( gbuffer>0 ) {
    std::vector<unsigned> neighbours; unsigned num_neighbours;
    std::vector<unsigned> ugrid_indices( ingrid->getDimension() );
    std::vector<bool> active( ingrid->getNumberOfPoints(), false );
    std::vector<unsigned> gbuffer_vec( ingrid->getDimension(), gbuffer );
    for(unsigned i=0; i<getCurrentNumberOfActiveTasks(); ++i) {
      // Get the point we are operating on
      unsigned ipoint = std::floor( getActiveTask(i) / ingrid->getDimension() );
      // Get the indices of this point
      ingrid->getIndices( ipoint, ugrid_indices );
      // Now activate buffer region
      ingrid->getNeighbors( ugrid_indices, gbuffer_vec, num_neighbours, neighbours );
      for(unsigned n=0; n<num_neighbours; ++n) active[ neighbours[n] ]=true;
    }
    ingrid->activateThesePoints( active );
  }
  std::vector<double> point( 1 + ingrid->getDimension() );
  of.printf("%u\n",mydata->getNumberOfStoredValues());
  of.printf("Points found on isocontour\n");
  for(unsigned i=0; i<mydata->getNumberOfStoredValues(); ++i) {
    mydata->retrieveSequentialValue( i, false, point ); of.printf("X");
    for(unsigned j=0; j<ingrid->getDimension(); ++j) of.printf( (" " + fmt_xyz).c_str(), lenunit*point[1+j] );
    of.printf("\n");
  }
}

}
}
