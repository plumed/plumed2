/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016 The plumed team
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

//+PLUMEDOC GRIDANALYSIS FIND_CONTOUR_SURFACE
/*
Find an isocontour by searching along either the x, y or direction.

\par Examples

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
public:
  static void registerKeywords( Keywords& keys );
  explicit FindContourSurface(const ActionOptions&ao);
  unsigned getNumberOfQuantities() const { return 2; }
  bool checkAllActive() const { return gbuffer==0; }
  void clearAverage();
  void prepareForAveraging();
  void compute( const unsigned& current, MultiValue& myvals ) const ;
  void finishAveraging();
};

PLUMED_REGISTER_ACTION(FindContourSurface,"FIND_CONTOUR_SURFACE")

void FindContourSurface::registerKeywords( Keywords& keys ){
  ContourFindingBase::registerKeywords( keys );
  keys.add("compulsory","SEARCHDIR","In which directions do you wish to search for the contour.");
  keys.add("compulsory","BUFFER","0","number of buffer grid points around location where grid was found on last step.  If this is zero the full grid is calculated on each step");
  keys.remove("FILE"); keys.remove("UNITS"); keys.remove("PRECISION");
}

FindContourSurface::FindContourSurface(const ActionOptions&ao):
Action(ao),
ContourFindingBase(ao),
firsttime(true),
ones(ingrid->getDimension(),1)
{
  if( ingrid->getDimension()<2 ) error("cannot find dividing surface if input grid is one dimensional");

  std::string dir; parse("SEARCHDIR",dir); parse("BUFFER",gbuffer);
  log.printf("  calculating location of contour on %d dimensional grid \n", ingrid->getDimension()-1 );
  if( gbuffer>0 ) log.printf("  after first step a subset of only %u grid points around where the countour was found will be checked\n",gbuffer);
  checkRead();

  unsigned n=0; gdirs.resize( ingrid->getDimension()-1 );
  for(unsigned i=0;i<ingrid->getDimension();++i){
      if( ingrid->getComponentName(i)==dir ){
          dir_n=i; 
      } else {
          if( n==gdirs.size() ) error("could not find " + dir + " direction in input grid");
          gdirs[n]=i; n++; 
      }
  }
  if( n!=(ingrid->getDimension()-1) ) error("output of grid is not understood");

  // Create the input from the old string
  std::string vstring = "COMPONENTS=" + getLabel() + " COORDINATES=" + ingrid->getComponentName( gdirs[0] );
  for(unsigned i=1;i<gdirs.size();++i) vstring += "," + ingrid->getComponentName( gdirs[i] );
  vstring += " PBC=";
  if( ingrid->isPeriodic(gdirs[0]) ) vstring+="T";
  else vstring+="F";
  for(unsigned i=1;i<gdirs.size();++i){
      if( ingrid->isPeriodic(gdirs[i]) ) vstring+=",T"; else vstring+=",F";
  }
  createGrid( "grid", vstring ); mygrid->setNoDerivatives(); 
  setAveragingAction( mygrid, true );
}

void FindContourSurface::clearAverage(){
  // Set the boundaries of the output grid
  std::vector<double> fspacing; std::vector<unsigned> snbins( ingrid->getDimension()-1 );
  std::vector<std::string> smin( ingrid->getDimension()-1 ), smax( ingrid->getDimension()-1 );
  for(unsigned i=0;i<gdirs.size();++i){
     smin[i]=ingrid->getMin()[gdirs[i]]; smax[i]=ingrid->getMax()[gdirs[i]];
     snbins[i]=ingrid->getNbin()[gdirs[i]];
  }   
  mygrid->setBounds( smin, smax, snbins, fspacing); resizeFunctions();
  ActionWithAveraging::clearAverage();
}

void FindContourSurface::prepareForAveraging(){
  // Create a task list if first time
  if( firsttime ){
      std::vector<unsigned> find( ingrid->getDimension() );
      std::vector<unsigned> ind( mygrid->getDimension() );
      for(unsigned i=0;i<mygrid->getNumberOfPoints();++i){
          find.assign( find.size(), 0 ); mygrid->getIndices( i, ind );
          for(unsigned j=0;j<gdirs.size();++j) find[gdirs[j]]=ind[j];
          // Current will be set equal to the start point for this grid index
          addTaskToList( ingrid->getIndex(find) );
      }
      // And prepare the task list
      deactivateAllTasks();
      for(unsigned i=0;i<getFullNumberOfTasks();++i) taskFlags[i]=1;
      lockContributors();
      // Set the direction in which to look for the contour
      direction.resize( ingrid->getDimension(), 0 );
      direction[dir_n] = 0.999999999*ingrid->getGridSpacing()[dir_n];
  }
}

void FindContourSurface::finishAveraging(){
  ContourFindingBase::finishAveraging();
  // And update the list of active grid points
  if( gbuffer>0 ){
      std::vector<double> dx( ingrid->getGridSpacing() );
      std::vector<double> point( ingrid->getDimension() );
      std::vector<double> lpoint( mygrid->getDimension() );
      std::vector<unsigned> neighbours; unsigned num_neighbours;
      std::vector<unsigned> ugrid_indices( ingrid->getDimension() );
      std::vector<bool> active( ingrid->getNumberOfPoints(), false ); 
      std::vector<unsigned> gbuffer_vec( ingrid->getDimension(), gbuffer );
      for(unsigned i=0;i<mygrid->getNumberOfPoints();++i){
          // Retrieve the coordinates of this grid point
          mygrid->getGridPointCoordinates( i, lpoint );
          point[dir_n] = mygrid->getGridElement( i, 0 );
          // 0.5*dx added here to prevent problems with flooring of grid points
          for(unsigned j=0;j<gdirs.size();++j) point[gdirs[j]]=lpoint[j] + 0.5*dx[gdirs[j]];
          ingrid->getIndices( point, ugrid_indices );
          // Now activate buffer region
          ingrid->getNeighbors( ugrid_indices , gbuffer_vec, num_neighbours, neighbours );
          for(unsigned n=0;n<num_neighbours;++n) active[ neighbours[n] ]=true; 
      }
      ingrid->activateThesePoints( active );
  }
  firsttime=false;
}

void FindContourSurface::compute( const unsigned& current, MultiValue& myvals ) const {
  std::vector<unsigned> neighbours; unsigned num_neighbours; unsigned nfound=0; double minv=0, minp;
  std::vector<unsigned> bins_n( ingrid->getNbin() ); unsigned shiftn=current;
  std::vector<unsigned> ind( ingrid->getDimension() ); std::vector<double> point( ingrid->getDimension() );
#ifndef DNDEBUG 
  std::vector<unsigned> oind( mygrid->getDimension() ); mygrid->getIndices( current, oind ); 
#endif
  for(unsigned i=0;i<bins_n[dir_n];++i){
#ifndef DNDEBUG
  std::vector<unsigned> base_ind( ingrid->getDimension() ); ingrid->getIndices( shiftn, base_ind );
  for(unsigned j=0;j<gdirs.size();++j) plumed_dbg_assert( base_ind[gdirs[j]]==oind[j] );
#endif
     // Ensure inactive grid points are ignored
     if( ingrid->inactive( shiftn ) ){ shiftn += ingrid->getStride()[dir_n]; continue; }
     // Get the index of the current grid point
     ingrid->getIndices( shiftn, ind );
     // Exit if we are at the edge of the grid
     if( !ingrid->isPeriodic(dir_n) && (ind[dir_n]+1)==bins_n[dir_n] ){
        shiftn += ingrid->getStride()[dir_n]; continue;
     }     

     // Ensure points with inactive neighbours are ignored
     ingrid->getNeighbors( ind, ones, num_neighbours, neighbours );
     bool cycle=false;
     for(unsigned j=0;j<num_neighbours;++j){
         if( ingrid->inactive( neighbours[j]) ){ cycle=true; break; }
     }
     if( cycle ){ shiftn += ingrid->getStride()[dir_n]; continue; }

     // Now get the function value at two points
     double val1=getFunctionValue( shiftn ) - contour; double val2;
     if( (ind[dir_n]+1)==bins_n[dir_n] ) val2 = getFunctionValue( current ) - contour;
     else val2=getFunctionValue( shiftn + ingrid->getStride()[dir_n] ) - contour;

     // Check if the minimum is bracketed 
     if( val1*val2<0 ){ 
         ingrid->getGridPointCoordinates( shiftn, point ); findContour( direction, point ); 
         minp=point[dir_n]; nfound++; break;
     }


     // This moves us on to the next point
     shiftn += ingrid->getStride()[dir_n];
  }
  if( nfound==0 ){
     std::string num; Tools::convert( getStep(), num );
     error("On step " + num + " failed to find required grid point"); 
  }
  myvals.setValue( 1, minp ); 
}

}
}
