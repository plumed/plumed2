/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "GridFunction.h"

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
  GridFunction* outgrid;
public:
  static void registerKeywords( Keywords& keys );
  explicit FindContourSurface(const ActionOptions&ao);
  void findContourSurface();
  bool checkAllActive() const { return gbuffer==0; }
  void performTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(FindContourSurface,"FIND_CONTOUR_SURFACE")

void FindContourSurface::registerKeywords( Keywords& keys ){
  ContourFindingBase::registerKeywords( keys );
  keys.add("compulsory","SEARCHDIR","In which directions do you wish to search for the contour.");
  keys.add("compulsory","BUFFER","0","number of buffer grid points around location where grid was found on last step.  If this is zero the full grid is calculated on each step");
  keys.reset_style("FILE","optional");
}

FindContourSurface::FindContourSurface(const ActionOptions&ao):
Action(ao),
ContourFindingBase(ao),
firsttime(true),
ones(mygrid->getDimension(),1)
{
  if( mygrid->getDimension()<2 ) error("cannot find dividing surface if input grid is one dimensional");

  std::string dir; parse("SEARCHDIR",dir); parse("BUFFER",gbuffer);
  log.printf("  calculating location of contour on %d dimensional grid \n", mygrid->getDimension()-1 );
  if( gbuffer>0 ) log.printf("  after first step a subset of only %u grid points around where the countour was found will be checked\n",gbuffer);
  checkRead();

  unsigned n=0; gdirs.resize( mygrid->getDimension()-1 );
  for(unsigned i=0;i<mygrid->getDimension();++i){
      if( mygrid->getComponentName(i)==dir ){
          dir_n=i; 
      } else {
          if( n==gdirs.size() ) error("could not find " + dir + " direction in input grid");
          gdirs[n]=i; n++; 
      }
  }
  if( n!=(mygrid->getDimension()-1) ) error("output of grid is not understood");

  // Create the input from the old string
  std::string vstring = "NOMEMORY COMPONENTS=" + getLabel() + " COORDINATES=" + mygrid->getComponentName( gdirs[0] );
  for(unsigned i=1;i<gdirs.size();++i) vstring += "," + mygrid->getComponentName( gdirs[i] );
  vstring += " PBC=";
  if( mygrid->isPeriodic(gdirs[0]) ) vstring+="T";
  else vstring+="F";
  for(unsigned i=1;i<gdirs.size();++i){
      if( mygrid->isPeriodic(gdirs[i]) ) vstring+=",T"; else vstring+=",F";
  }
  
  // Create a grid on which to store the output grid
  vesselbase::VesselOptions da("mygrid","",-1,vstring,this);
  Keywords keys; GridFunction::registerKeywords( keys );
  vesselbase::VesselOptions dar( da, keys );
  outgrid = new GridFunction(dar); addVessel( outgrid );
  outgrid->setNoDerivatives();
}

void FindContourSurface::findContourSurface(){
  // Set the boundaries of the output grid
  std::vector<double> fspacing; std::vector<unsigned> snbins( mygrid->getDimension()-1 );
  std::vector<std::string> smin( mygrid->getDimension()-1 ), smax( mygrid->getDimension()-1 );
  for(unsigned i=0;i<gdirs.size();++i){
     smin[i]=mygrid->getMin()[gdirs[i]]; smax[i]=mygrid->getMax()[gdirs[i]];
     snbins[i]=mygrid->getNbin()[gdirs[i]];
  }
  outgrid->setBounds( smin, smax, snbins, fspacing); resizeFunctions();
  // Create a task list if first time
  if( firsttime ){
      std::vector<unsigned> find( mygrid->getDimension() );
      std::vector<unsigned> ind( outgrid->getDimension() );
      for(unsigned i=0;i<outgrid->getNumberOfPoints();++i){
          find.assign( find.size(), 0 ); outgrid->getIndices( i, ind );
          for(unsigned j=0;j<gdirs.size();++j) find[gdirs[j]]=ind[j];
          // Current will be set equal to the start point for this grid index
          addTaskToList( mygrid->getIndex(find) );
      }
      // And prepare the task list
      deactivateAllTasks();
      for(unsigned i=0;i<getFullNumberOfTasks();++i) taskFlags[i]=1;
      lockContributors();
      // Set the direction in which to look for the contour
      direction.resize( mygrid->getDimension(), 0 );
      direction[dir_n] = 0.999999999*mygrid->getGridSpacing()[dir_n];
  }
  // This finds the contour but it is parallized
  runAllTasks(); firsttime=false; 

  // And update the list of active grid points
  if( gbuffer>0 ){
      std::vector<double> dx( mygrid->getGridSpacing() );
      std::vector<double> point( mygrid->getDimension() );
      std::vector<double> lpoint( outgrid->getDimension() );
      std::vector<unsigned> neighbours; unsigned num_neighbours;
      std::vector<unsigned> ugrid_indices( mygrid->getDimension() );
      std::vector<bool> active( mygrid->getNumberOfPoints(), false ); 
      std::vector<unsigned> gbuffer_vec( mygrid->getDimension(), gbuffer );
      for(unsigned i=0;i<outgrid->getNumberOfPoints();++i){
          // Retrieve the coordinates of this grid point
          outgrid->getGridPointCoordinates( i, lpoint );
          point[dir_n] = outgrid->getGridElement( i, 0 );
          // 0.5*dx added here to prevent problems with flooring of grid points
          for(unsigned j=0;j<gdirs.size();++j) point[gdirs[j]]=lpoint[j] + 0.5*dx[gdirs[j]];
          mygrid->getIndices( point, ugrid_indices );
          // Now activate buffer region
          mygrid->getNeighbors( ugrid_indices , gbuffer_vec, num_neighbours, neighbours );
          for(unsigned n=0;n<num_neighbours;++n) active[ neighbours[n] ]=true; 
      }
      mygrid->activateThesePoints( active );
  }
}

void FindContourSurface::performTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const {
  std::vector<unsigned> neighbours; unsigned num_neighbours; unsigned nfound=0; double minv=0, minp;
  std::vector<unsigned> bins_n( mygrid->getNbin() ); unsigned shiftn=current;
  std::vector<unsigned> ind( mygrid->getDimension() ); std::vector<double> point( mygrid->getDimension() );
#ifndef DNDEBUG 
  std::vector<unsigned> oind( outgrid->getDimension() ); outgrid->getIndices( task_index, oind ); 
#endif
  for(unsigned i=0;i<bins_n[dir_n];++i){
#ifndef DNDEBUG
  std::vector<unsigned> base_ind( mygrid->getDimension() ); mygrid->getIndices( shiftn, base_ind );
  for(unsigned j=0;j<gdirs.size();++j) plumed_dbg_assert( base_ind[gdirs[j]]==oind[j] );
#endif
     // Ensure inactive grid points are ignored
     if( mygrid->inactive( shiftn ) ){ shiftn += mygrid->getStride()[dir_n]; continue; }
     // Get the index of the current grid point
     mygrid->getIndices( shiftn, ind );
     // Exit if we are at the edge of the grid
     if( !mygrid->isPeriodic(dir_n) && (ind[dir_n]+1)==bins_n[dir_n] ){
        shiftn += mygrid->getStride()[dir_n]; continue;
     }     

     // Ensure points with inactive neighbours are ignored
     mygrid->getNeighbors( ind, ones, num_neighbours, neighbours );
     bool cycle=false;
     for(unsigned j=0;j<num_neighbours;++j){
         if( mygrid->inactive( neighbours[j]) ){ cycle=true; break; }
     }
     if( cycle ){ shiftn += mygrid->getStride()[dir_n]; continue; }

     // Now get the function value at two points
     double val1=getFunctionValue( shiftn ) - contour; double val2;
     if( (ind[dir_n]+1)==bins_n[dir_n] ) val2 = getFunctionValue( current ) - contour;
     else val2=getFunctionValue( shiftn + mygrid->getStride()[dir_n] ) - contour;

     // Check if the minimum is bracketed 
     if( val1*val2<0 ){ 
         mygrid->getGridPointCoordinates( shiftn, point ); findContour( direction, point ); 
         for(unsigned j=0;j<gdirs.size();++j) myvals.setValue( 1+gdirs[j], point[gdirs[j]] );
         if( firsttime ){
             if( nfound>1 ) error("in first frame found more than one location of dividing surface");
             minp=point[dir_n];
         } else if( nfound==0 ){
             minp=point[dir_n]; minv=std::fabs( point[dir_n] - outgrid->getGridElement( i, 0 ) ); 
         } else {
             double tmp = std::fabs( point[dir_n] - outgrid->getGridElement( i, 0 ) );
             if( tmp<minv ){ minv=tmp; minp=point[dir_n]; }
         }
         nfound++;
     }

     // This moves us on to the next point
     shiftn += mygrid->getStride()[dir_n];
  }
  if( nfound==0 ){
     std::string num; Tools::convert( getStep(), num );
     error("On step " + num + " failed to find required grid point"); 
  }
  outgrid->setGridElement( task_index, 0, minp ); myvals.setValue( 1+dir_n, minp ); myvals.setValue(0, 1.0 );
}

}
}
