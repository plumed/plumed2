/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015,2016 The plumed team
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

//+PLUMEDOC GRIDANALYSIS FIND_CONTOUR
/*
Find an isocontour in a smooth function.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class FindContour : public ContourFindingBase {
private:
  bool firsttime;
  unsigned gbuffer;
public:
  static void registerKeywords( Keywords& keys );
  explicit FindContour(const ActionOptions&ao);
  bool checkAllActive() const { return gbuffer==0; }
  void prepareForAveraging();
  bool isPeriodic(){ return false; }
  void compute( const unsigned& current, MultiValue& myvals ) const ;
  void finishAveraging();
};

PLUMED_REGISTER_ACTION(FindContour,"FIND_CONTOUR")

void FindContour::registerKeywords( Keywords& keys ){
  ContourFindingBase::registerKeywords( keys );
// We want a better way of doing this bit
  keys.add("compulsory","BUFFER","0","number of buffer grid points around location where grid was found on last step.  If this is zero the full grid is calculated on each step");
}

FindContour::FindContour(const ActionOptions&ao):
Action(ao),
ContourFindingBase(ao),
firsttime(true)
{

  parse("BUFFER",gbuffer);
  if( gbuffer>0 ) log.printf("  after first step a subset of only %u grid points around where the countour was found will be checked\n",gbuffer);
  checkRead();
}

void FindContour::prepareForAveraging(){
  // Create a task list if first time
  if( firsttime ){
      for(unsigned i=0;i<ingrid->getDimension()*ingrid->getNumberOfPoints();++i) addTaskToList( i );
  }
  firsttime=false; deactivateAllTasks();

  // We now need to identify the grid points that we need to search through
  std::vector<unsigned> nbin( ingrid->getNbin() );
  std::vector<unsigned> ind( ingrid->getDimension() );
  std::vector<unsigned> ones( ingrid->getDimension(), 1 );
  unsigned num_neighbours; std::vector<unsigned> neighbours;
  for(unsigned i=0;i<ingrid->getNumberOfPoints();++i){
     // Ensure inactive grid points are ignored
     if( ingrid->inactive(i) ) continue;

     // Get the index of the current grid point
     ingrid->getIndices( i, ind );
     ingrid->getNeighbors( ind, ones, num_neighbours, neighbours );
     bool cycle=false;
     for(unsigned j=0;j<num_neighbours;++j){
         if( ingrid->inactive( neighbours[j]) ){ cycle=true; break; }
     }
     if( cycle ) continue;

     // Get the value of a point on the grid
     double val1=getFunctionValue( i ) - contour;
     bool edge=false;
     for(unsigned j=0;j<ingrid->getDimension();++j){
         // Make sure we don't search at the edge of the grid
         if( !ingrid->isPeriodic(j) && (ind[j]+1)==nbin[j] ) continue;
         else if( (ind[j]+1)==nbin[j] ){ edge=true; ind[j]=0; }
         else ind[j]+=1;
         double val2=getFunctionValue( ind ) - contour;
         if( val1*val2<0 ) taskFlags[ ingrid->getDimension()*i + j ] = 1;
         if( ingrid->isPeriodic(j) && edge ){ edge=false; ind[j]=nbin[j]-1; }
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
  std::vector<double> direction( ingrid->getDimension() , 0 );
  direction[gdir] = 0.999999999*ingrid->getGridSpacing()[gdir];

  // Now find the contour
  findContour( direction, point );
  // And transfer to the store data vessel
  for(unsigned i=0;i<ingrid->getDimension();++i) myvals.setValue( 1+i, point[i] );
}

void FindContour::finishAveraging(){
  ContourFindingBase::finishAveraging();
  // And update the list of active grid points
  if( gbuffer>0 ){
      std::vector<unsigned> neighbours; unsigned num_neighbours;
      std::vector<unsigned> ugrid_indices( ingrid->getDimension() );
      std::vector<bool> active( ingrid->getNumberOfPoints(), false );  
      std::vector<unsigned> gbuffer_vec( ingrid->getDimension(), gbuffer );
      for(unsigned i=0;i<getCurrentNumberOfActiveTasks();++i){
          // Get the point we are operating on
          unsigned ipoint = std::floor( getActiveTask(i) / ingrid->getDimension() );
          // Get the indices of this point
          ingrid->getIndices( ipoint, ugrid_indices );
          // Now activate buffer region
          ingrid->getNeighbors( ugrid_indices , gbuffer_vec, num_neighbours, neighbours );
          for(unsigned n=0;n<num_neighbours;++n) active[ neighbours[n] ]=true;  
      }
      ingrid->activateThesePoints( active );
  }
}

}
}
