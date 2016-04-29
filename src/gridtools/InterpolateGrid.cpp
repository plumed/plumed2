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
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "ActionWithInputGrid.h"
#include "GridFunction.h"

//+PLUMEDOC GRIDANALYSIS INTERPOLATE_GRID
/*
Interpolate a smooth function stored on a grid onto a grid with a smaller grid spacing.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class InterpolateGrid : public ActionWithInputGrid {
private:
  std::vector<unsigned> nbin;
  std::vector<double> gspacing;
  GridFunction* outgrid;
public:
  static void registerKeywords( Keywords& keys );
  explicit InterpolateGrid(const ActionOptions&ao);
  void performOperationsWithGrid( const bool& from_update );
  unsigned getNumberOfDerivatives(){ return 0; }
  unsigned getNumberOfQuantities() const ;
  void performTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const ;
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(InterpolateGrid,"INTERPOLATE_GRID")

void InterpolateGrid::registerKeywords( Keywords& keys ){
  ActionWithInputGrid::registerKeywords( keys );
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
}

InterpolateGrid::InterpolateGrid(const ActionOptions&ao):
Action(ao),
ActionWithInputGrid(ao)
{
  plumed_assert( mygrid->getNumberOfComponents()==1 );
  if( mygrid->noDerivatives() ) error("cannot interpolate a grid that does not have derivatives"); 
  // Create the input from the old string
  std::string vstring = "COMPONENTS=" + getLabel() + " " + mygrid->getInputString();

  // Create a grid
  vesselbase::VesselOptions da("mygrid","",-1,vstring,this);
  Keywords keys; GridFunction::registerKeywords( keys );
  vesselbase::VesselOptions dar( da, keys );
  outgrid = new GridFunction(dar); addVessel( outgrid ); 

  parseVector("GRID_BIN",nbin); parseVector("GRID_SPACING",gspacing);
  if( nbin.size()!=mygrid->getDimension() && gspacing.size()!=mygrid->getDimension() ){
      error("GRID_BIN or GRID_SPACING must be set");
  } 
  outgrid->setBounds( mygrid->getMin(), mygrid->getMax(), nbin, gspacing ); 
  resizeFunctions();

  // Now create task list
  for(unsigned i=0;i<outgrid->getNumberOfPoints();++i) addTaskToList(i);
  // And activate all tasks
  deactivateAllTasks(); 
  for(unsigned i=0;i<outgrid->getNumberOfPoints();++i) taskFlags[i]=1;
  lockContributors();
}

unsigned InterpolateGrid::getNumberOfQuantities() const {
  return 2 + mygrid->getDimension();
}

void InterpolateGrid::performOperationsWithGrid( const bool& from_update ){
  outgrid->setBounds( mygrid->getMin(), mygrid->getMax(), nbin, gspacing );
  outgrid->clear(); outgrid->setNorm( mygrid->getNorm() ); runAllTasks(); outgrid->reset();
}

void InterpolateGrid::performTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const {
  std::vector<double> pos( mygrid->getDimension() ); outgrid->getGridPointCoordinates( current, pos );
  std::vector<double> der( mygrid->getDimension() ); double val = getFunctionValueAndDerivatives( pos, der );
  myvals.setValue( 0, 1.0 ); myvals.setValue(1, val ); 
  for(unsigned i=0;i<mygrid->getDimension();++i) myvals.setValue( 2+i, der[i] ); 
}

}
}
