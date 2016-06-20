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
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "ActionWithInputGrid.h"

//+PLUMEDOC GRIDANALYSIS INTERPOLATE_GRID
/*
Interpolate a smooth function stored on a grid onto a grid with a smaller grid spacing.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class InterpolateGrid : public ActionWithInputGrid {
public:
  static void registerKeywords( Keywords& keys );
  explicit InterpolateGrid(const ActionOptions&ao);
  unsigned getNumberOfQuantities() const ;
  void compute( const unsigned& current, MultiValue& myvals ) const ;
  bool isPeriodic(){ return false; }
};

PLUMED_REGISTER_ACTION(InterpolateGrid,"INTERPOLATE_GRID")

void InterpolateGrid::registerKeywords( Keywords& keys ){
  ActionWithInputGrid::registerKeywords( keys );
  keys.add("optional","GRID_BIN","the number of bins for the grid");
  keys.add("optional","GRID_SPACING","the approximate grid spacing (to be used as an alternative or together with GRID_BIN)");
  keys.remove("KERNEL"); keys.remove("BANDWIDTH"); 
}

InterpolateGrid::InterpolateGrid(const ActionOptions&ao):
Action(ao),
ActionWithInputGrid(ao)
{
  plumed_assert( ingrid->getNumberOfComponents()==1 );
  if( ingrid->noDerivatives() ) error("cannot interpolate a grid that does not have derivatives"); 
  // Create the input from the old string
  createGrid( "grid", "COMPONENTS=" + getLabel() + " " + ingrid->getInputString()  );

  std::vector<unsigned> nbin; parseVector("GRID_BIN",nbin); 
  std::vector<double> gspacing; parseVector("GRID_SPACING",gspacing);
  if( nbin.size()!=ingrid->getDimension() && gspacing.size()!=ingrid->getDimension() ){
      error("GRID_BIN or GRID_SPACING must be set");
  } 

  // Need this for creation of tasks
  mygrid->setBounds( ingrid->getMin(), ingrid->getMax(), nbin, gspacing ); 
  setAveragingAction( mygrid, true );

  // Now create task list
  for(unsigned i=0;i<mygrid->getNumberOfPoints();++i) addTaskToList(i);
  // And activate all tasks
  deactivateAllTasks(); 
  for(unsigned i=0;i<mygrid->getNumberOfPoints();++i) taskFlags[i]=1;
  lockContributors();
}

unsigned InterpolateGrid::getNumberOfQuantities() const {
  return 2 + ingrid->getDimension();
}

void InterpolateGrid::compute( const unsigned& current, MultiValue& myvals ) const {
  std::vector<double> pos( mygrid->getDimension() ); mygrid->getGridPointCoordinates( current, pos );
  std::vector<double> der( mygrid->getDimension() ); double val = getFunctionValueAndDerivatives( pos, der );
  myvals.setValue( 0, 1.0 ); myvals.setValue(1, val ); 
  for(unsigned i=0;i<mygrid->getDimension();++i) myvals.setValue( 2+i, der[i] ); 
}

}
}
