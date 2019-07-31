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
#include "core/PlumedMain.h"
#include "ActionWithInputGrid.h"

//+PLUMEDOC GRIDANALYSIS INTERPOLATE_GRID
/*
Interpolate a smooth function stored on a grid onto a grid with a smaller grid spacing.

This action takes a function evaluated on a grid as input and can be used to interpolate the values of that
function on to a finer grained grid.  The interpolation within this algorithm is done using splines.

\par Examples

The input below can be used to post process a trajectory.  It calculates a \ref HISTOGRAM as a function the
distance between atoms 1 and 2 using kernel density estimation.  During the calculation the values of the kernels
are evaluated at 100 points on a uniform grid between 0.0 and 3.0.  Prior to outputting this function at the end of the
simulation this function is interpolated onto a finer grid of 200 points between 0.0 and 3.0.

\plumedfile
x: DISTANCE ATOMS=1,2
hA1: HISTOGRAM ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1
ii: INTERPOLATE_GRID GRID=hA1 GRID_BIN=200
DUMPGRID GRID=ii FILE=histo.dat
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class InterpolateGrid : public ActionWithInputGrid {
public:
  static void registerKeywords( Keywords& keys );
  explicit InterpolateGrid(const ActionOptions&ao);
  unsigned getNumberOfQuantities() const override;
  void compute( const unsigned& current, MultiValue& myvals ) const override;
  bool isPeriodic() override { return false; }
};

PLUMED_REGISTER_ACTION(InterpolateGrid,"INTERPOLATE_GRID")

void InterpolateGrid::registerKeywords( Keywords& keys ) {
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
  auto grid=createGrid( "grid", "COMPONENTS=" + getLabel() + " " + ingrid->getInputString()  );
  // notice that createGrid also sets mygrid=grid.get()

  std::vector<unsigned> nbin; parseVector("GRID_BIN",nbin);
  std::vector<double> gspacing; parseVector("GRID_SPACING",gspacing);
  if( nbin.size()!=ingrid->getDimension() && gspacing.size()!=ingrid->getDimension() ) {
    error("GRID_BIN or GRID_SPACING must be set");
  }

  // Need this for creation of tasks
  grid->setBounds( ingrid->getMin(), ingrid->getMax(), nbin, gspacing );
  setAveragingAction( std::move(grid), true );

  // Now create task list
  for(unsigned i=0; i<mygrid->getNumberOfPoints(); ++i) addTaskToList(i);
  // And activate all tasks
  deactivateAllTasks();
  for(unsigned i=0; i<mygrid->getNumberOfPoints(); ++i) taskFlags[i]=1;
  lockContributors();
}

unsigned InterpolateGrid::getNumberOfQuantities() const {
  return 2 + ingrid->getDimension();
}

void InterpolateGrid::compute( const unsigned& current, MultiValue& myvals ) const {
  std::vector<double> pos( mygrid->getDimension() ); mygrid->getGridPointCoordinates( current, pos );
  std::vector<double> der( mygrid->getDimension() ); double val = getFunctionValueAndDerivatives( pos, der );
  myvals.setValue( 0, 1.0 ); myvals.setValue(1, val );
  for(unsigned i=0; i<mygrid->getDimension(); ++i) myvals.setValue( 2+i, der[i] );
}

}
}
