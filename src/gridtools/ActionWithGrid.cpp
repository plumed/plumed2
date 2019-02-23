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
#include "ActionWithGrid.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace gridtools {

void ActionWithGrid::registerKeywords( Keywords& keys ) {
  vesselbase::ActionWithAveraging::registerKeywords( keys );
  keys.add("compulsory","BANDWIDTH","the bandwidths for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian","the kernel function you are using.  More details on  the kernels available "
           "in plumed plumed can be found in \\ref kernelfunctions.");
  keys.add("optional","CONCENTRATION","the concentration parameter for Von Mises-Fisher distributions");
}

ActionWithGrid::ActionWithGrid( const ActionOptions& ao):
  Action(ao),
  ActionWithAveraging(ao),
  mygrid(NULL)
{
}

std::unique_ptr<GridVessel> ActionWithGrid::createGrid( const std::string& type, const std::string& inputstr ) {
  // Start creating the input for the grid
  std::string vstring = inputstr;
  if( keywords.exists("KERNEL") ) {
    std::string vconc; parse("CONCENTRATION",vconc);
    if( vconc.length()>0 ) {
      vstring += " TYPE=fibonacci CONCENTRATION=" + vconc;
    } else {
      std::string kstring; parse("KERNEL",kstring);
      if( kstring=="DISCRETE" ) vstring += " KERNEL=" + kstring;
      else vstring += " KERNEL=" + kstring + " " + getKeyword("BANDWIDTH");
    }
  }
  vesselbase::VesselOptions da("mygrid","",-1,vstring,this);
  Keywords keys; gridtools::AverageOnGrid::registerKeywords( keys );
  vesselbase::VesselOptions dar( da, keys );
  std::unique_ptr<GridVessel> grid;
  if( type=="histogram" ) {
    grid.reset( new HistogramOnGrid(dar) );
  } else if( type=="average" ) {
    grid.reset( new AverageOnGrid(dar) );
  } else if( type=="grid" ) {
    grid.reset( new GridVessel(dar) );
  } else {
    plumed_merror("no way to create grid of type " + type );
  }
  mygrid=grid.get();
  return grid;
}

void ActionWithGrid::turnOnDerivatives() {
  needsDerivatives(); ActionWithValue::turnOnDerivatives();
  if( getStride()==1 ) setStride(0);
  else if( getStride()!=0 ) error("conflicting instructions for grid - stride was set but must be evaluated on every step for derivatives - remove STRIDE keyword");
  if( clearstride>1 ) error("conflicting instructions for grid - CLEAR was set but grid must be reset on every step for derivatives - remove CLEAR keyword" );
  if( weights.size()>0 ) error("conflicting instructions for grid - LOGWEIGHTS was set but weights are not considered when derivatives of grid are evaluated - remove LOGWEIGHTS keyword");
}

void ActionWithGrid::calculate() {
  // Do nothing if derivatives are not required
  if( doNotCalculateDerivatives() ) return;
  // Clear on every step
  if( mygrid ) clearAverage();
  // Should not be any reweighting so just set these accordingly
  lweight=0; cweight=1.0;
  // Prepare to do the averaging
  prepareForAveraging();
  // Run all the tasks (if required
  if( useRunAllTasks ) runAllTasks();
  // This the averaging if it is not done using task list
  else performOperations( true );
  // Update the norm
  if( mygrid ) mygrid->setNorm( cweight );
  // Finish the averaging
  finishAveraging();
  // And reset for next step
  if( mygrid ) mygrid->reset();
}

void ActionWithGrid::runTask( const unsigned& current, MultiValue& myvals ) const {
  // Set the weight of this point
  myvals.setValue( 0, cweight ); compute( current, myvals );
}

}
}
