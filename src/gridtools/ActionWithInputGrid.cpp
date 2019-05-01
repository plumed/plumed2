/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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
#include "ActionWithInputGrid.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace gridtools {

void ActionWithInputGrid::registerKeywords( Keywords& keys ) {
  ActionWithGrid::registerKeywords( keys );
  keys.add("compulsory","GRID","the action that creates the input grid you would like to use");
  keys.add("optional","COMPONENT","if your input is a vector field use this to specify the component of the input vector field for which you wish to use");
}

ActionWithInputGrid::ActionWithInputGrid(const ActionOptions&ao):
  Action(ao),
  ActionWithGrid(ao),
  ingrid(NULL)
{
  std::string mlab; parse("GRID",mlab);
  vesselbase::ActionWithVessel* mves= plumed.getActionSet().selectWithLabel<vesselbase::ActionWithVessel*>(mlab);
  if(!mves) error("action labelled " +  mlab + " does not exist or does not have vessels");
  addDependency(mves);

  for(unsigned i=0; i<mves->getNumberOfVessels(); ++i) {
    ingrid=dynamic_cast<GridVessel*>( mves->getPntrToVessel(i) );
    if( ingrid ) break;
  }
  if( !ingrid ) error("input action does not calculate a grid");

  if( ingrid->getNumberOfComponents()==1 ) {
    mycomp=0;
  } else {
    int tcomp=-1; parse("COMPONENT",tcomp);
    if( tcomp<0 ) error("component of vector field was not specified - use COMPONENT keyword");
    mycomp=tcomp;
  }
  log.printf("  using %uth component of grid calculated by action %s \n",mycomp,mves->getLabel().c_str() );
}

void ActionWithInputGrid::clearAverage() {
  if( mygrid->getType()=="flat" ) mygrid->setBounds( ingrid->getMin(), ingrid->getMax(), mygrid->getNbin(), mygrid->getGridSpacing() );
  ActionWithAveraging::clearAverage();
}

void ActionWithInputGrid::prepareForAveraging() {
  if( checkAllActive() ) {
    for(unsigned i=0; i<ingrid->getNumberOfPoints(); ++i) {
      if( ingrid->inactive(i) ) error("if FIND_CONTOUR is used with BUFFER option then other actions cannot be performed with grid");
    }
  }
}

void ActionWithInputGrid::performOperations( const bool& from_update ) {
  prepareForAveraging(); runAllTasks();
}

}
}

