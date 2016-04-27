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
#include "ActionWithInputGrid.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace gridtools {

void ActionWithInputGrid::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  vesselbase::ActionWithVessel::registerKeywords( keys );
  keys.add("compulsory","GRID","the action that creates the input grid you would like to use");
  keys.add("compulsory","STRIDE","0","the frequency with which to output the grid.  Default of zero means analyse all data");
  keys.add("optional","COMPONENT","if your input is a vector field use this to specifiy the component of the input vector field for which you wish to use");
}

ActionWithInputGrid::ActionWithInputGrid(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionWithVessel(ao),
mygrid(NULL)
{
  std::string mlab; parse("GRID",mlab);
  mves= plumed.getActionSet().selectWithLabel<vesselbase::ActionWithVessel*>(mlab);
  if(!mves) error("action labelled " +  mlab + " does not exist or does not have vessels");
  addDependency(mves);

  ActionPilot* ap=dynamic_cast<ActionPilot*>( mves );
  if( ap ){
     // This is perhaps necessary check what happens if Print is called between kernel adding GAT  
     if( getStride()>0 && getStride()%ap->getStride()!=0 ) error("mismatch between strides in " + ap->getLabel() + " and " +  getLabel() );
  }
  for(unsigned i=0;i<mves->getNumberOfVessels();++i){
      mygrid=dynamic_cast<GridVessel*>( mves->getPntrToVessel(i) );
      if( mygrid ) break; 
  }
  if( !mygrid ) error("input action does not calculate a grid");

  if( mygrid->getNumberOfComponents()==1 ){
     mycomp=0;
  } else {
     int tcomp=-1; parse("COMPONENT",tcomp);
     if( tcomp<0 ) error("component of vector field was not specified - use COMPONENT keyword");
     mycomp=tcomp;
  }
  log.printf("  using %uth component of grid calculated by action %s \n",mycomp,mves->getLabel().c_str() );

  if( getStride()>0 ) log.printf("  outputting every %u steps \n", getStride() );
}

void ActionWithInputGrid::update(){
  // Don't analyse the first frame in the trajectory
  if( getStep()==0 ) return;

  if( checkAllActive() ){
     for(unsigned i=0;i<mygrid->getNumberOfPoints();++i){
         if( mygrid->inactive(i) ) error("if FIND_CONTOUR is used with BUFFER option then other actions cannot be performed with grid");
     }
  }
  performOperationsWithGrid( true );
}

void ActionWithInputGrid::runFinalJobs(){
  if( getStride()>0 ) return;
  performOperationsWithGrid( false );
}

void ActionWithInputGrid::invertTask( const std::vector<double>& indata, std::vector<double>& outdata ){
  plumed_merror("Shouldn't appear here");
}

}
}

