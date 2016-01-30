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
  keys.add("optional","STRIDE","the frequency with which to output the grid");
  keys.addFlag("USE_ALL_DATA",false,"use the data from the entire trajectory to perform the analysis");
  keys.addFlag("UNNORMALIZED",false,"output an unormalised grid"); 
}

ActionWithInputGrid::ActionWithInputGrid(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionWithVessel(ao),
norm(0.0),
mygrid(NULL)
{
  std::string mlab; parse("GRID",mlab);
  vesselbase::ActionWithVessel* mves= plumed.getActionSet().selectWithLabel<vesselbase::ActionWithVessel*>(mlab);
  if(!mves) error("action labelled " +  mlab + " does not exist or does not have vessels");
  addDependency(mves);

  log.printf("  using grid calculated by action %s \n",mves->getLabel().c_str() );
  for(unsigned i=0;i<mves->getNumberOfVessels();++i){
      mygrid=dynamic_cast<GridVessel*>( mves->getPntrToVessel(i) );
      if( mygrid ) break; 
  }
  if( !mygrid ) error("input action does not calculate a grid");

  parseFlag("UNNORMALIZED",unormalised);
  if( unormalised ){
     log.printf("  working with unormalised grid \n");
     mygrid->switchOffNormalisation(); 
  }

  if( keywords.exists("USE_ALL_DATA") ){
     parseFlag("USE_ALL_DATA",single_run);
     if( !single_run ){
        mves->setAnalysisStride( false, getStride() ); 
        log.printf("  outputting grid every %u steps \n", getStride() );
     } else log.printf("  outputting grid at end of calculation\n");
  }
}

void ActionWithInputGrid::setAnalysisStride( const bool& use_all, const unsigned& astride ){
  single_run=use_all;
  if( !single_run ) setStride( astride );
}

void ActionWithInputGrid::update(){
  if( unormalised ) norm = 1.0;
  else norm=1.0/mygrid->getNorm();

  if( checkAllActive() ){
     for(unsigned i=0;i<mygrid->getNumberOfPoints();++i){
         if( mygrid->inactive(i) ) error("if FIND_CONTOUR is used with BUFFER option then other actions cannot be performed with grid");
     }
  }
  performOperationsWithGrid( true );
}

void ActionWithInputGrid::runFinalJobs(){
  if( !single_run ) return ;
  if( unormalised ) norm = 1.0;
  else norm=1.0/mygrid->getNorm();
  performOperationsWithGrid( false );
}

}
}

