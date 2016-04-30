/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014,2015 The plumed team
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
#include "ActionWithGrid.h"

namespace PLMD {
namespace gridtools {

void ActionWithGrid::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  ActionWithVessel::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which the data should be collected and added to the grid");
  keys.add("compulsory","CLEAR","0","the frequency with which to clear all the accumulated data.  The default value "
                                    "of 0 implies that all the data will be used and that the grid will never be cleared");
  keys.add("compulsory","BANDWIDTH","the bandwidths for kernel density esimtation");
  keys.add("compulsory","KERNEL","gaussian","the kernel function you are using.  More details on  the kernels available "
                                            "in plumed plumed can be found in \\ref kernelfunctions.");
  keys.addFlag("UNORMALIZED",false,"output the unormalized density on the grid.");
}

ActionWithGrid::ActionWithGrid( const ActionOptions& ao):
Action(ao),
ActionPilot(ao),
ActionAtomistic(ao),
ActionWithArguments(ao),
ActionWithVessel(ao),
requiresNorm(false),
mygrid(NULL)
{
  if( keywords.exists("CLEAR") ){
      parse("CLEAR",clearstride);
      if( clearstride>0 ){
          if( clearstride%getStride()!=0 ) error("CLEAR parameter must be a multiple of STRIDE");
          log.printf("  clearing grid every %d steps \n",clearstride);
      }
  }
}

void ActionWithGrid::createGrid( const std::string& type, const std::string& inputstr ){
  // Start creating the input for the grid
  std::string vstring = inputstr + " "; 
  if( keywords.exists("KERNEL") ){
      vstring += getKeyword("KERNEL") + " " + getKeyword("BANDWIDTH");
  }
  bool sumflag; parseFlag("UNORMALIZED",sumflag);
  if( sumflag ) vstring += " UNORMALIZED";
  if( clearstride>0 ) vstring += " NOMEMORY";

  vesselbase::VesselOptions da("mygrid","",-1,vstring,this);
  Keywords keys; gridtools::AverageOnGrid::registerKeywords( keys );
  vesselbase::VesselOptions dar( da, keys );
  if( type=="histogram" ){
     mygrid = new HistogramOnGrid(dar); requiresNorm=true; mygrid->setNorm(0);
  } else if( type=="average" ){
     mygrid = new AverageOnGrid(dar); requiresNorm=false;
  } else if( type=="grid" ){
     mygrid = new GridVessel(dar); requiresNorm=false;
  } else {
     plumed_merror("no way to create grid of type " + type );
  } 
}

void ActionWithGrid::finishGridSetup(){
  addVessel( mygrid );
}

void ActionWithGrid::lockRequests(){
  ActionAtomistic::lockRequests();
  ActionWithArguments::lockRequests();
}

void ActionWithGrid::unlockRequests(){
  ActionAtomistic::unlockRequests();
  ActionWithArguments::unlockRequests();
}

void ActionWithGrid::calculateNumericalDerivatives( PLMD::ActionWithValue* vv ){
  error("not possible to compute numerical derivatives for this action");
}

void ActionWithGrid::clearGrid(){
  mygrid->clear();
}

void ActionWithGrid::update(){
  if( getStep()==0 || !onStep() ) return;
  // Clear if it is time to reset
  if( mygrid ){
      if( mygrid->wasreset() ) clearGrid(); 
  }
  // Run all the tasks (if required
  if( prepareForTasks() ) runAllTasks(); 
  // This runs the jobs if we do not have a task list
  else performGridOperations( true );
  // Update the norm if this is necessary
  if( mygrid && requiresNorm ) mygrid->setNorm( 1 + mygrid->getNorm() );
  // Finish the tasks
  finishTaskSet();
  // By resetting here we are ensuring that the grid will be cleared at the start of the next step
  if( mygrid ){
      if( getStride()==0 || (clearstride>0 && getStep()%clearstride==0) ) mygrid->reset();
  }
}

void ActionWithGrid::performTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const {
  // Set the weight of this point
  myvals.setValue( 0, 1.0 ); compute( current, myvals );
}

void ActionWithGrid::performGridOperations( const bool& from_update ){
  plumed_error();
}

}
}
