/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016,2017 The plumed team
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
#include "ActionWithAveraging.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

namespace PLMD {
namespace vesselbase {

void ActionWithAveraging::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys ); ActionPilot::registerKeywords( keys ); ActionAtomistic::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys ); ActionWithValue::registerKeywords( keys ); ActionWithVessel::registerKeywords( keys );
  keys.add("compulsory","STRIDE","1","the frequency with which the data should be collected and added to the quantity being averaged");
  keys.add("compulsory","CLEAR","0","the frequency with which to clear all the accumulated data.  The default value "
           "of 0 implies that all the data will be used and that the grid will never be cleared");
  keys.add("optional","LOGWEIGHTS","list of actions that calculates log weights that should be used to weight configurations when calculating averages");
  keys.add("compulsory","NORMALIZATION","true","This controls how the data is normalized it can be set equal to true, false or ndata.  The differences between these options are explained in the manual page for \\ref HISTOGRAM");
  keys.remove("NUMERICAL_DERIVATIVES");
}

ActionWithAveraging::ActionWithAveraging( const ActionOptions& ao ):
  Action(ao),
  ActionPilot(ao),
  ActionAtomistic(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao),
  ActionWithVessel(ao),
  myaverage(NULL),
  useRunAllTasks(false),
  clearstride(0),
  lweight(0),cweight(0)
{
  if( keywords.exists("CLEAR") ) {
    parse("CLEAR",clearstride);
    if( clearstride>0 ) {
      if( clearstride%getStride()!=0 ) error("CLEAR parameter must be a multiple of STRIDE");
      log.printf("  clearing grid every %u steps \n",clearstride);
    }
  }
  if( keywords.exists("LOGWEIGHTS") ) {
    std::vector<std::string> wwstr; parseVector("LOGWEIGHTS",wwstr);
    if( wwstr.size()>0 ) log.printf("  reweighting using weights from ");
    std::vector<Value*> arg( getArguments() );
    for(unsigned i=0; i<wwstr.size(); ++i) {
      ActionWithValue* val = plumed.getActionSet().selectWithLabel<ActionWithValue*>(wwstr[i]);
      if( !val ) error("could not find value named");
      weights.push_back( val->copyOutput(val->getLabel()) );
      arg.push_back( val->copyOutput(val->getLabel()) );
      log.printf("%s ",wwstr[i].c_str() );
    }
    if( wwstr.size()>0 ) log.printf("\n");
    else log.printf("  weights are all equal to one\n");
    requestArguments( arg );
  }
  if( keywords.exists("NORMALIZATION") ) {
    std::string normstr; parse("NORMALIZATION",normstr);
    if( normstr=="true" ) normalization=t;
    else if( normstr=="false" ) normalization=f;
    else if( normstr=="ndata" ) normalization=ndata;
    else error("invalid instruction for NORMALIZATION flag should be true, false, or ndata");
  }
}

void ActionWithAveraging::setAveragingAction( AveragingVessel* av_vessel, const bool& usetasks ) {
  myaverage=av_vessel; addVessel( myaverage );
  useRunAllTasks=usetasks; resizeFunctions();
}

void ActionWithAveraging::lockRequests() {
  ActionAtomistic::lockRequests();
  ActionWithArguments::lockRequests();
}

void ActionWithAveraging::unlockRequests() {
  ActionAtomistic::unlockRequests();
  ActionWithArguments::unlockRequests();
}

void ActionWithAveraging::calculateNumericalDerivatives(PLMD::ActionWithValue*) {
  error("not possible to compute numerical derivatives for this action");
}

void ActionWithAveraging::update() {
  if( (clearstride!=1 && getStep()==0) || !onStep() ) return;
  // Clear if it is time to reset
  if( myaverage ) {
    if( myaverage->wasreset() ) clearAverage();
  }
  // Calculate the weight for all reweighting
  if ( weights.size()>0 ) {
    double sum=0; for(unsigned i=0; i<weights.size(); ++i) sum+=weights[i]->get();
    lweight=sum; cweight = exp( sum );
  } else {
    lweight=0; cweight=1.0;
  }
  // Prepare to do the averaging
  prepareForAveraging();
  // Run all the tasks (if required
  if( useRunAllTasks ) runAllTasks();
  // This the averaging if it is not done using task list
  else performOperations( true );
  // Update the norm
  double normt = cweight; if( normalization==ndata ) normt = 1;
  if( myaverage ) myaverage->setNorm( normt + myaverage->getNorm() );
  // Finish the averaging
  finishAveraging();
  // By resetting here we are ensuring that the grid will be cleared at the start of the next step
  if( myaverage ) {
    if( getStride()==0 || (clearstride>0 && getStep()%clearstride==0) ) myaverage->reset();
  }
}

void ActionWithAveraging::clearAverage() { plumed_assert( myaverage->wasreset() ); myaverage->clear(); }

void ActionWithAveraging::performOperations( const bool& from_update ) { plumed_error(); }


}
}
