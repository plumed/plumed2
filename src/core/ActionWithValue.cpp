/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#include "ActionWithValue.h"
#include "ActionWithArguments.h"
#include "ActionAtomistic.h"
#include "ActionWithVirtualAtom.h"
#include "PlumedMain.h"
#include "ActionSet.h"
#include "ActionRegister.h"
#include "tools/Stopwatch.h"
#include "tools/Exception.h"
#include "tools/OpenMP.h"

namespace PLMD {

void ActionWithValue::registerKeywords(Keywords& keys) {
  keys.setComponentsIntroduction("By default the value of the calculated quantity can be referenced elsewhere in the "
                                 "input file by using the label of the action.  Alternatively this Action can be used "
                                 "to calculate the following quantities by employing the keywords listed "
                                 "below.  These quantities can be referenced elsewhere in the input by using this Action's "
                                 "label followed by a dot and the name of the quantity required from the list below.");
  keys.addFlag("NUMERICAL_DERIVATIVES", false, "calculate the derivatives for these quantities numerically");
  keys.addFlag("SERIAL",false,"do the calculation in serial.  Do not parallelize");
  keys.addFlag("TIMINGS",false,"output information on the timings of the various parts of the calculation");
}

void ActionWithValue::noAnalyticalDerivatives(Keywords& keys) {
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.addFlag("NUMERICAL_DERIVATIVES",true,"analytical derivatives are not implemented for this keyword so numerical derivatives are always used");
}

void ActionWithValue::componentsAreNotOptional(Keywords& keys) {
  keys.setComponentsIntroduction("By default this Action calculates the following quantities. These quantities can "
                                 "be referenced elsewhere in the input by using this Action's label followed by a "
                                 "dot and the name of the quantity required from the list below.");
}

void ActionWithValue::useCustomisableComponents(Keywords& keys) {
  keys.setComponentsIntroduction("The names of the components in this action can be customized by the user in the "
                                 "actions input file.  However, in addition to the components that can be customized the "
                                 "following quantities will always be output");
}

ActionWithValue::ActionWithValue(const ActionOptions&ao):
  Action(ao),
  noderiv(true),
  numericalDerivatives(false),
  no_openmp(plumed.getMDEngine()=="plumed"),
  serial(false),
  timers(false),
  nactive_tasks(0),
  thisAsActionWithArguments(NULL),
  thisAsActionWithVatom(NULL),
  action_to_do_before(NULL),
  action_to_do_after(NULL),
  atom_action_to_do_after(NULL)
{
  if( keywords.exists("NUMERICAL_DERIVATIVES") ) parseFlag("NUMERICAL_DERIVATIVES",numericalDerivatives);
  if(numericalDerivatives) log.printf("  using numerical derivatives\n");
  if( keywords.exists("SERIAL") ) parseFlag("SERIAL",serial);
  else serial=true;
  if( keywords.exists("TIMINGS") ) {
    parseFlag("TIMINGS",timers);
    if( timers ) { stopwatch.start(); stopwatch.pause(); }
  }
}

ActionWithValue::~ActionWithValue() {
  stopwatch.start(); stopwatch.stop();
  if(timers) {
    stopwatch.start(); stopwatch.stop();
    log.printf("timings for action %s with label %s \n", getName().c_str(), getLabel().c_str() );
    log<<stopwatch;
  }
}

ActionWithValue* ActionWithValue::getActionThatCalculates() {
  // Return this if we have no dependencies
  if( !action_to_do_before ) return this;
  // Recursively go through actiosn before
  return action_to_do_before->getActionThatCalculates();
}

void ActionWithValue::getAllActionLabelsInChain( std::vector<std::string>& mylabels ) const {
  bool found = false ;
  for(unsigned i=0; i<mylabels.size(); ++i) {
    if( getLabel()==mylabels[i] ) { found=true; }
  }
  if( !found ) mylabels.push_back( getLabel() );
  if( action_to_do_after ) action_to_do_after->getAllActionLabelsInChain( mylabels );
}

void ActionWithValue::getAllActionsRequired( std::vector<const ActionWithValue*>& allvals ) const {
  const ActionWithArguments* aa = dynamic_cast<const ActionWithArguments*>( this );
  if( !aa ) {
      bool found = false;
      for(unsigned k=0;k<allvals.size();++k) {
          if( allvals[k]==this ) { found=true; break; }
      }
      if( !found ) allvals.push_back( this );
  } else {
     for(unsigned j=0;j<aa->getNumberOfArguments();++j) {
        ((aa->getPntrToArgument(j))->getPntrToAction())->getAllActionsRequired( allvals );
        bool found = false;
        for(unsigned k=0;k<allvals.size();++k) {
            if( allvals[k]==(aa->getPntrToArgument(j))->getPntrToAction() ) { found=true; break; }
        }
        if( !found ) allvals.push_back( (aa->getPntrToArgument(j))->getPntrToAction() );
     }
  }
}

bool ActionWithValue::addActionToChain( const std::vector<std::string>& alabels, ActionWithValue* act ) {
  if( action_to_do_after ) { bool state=action_to_do_after->addActionToChain( alabels, act ); return state; }

  // Check action is not already in chain
  std::vector<std::string> mylabels; getActionThatCalculates()->getAllActionLabelsInChain( mylabels );
  for(unsigned i=0; i<mylabels.size(); ++i) {
    if( act->getLabel()==mylabels[i] ) return true;
  }

  // Check that everything that is required has been calculated
  for(unsigned i=0; i<alabels.size(); ++i) {
    bool found=false;
    for(unsigned j=0; j<mylabels.size(); ++j) {
      if( alabels[i]==mylabels[j] ) { found=true; break; }
    }
    if( !found ) {
      ActionWithValue* av=plumed.getActionSet().selectWithLabel<ActionWithValue*>( alabels[i] );
      bool storingall=true;
      for(unsigned j=0; j<av->getNumberOfComponents(); ++j) {
        if( !(av->getPntrToOutput(j))->storedata ) storingall=false;
      }
      if( !storingall ) return false;
    }
  }
  // This checks that there is nothing that will cause problems in the chain
  mylabels.resize(0); getActionThatCalculates()->getAllActionLabelsInChain( mylabels );
  for(unsigned i=0;i<mylabels.size();++i) {
      ActionWithValue* av1=plumed.getActionSet().selectWithLabel<ActionWithValue*>( mylabels[i] ); 
      for(unsigned j=0;j<i;++j) {
          ActionWithValue* av2=plumed.getActionSet().selectWithLabel<ActionWithValue*>( mylabels[j] );
          if( !av1->canBeAfterInChain( av2 ) ) error("must calculate " + mylabels[j] + " before " + mylabels[i] );
      }
  }
  action_to_do_after=act; 
  atom_action_to_do_after=dynamic_cast<ActionAtomistic*>( action_to_do_after ); 
  act->action_to_do_before=this;
  return true;
}

void ActionWithValue::clearInputForces() {
  for(unsigned i=0; i<values.size(); i++) values[i]->clearInputForce();
}

void ActionWithValue::clearDerivatives( const bool& force ) {
  // This ensures we do not clear derivatives calculated in a chain until the next time the first
  // action in the chain is called
  if( !force && action_to_do_before ) return ;
  unsigned nt = OpenMP::getNumThreads();
  #pragma omp parallel num_threads(nt)
  {
    #pragma omp for
    for(unsigned i=0; i<values.size(); i++) values[i]->clearDerivatives();
  }
  if( action_to_do_after ) action_to_do_after->clearDerivatives( true );
}

// -- These are the routine for copying the value pointers to other classes -- //

bool ActionWithValue::exists( const std::string& name ) const {
  for(unsigned i=0; i<values.size(); ++i) {
    if (values[i]->name==name) return true;
  }
  return false;
}

Value* ActionWithValue::copyOutput( const std::string& name ) const {
  for(unsigned i=0; i<values.size(); ++i) {
    if (values[i]->name==name) return values[i].get();
  }
  plumed_merror("there is no pointer with name " + name);
  return NULL;
}

Value* ActionWithValue::copyOutput( const unsigned& n ) const {
  plumed_massert(n<values.size(),"you have requested a pointer that is out of bounds");
  return values[n].get();
}

// -- HERE WE HAVE THE STUFF FOR THE DEFAULT VALUE -- //

void ActionWithValue::addValue( const std::vector<unsigned>& shape ) {
  plumed_massert(values.empty(),"You have already added the default value for this action");
  values.emplace_back(Tools::make_unique<Value>(this,getLabel(), false, shape ) );
}

void ActionWithValue::addValueWithDerivatives( const std::vector<unsigned>& shape ) {
  plumed_massert(values.empty(),"You have already added the default value for this action");
  values.emplace_back(Tools::make_unique<Value>(this,getLabel(), true, shape ) );
}

void ActionWithValue::setNotPeriodic() {
  plumed_massert(values.size()==1,"The number of components is not equal to one");
  plumed_massert(values[0]->name==getLabel(), "The value you are trying to set is not the default");
  values[0]->min=0; values[0]->max=0;
  values[0]->setupPeriodicity();
}

void ActionWithValue::setPeriodic( const std::string& min, const std::string& max ) {
  plumed_massert(values.size()==1,"The number of components is not equal to one");
  plumed_massert(values[0]->name==getLabel(), "The value you are trying to set is not the default");
  values[0]->setDomain( min, max );
}

Value* ActionWithValue::getPntrToValue() {
  plumed_dbg_massert(values.size()==1,"The number of components is not equal to one");
  plumed_dbg_massert(values[0]->name==getLabel(), "The value you are trying to retrieve is not the default");
  return values[0].get();
}

// -- HERE WE HAVE THE STUFF FOR NAMED VALUES / COMPONENTS -- //

void ActionWithValue::addComponent( const std::string& name, const std::vector<unsigned>& shape ) {
  if( !keywords.outputComponentExists(name,true) ) {
    plumed_merror("a description of component " + name + " has not been added to the manual. Components should be registered like keywords in "
                  "registerKeywords as described in the developer docs.");
  }
  std::string thename; thename=getLabel() + "." + name;
  for(unsigned i=0; i<values.size(); ++i) {
    if( !allowComponentsAndValue() ) plumed_massert(values[i]->name!=getLabel(),"Cannot mix single values with components");
    plumed_massert(values[i]->name!=thename,"there is already a value with this name: "+thename);
    plumed_massert(values[i]->name!=thename&&name!="bias","Since PLUMED 2.3 the component 'bias' is automatically added to all biases by the general constructor!\n"
                   "Remove the line addComponent(\"bias\") from your bias.");
  }
  values.emplace_back(Tools::make_unique<Value>(this,thename, false, shape ) );
  std::string msg="  added component to this action:  "+thename+" \n";
  log.printf(msg.c_str());
}

void ActionWithValue::addComponentWithDerivatives( const std::string& name, const std::vector<unsigned>& shape ) {
  if( !keywords.outputComponentExists(name,true) ) {
    plumed_merror("a description of component " + name + " has not been added to the manual. Components should be registered like keywords in "
                  "registerKeywords as described in the developer doc.");
  }
  std::string thename; thename=getLabel() + "." + name;
  for(unsigned i=0; i<values.size(); ++i) {
    plumed_massert(values[i]->name!=getLabel(),"Cannot mix single values with components");
    plumed_massert(values[i]->name!=thename,"there is already a value with this name: "+thename);
    plumed_massert(values[i]->name!=thename&&name!="bias","Since PLUMED 2.3 the component 'bias' is automatically added to all biases by the general constructor!\n"
                   "Remove the line addComponentWithDerivatives(\"bias\") from your bias.");
  }
  values.emplace_back(Tools::make_unique<Value>(this,thename, true, shape ) );
  std::string msg="  added component to this action:  "+thename+" \n";
  log.printf(msg.c_str());
}

int ActionWithValue::getComponent( const std::string& name ) const {
  if( !allowComponentsAndValue() ) plumed_massert( !exists( getLabel() ), "You should not be calling this routine if you are using a value");
  std::string thename; thename=getLabel() + "." + name;
  for(unsigned i=0; i<values.size(); ++i) {
    if (values[i]->name==thename) return i;
  }
  plumed_merror("there is no component with name " + name);
  return -1;
}

std::string ActionWithValue::getComponentsList( ) const {
  std::string complist;
  for(unsigned i=0; i<values.size(); ++i) {
    complist+=values[i]->name+" ";
  }
  return complist;
}

std::vector<std::string> ActionWithValue::getComponentsVector( ) const {
  std::vector<std::string> complist;
  for(unsigned i=0; i<values.size(); ++i) {
    complist.push_back(values[i]->name);
  }
  return complist;
}

void ActionWithValue::componentIsNotPeriodic( const std::string& name ) {
  int kk=getComponent(name);
  values[kk]->min=0; values[kk]->max=0;
  values[kk]->setupPeriodicity();
}

void ActionWithValue::componentIsPeriodic( const std::string& name, const std::string& min, const std::string& max ) {
  int kk=getComponent(name);
  values[kk]->setDomain(min,max);
}

void ActionWithValue::setGradientsIfNeeded() {
  if(isOptionOn("GRADIENTS")) {
    for(unsigned i=0; i<values.size(); i++) values[i]->setGradients();
  }
}

void ActionWithValue::turnOnDerivatives() {
  // Turn on the derivatives in all actions on which we are dependent
  for(const auto & p : getDependencies() ) {
    ActionWithValue* vv=dynamic_cast<ActionWithValue*>(p);
    if(vv) vv->turnOnDerivatives();
  }
  // Turn on the derivatives
  noderiv=false;
  // Resize the derivatives
  for(unsigned i=0; i<values.size(); ++i) values[i]->resizeDerivatives( getNumberOfDerivatives() );
}

Value* ActionWithValue::getPntrToOutput( const unsigned& ind ) const {
  plumed_dbg_massert(ind<values.size(),"you have requested a pointer that is out of bounds");
  return values[ind].get();
}

Value* ActionWithValue::getPntrToComponent( const std::string& name ) {
  int kk=getComponent(name);
  return values[kk].get();
}

Value* ActionWithValue::getPntrToComponent( int n ) {
  plumed_dbg_massert(n<values.size(),"you have requested a pointer that is out of bounds");
  return values[n].get();
}

unsigned ActionWithValue::getGridArgumentIndex( const ActionWithArguments* aa ) const { 
  unsigned ng=0; 
  for(unsigned i=0;i<aa->getNumberOfArguments();++i) {
      Value* myarg = aa->getPntrToArgument(i); if( myarg->getRank()==0 ) continue;
      if( myarg->hasDerivatives() || myarg->isTimeSeries() ) ng=i;
  }
  return ng;
}

void ActionWithValue::getInfoForGridHeader( std::string& gtype, std::vector<std::string>& argn, std::vector<std::string>& min,
                                            std::vector<std::string>& max, std::vector<unsigned>& nbin,
                                            std::vector<double>& spacing, std::vector<bool>& pbc, const bool& dumpcube ) const {
  const ActionWithArguments* aa = dynamic_cast<const ActionWithArguments*>( this );
  if( aa ) {
    unsigned ng = getGridArgumentIndex(aa);
    ((aa->getPntrToArgument(ng))->getPntrToAction())->getInfoForGridHeader( gtype, argn, min, max, nbin, spacing, pbc, dumpcube );
    return;
  }
  plumed_merror( "problem in getting grid data for " + getLabel() );
}

void ActionWithValue::getGridPointIndicesAndCoordinates( const unsigned& ind, std::vector<unsigned>& indices, std::vector<double>& coords ) const { 
  const ActionWithArguments* aa = dynamic_cast<const ActionWithArguments*>( this );
  if( aa ) {
    unsigned ng = getGridArgumentIndex(aa); 
    ((aa->getPntrToArgument(ng))->getPntrToAction())->getGridPointIndicesAndCoordinates( ind, indices, coords );
    return;
  }
  plumed_merror("problem in getting grid data for " + getLabel() ); 
}

void ActionWithValue::getGridPointAsCoordinate( const unsigned& ind, const bool& setlength, std::vector<double>& coords ) const { 
  const ActionWithArguments* aa = dynamic_cast<const ActionWithArguments*>( this );
  if( aa ) {
    unsigned ng = getGridArgumentIndex(aa);
    ((aa->getPntrToArgument(ng))->getPntrToAction())->getGridPointAsCoordinate( ind, false, coords );
    if( coords.size()==(getPntrToOutput(0)->getRank()+1) ) coords[getPntrToOutput(0)->getRank()] = getPntrToOutput(0)->get(ind);
    else if( setlength ) { double val=getPntrToOutput(0)->get(ind); for(unsigned i=0; i<coords.size(); ++i) coords[i] = val*coords[i]; }
    return;
  }
  plumed_merror("problem in getting grid data for " + getLabel() ); 
}

void ActionWithValue::interpretDataLabel( const std::string& mystr, Action* myuser, unsigned& nargs, std::vector<Value*>& args ) {
  // Check for streams
  unsigned nstr=0; for(unsigned i=0; i<values.size(); ++i) { if( values[i]->getRank()>0 ) nstr++; }

  if( mystr=="" || mystr==getLabel() ) {
    // Retrieve value with name of action
    if( values[0]->name!=getLabel() ) myuser->error("action " + getLabel() + " does not have a value");
    // args.push_back( values[0] );
    values[0]->interpretDataRequest( myuser->getLabel(), nargs, args, "" );
  } else if( mystr==getLabel() + ".*" ) {
    for(unsigned i=0; i<values.size(); ++i ) values[i]->interpretDataRequest( myuser->getLabel(), nargs, args, "" );
  } else if( mystr.find(".")!=std::string::npos && exists( mystr ) ) {
    // Retrieve value with specific name
    copyOutput( mystr )->interpretDataRequest( myuser->getLabel(), nargs, args, "" );
  } else {
    std::size_t dot1 = mystr.find_first_of('.'); std::string thelab = mystr.substr(0,dot1);
    plumed_assert( thelab==getLabel() ); std::string rest = mystr.substr(dot1+1);
    if( values.size()==1 ) {
      plumed_assert( values[0]->getRank()>0 && values[0]->getName()==getLabel() );
      values[0]->interpretDataRequest( myuser->getLabel(), nargs, args, rest );
    } else {
      std::size_t dot2 = rest.find_first_of('.'); std::string thecomp = rest.substr(0,dot2);
      if( !exists( thelab + "." + thecomp ) ) myuser->error("could not find component with label " + thelab + "." + thecomp );
      getPntrToComponent( thecomp )->interpretDataRequest( myuser->getLabel(), nargs, args, rest.substr(dot2+1) );
    }
  }
}

void ActionWithValue::addTaskToList( const unsigned& taskCode ) {
  if( !thisAsActionWithArguments ) thisAsActionWithArguments = dynamic_cast<const ActionWithArguments*>( this );
  if( !thisAsActionWithVatom ) thisAsActionWithVatom=dynamic_cast<const ActionWithVirtualAtom*>( this );
  fullTaskList.push_back( taskCode ); taskFlags.push_back(0);
  plumed_assert( fullTaskList.size()==taskFlags.size() );
}

bool ActionWithValue::checkUsedOutsideOfChain( const std::vector<std::string>& actionLabelsInChain, const std::string& parent, 
                                               std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags ) {
  bool usedOutsideOfChain=false;
  for(unsigned i=0;i<values.size();++i) {
      if( values[i]->getRank()==0 ) continue;
      // Check for users of full list of values 
      for(const auto & p : values[i]->userdata) {
          // Check if the action is only being used within the chain
          bool inchain=false;
          for(unsigned j=0;j<actionLabelsInChain.size();++j) {
              if( p.first==actionLabelsInChain[j] ){ inchain=true; break; }
          }
          // If the action we are using is not in the chain check if it is safe to deactivate some stuff 
          if( !inchain ) {
              ActionWithArguments* user=plumed.getActionSet().selectWithLabel<ActionWithArguments*>( p.first );
              if( user ) user->getTasksForParent( parent, actionsThatSelectTasks, tflags );
              usedOutsideOfChain=true;
          }
      }
  } 
  return usedOutsideOfChain;
}

void ActionWithValue::selectActiveTasks( const std::vector<std::string>& actionLabelsInChain, bool& forceAllTasks,
                                         std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags ) {
  buildCurrentTaskList( forceAllTasks, actionsThatSelectTasks, tflags );
  // Check which actions are using the values calculated by this action
  bool usedOutsideOfChain=checkUsedOutsideOfChain( actionLabelsInChain, getLabel(), actionsThatSelectTasks, tflags );
  // Now check if we can deactivate tasks with this action by checking if it is one of the actions that 
  // allows deactivated tasks
  if( usedOutsideOfChain ) {
      bool OKToDeactivate=false;
      for(unsigned i=0;i<actionsThatSelectTasks.size();++i) {
          if( getLabel()==actionsThatSelectTasks[i] ){ OKToDeactivate=true; break; }
      }
      if( !OKToDeactivate ) forceAllTasks=true;
  }
  // And now do stuff for the next action in the chain
  if( action_to_do_after ) {
      // Retrieve the atoms for tasks in the stream
      // ActionAtomistic* aa=dynamic_cast<ActionAtomistic*>( action_to_do_after ); 
      if( atom_action_to_do_after ) atom_action_to_do_after->retrieveAtoms();
      // And set up the task list for the next one
      action_to_do_after->selectActiveTasks( actionLabelsInChain, forceAllTasks, actionsThatSelectTasks, tflags );
  }
}

void ActionWithValue::prepareForTaskLoop( const unsigned& nactive, const std::vector<unsigned>& pTaskList ) {
  prepareForTasks( nactive, pTaskList ); 
  if( action_to_do_after ) action_to_do_after->prepareForTaskLoop( nactive, pTaskList );
}

unsigned ActionWithValue::setTaskFlags( std::vector<unsigned>& tflags, std::vector<unsigned>&  pTaskList, std::vector<unsigned>& pIndexList ) {
  if( actionsLabelsInChain.size()==0 ) getAllActionLabelsInChain( actionsLabelsInChain );
  std::vector<std::string> actionsThatSelectTasks; bool forceAllTasks=false;
  tflags.assign(tflags.size(),0); selectActiveTasks( actionsLabelsInChain, forceAllTasks, actionsThatSelectTasks, tflags );
  if( forceAllTasks || actionsThatSelectTasks.size()==0 ) tflags.assign(tflags.size(),1);
  // Now actually build the list of active tasks
  unsigned n_active = 0;
  for(unsigned i=0; i<fullTaskList.size(); ++i) {
    if( tflags[i]>0 ) n_active++;
  }
  // Now create the partial task list
  unsigned n=0; pTaskList.resize( n_active ); pIndexList.resize( n_active );
  for(unsigned i=0; i<fullTaskList.size(); ++i) {
    // Deactivate sets inactive tasks to number not equal to zero
    if( tflags[i]>0 ) { pTaskList[n] = fullTaskList[i]; pIndexList[n]=i; n++; } 
  }
  return n_active;
}

void ActionWithValue::runAllTasks() {
// Skip this if this is done elsewhere
  if( action_to_do_before ) return;

  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial) { stride=1; rank=0; }

  // Get number of threads for OpenMP
  unsigned nt=OpenMP::getNumThreads();
  if( nt*stride*10>nactive_tasks ) nt=nactive_tasks/stride/10;
  if( nt==0 || no_openmp ) nt=1;

  // Determine if some tasks can be switched off
  nactive_tasks = setTaskFlags( taskFlags, partialTaskList, indexOfTaskInFullList );
  // Now do all preparations required to run all the tasks
  prepareForTaskLoop( nactive_tasks, partialTaskList );

  // Get the total number of streamed quantities that we need
  unsigned nquantities = 0, ncols=0, nmatrices=0;
  getNumberOfStreamedQuantities( nquantities, ncols, nmatrices );
  // This needs to be separate from StreamedQuntities as the above
  // has to be const in order to be called from rerunTask
  setupVirtualAtomStashes( nquantities );
  // Get size for buffer
  unsigned bufsize=0; getSizeOfBuffer( nactive_tasks, bufsize );
  if( buffer.size()!=bufsize ) buffer.resize( bufsize );
  // Clear buffer
  buffer.assign( buffer.size(), 0.0 );

  // Recover the number of derivatives we require
  unsigned nderivatives = 0; bool gridsInStream=checkForGrids();
  if( !noderiv || gridsInStream ) getNumberOfStreamedDerivatives( nderivatives );

  if(timers) stopwatch.start("2 Loop over tasks");
  #pragma omp parallel num_threads(nt)
  {
    std::vector<double> omp_buffer;
    if( nt>1 ) omp_buffer.resize( bufsize, 0.0 );
    MultiValue myvals( nquantities, nderivatives, ncols, nmatrices );
    myvals.clearAll();

    #pragma omp for nowait
    for(unsigned i=rank; i<nactive_tasks; i+=stride) {
      // Calculate the stuff in the loop for this action
      runTask( indexOfTaskInFullList[i], partialTaskList[i], myvals );

      // Now transfer the data to the actions that accumulate values from the calculated quantities
      if( nt>1 ) {
        gatherAccumulators( indexOfTaskInFullList[i], myvals, omp_buffer );
      } else {
        gatherAccumulators( indexOfTaskInFullList[i], myvals, buffer );
      }

      // Clear the value
      myvals.clearAll();
    }
    #pragma omp critical
    if(nt>1) for(unsigned i=0; i<bufsize; ++i) buffer[i]+=omp_buffer[i];
  }
  if(timers) stopwatch.stop("2 Loop over tasks");

  if(timers) stopwatch.start("3 MPI gather");
  // MPI Gather everything
  if( !serial && buffer.size()>0 ) comm.Sum( buffer );
  // Update the elements that are makign contributions to the sum here
  // this causes problems if we do it in prepare
  if(timers) stopwatch.stop("3 MPI gather");

  if(timers) stopwatch.start("4 Finishing computations");
  finishComputations( buffer );
  if(timers) stopwatch.stop("4 Finishing computations");
}

bool ActionWithValue::checkForGrids() const {
  for(unsigned i=0; i<values.size(); ++i) {
    if( values[i]->getRank()>0 && values[i]->hasDerivatives() ) return true;
  }
  if( action_to_do_after ) return action_to_do_after->checkForGrids();
  return false;
}

void ActionWithValue::getNumberOfStreamedDerivatives( unsigned& nderivatives ) const {
  unsigned nnd=0;
  if( !noderiv ) {
    nnd = getNumberOfDerivatives();
  } else {
    for(unsigned i=0; i<values.size(); ++i) {
      if( values[i]->getRank()>0 && values[i]->hasDerivatives() ) { nnd = getNumberOfDerivatives(); break; }
    }
  }
  if( nnd>nderivatives ) nderivatives = nnd;
  if( action_to_do_after ) action_to_do_after->getNumberOfStreamedDerivatives( nderivatives );
}

void ActionWithValue::setupVirtualAtomStashes( unsigned& nquants ) {
  ActionWithVirtualAtom* av = dynamic_cast<ActionWithVirtualAtom*>( this );
  if( av ) av->setStashIndices( nquants );
  if( action_to_do_after ) action_to_do_after->setupVirtualAtomStashes( nquants );
}

void ActionWithValue::getNumberOfStreamedQuantities( unsigned& nquants, unsigned& ncols, unsigned& nmat ) const {
  // const ActionWithArguments* aa = dynamic_cast<const ActionWithArguments*>( this );
  if( thisAsActionWithArguments ) thisAsActionWithArguments->getNumberOfStashedInputArguments( nquants );

  for(unsigned i=0; i<values.size(); ++i) {
    if( values[i]->getRank()==2 && !values[i]->hasDerivatives() ) {
      if( values[i]->getShape()[1]>ncols ) { ncols = values[i]->getShape()[1]; }
      values[i]->matpos=nmat; nmat++;
      // Matrices store is reshaped every step as we do a sparse storage of data to minimise memory requirements
      if( values[i]->storedata ) values[i]->reshapeMatrixStore();
    } else if( values[i]->getRank()==1 && values[i]->columnsums ) {
      values[i]->matpos=nmat; nmat++;
    }
    values[i]->streampos=nquants; nquants++;
  }
  if( action_to_do_after ) action_to_do_after->getNumberOfStreamedQuantities( nquants, ncols, nmat );
}

void ActionWithValue::getSizeOfBuffer( const unsigned& nactive_tasks, unsigned& bufsize ) {
  for(unsigned i=0; i<values.size(); ++i) {
    values[i]->bufstart=bufsize; bufsize += values[i]->data.size();
  }
  if( action_to_do_after ) action_to_do_after->getSizeOfBuffer( nactive_tasks, bufsize );
}

void ActionWithValue::runTask( const std::string& controller, const unsigned& task_index, const unsigned& current, const unsigned colno, MultiValue& myvals ) const {
  // Do matrix element task
  bool wasperformed=false; myvals.setTaskIndex(task_index); myvals.setSecondTaskIndex( colno );
  if( isActive() ) wasperformed=performTask( controller, current, colno, myvals );
  // const ActionWithArguments* aa = dynamic_cast<const ActionWithArguments*>( this );
  if( thisAsActionWithArguments ) {
    if( actionInChain() && thisAsActionWithArguments->getNumberOfArguments()>0 ) {
      // Now check if the task takes a matrix as input - if it does do it
      bool do_this_task = ((thisAsActionWithArguments->getPntrToArgument(0))->getRank()==2 && !(thisAsActionWithArguments->getPntrToArgument(0))->hasDerivatives() );
#ifdef DNDEBUG
      if( do_this_task ) {
        for(unsigned i=1; i<thisAsActionWithArguments->getNumberOfArguments(); ++i) {
          plumed_dbg_assert( (thisAsActionWithArguments->getPntrToArgument(i))->getRank()==2 && !(thisAsActionWithArguments->getPntrToArgument(0))->hasDerivatives() );
        }
      }
#endif
      if( do_this_task && isActive() ) { myvals.vector_call=false; myvals.setTaskIndex(task_index); performTask( current, myvals ); }
    }
  }

  // Check if we need to store stuff
  bool matrix=wasperformed;
  for(unsigned i=0; i<values.size(); ++i) {
    if( values[i]->getRank()!=2 || values[i]->hasDerivatives() ) { matrix=false; break; }
  }
  if( matrix ) {
    unsigned col_stash_index = colno; if( colno>=getFullNumberOfTasks() ) col_stash_index = colno - getFullNumberOfTasks();
    for(unsigned i=0; i<values.size(); ++i) {
      unsigned matindex = values[i]->getPositionInMatrixStash();;
      if( values[i]->hasForce ) {
        unsigned sind = values[i]->streampos, find = col_stash_index;
        if( values[i]->getNumberOfColumns()<values[i]->shape[1] ) find=myvals.getNumberOfStashedMatrixElements(matindex);
        double fforce = values[i]->getForce( myvals.getTaskIndex()*getNumberOfColumns() + find );
        for(unsigned j=0; j<myvals.getNumberActive(sind); ++j) {
          unsigned kindex = myvals.getActiveIndex(sind,j); myvals.addMatrixForce( matindex, kindex, fforce*myvals.getDerivative(sind,kindex ) );
        }
      } 
      if( values[i]->storedata ) myvals.stashMatrixElement( matindex, col_stash_index, myvals.get( values[i]->getPositionInStream() ) );
    }
  }

  // Now continue on with the stream
  if( action_to_do_after ) action_to_do_after->runTask( controller, task_index, current, colno, myvals );
}

void ActionWithValue::runTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const {
  if( isActive() ) {
    myvals.setTaskIndex(task_index); myvals.vector_call=true; performTask( current, myvals );
  }
  if( action_to_do_after ) action_to_do_after->runTask( task_index, current, myvals );
}

void ActionWithValue::clearMatrixElements( MultiValue& myvals ) const {
  if( isActive() ) {
    for(unsigned i=0; i<values.size(); ++i) {
      if( values[i]->getRank()==2 && !values[i]->hasDerivatives() ) myvals.clear( values[i]->getPositionInStream() );
    }
  }
  if( action_to_do_after ) action_to_do_after->clearMatrixElements( myvals );
}

void ActionWithValue::rerunTask( const unsigned& task_index, MultiValue& myvals ) const {
  if( !action_to_do_before ) {
    unsigned ncol=0, nmat=0, nquantities = 0; getNumberOfStreamedQuantities( nquantities, ncol, nmat );
    unsigned nderivatives = 0; getNumberOfStreamedDerivatives( nderivatives );
    if( myvals.getNumberOfValues()!=nquantities || myvals.getNumberOfDerivatives()!=nderivatives ) myvals.resize( nquantities, nderivatives, ncol, nmat );
    runTask( task_index, fullTaskList[task_index], myvals ); return;
  }
  action_to_do_before->rerunTask( task_index, myvals );

}

void ActionWithValue::gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals,
                                         const unsigned& bufstart, std::vector<double>& buffer ) const {
  plumed_dbg_assert( !values[valindex]->hasDeriv && values[valindex]->getRank()>0 );
  unsigned sind = values[valindex]->streampos;
  // This looks after storing for matrices
  if( values[valindex]->getRank()==2 && !values[valindex]->alwaysstore ) {
    unsigned rstart = (1+values[valindex]->getNumberOfColumns())*code;
    unsigned matind = values[valindex]->getPositionInMatrixStash();
    unsigned vindex = bufstart + code*values[valindex]->getNumberOfColumns();
    // Only matrix elements that are definitely non-zero are stored here
    // This is used with contact matrices to reduce the ammount of memory required
    if( values[valindex]->getNumberOfColumns()<values[valindex]->getShape()[1] ) {
        values[valindex]->matindexes[rstart]=0; 
        for(unsigned j=0; j<myvals.getNumberOfStashedMatrixElements(matind); ++j) {
          unsigned jind = myvals.getStashedMatrixIndex(matind,j); 
          values[valindex]->matindexes[rstart+1+values[valindex]->matindexes[rstart]]=jind; values[valindex]->matindexes[rstart]++;
          plumed_dbg_massert( vindex+j<buffer.size(), "failing in " + getLabel() + " on value " + values[valindex]->getName() );
          buffer[vindex + j] += myvals.getStashedMatrixElement( matind, jind );
        }
    // This stores all matrix elements in the expected places for when the matrix is dense
    } else {
        values[valindex]->matindexes[rstart]=values[valindex]->getShape()[1]; unsigned k=rstart+1;
        for(unsigned j=0;j<values[valindex]->matindexes[rstart];++j) { values[valindex]->matindexes[k]=j; k++; }
        for(unsigned j=0; j<myvals.getNumberOfStashedMatrixElements(matind); ++j) {
           unsigned jind = myvals.getStashedMatrixIndex(matind,j);
           plumed_dbg_massert( vindex+jind<buffer.size(), "failing in " + getLabel() + " on value " + values[valindex]->getName() );
           buffer[vindex + jind] += myvals.getStashedMatrixElement( matind, jind );    
        } 
    } 
    // This looks after sums over columns of matrix
  } else if ( values[valindex]->getRank()==1 && values[valindex]->columnsums ) {
    unsigned vindex = bufstart; unsigned matind = values[valindex]->getPositionInMatrixStash();
    for(unsigned j=0; j<myvals.getNumberOfStashedMatrixElements(matind); ++j) {
      unsigned jind = myvals.getStashedMatrixIndex(matind,j);
      buffer[vindex + jind] += myvals.getStashedMatrixElement( matind, jind );
    } 
 // This looks after storing in all other cases
  } else if( getFullNumberOfTasks()==values[valindex]->getNumberOfValues(getLabel()) ) {
    unsigned nspace=1; if( values[valindex]->hasDeriv ) nspace=(1 + values[valindex]->getNumberOfDerivatives() );
    unsigned vindex = bufstart + code*nspace; plumed_dbg_massert( vindex<buffer.size(), "failing in " + getLabel() );
    buffer[vindex] += myvals.get(sind);
  }
} 

void ActionWithValue::gatherAccumulators( const unsigned& taskCode, const MultiValue& myvals, std::vector<double>& buffer ) const {
  if( isActive() ) {
    for(unsigned i=0; i<values.size(); ++i) {
      unsigned bufstart = values[i]->bufstart;
      if( values[i]->getRank()==0 ) {
        plumed_dbg_massert( bufstart<buffer.size(), "problem in " + getLabel() );
        unsigned sind = values[i]->streampos; buffer[bufstart] += myvals.get(sind);
        if( values[i]->hasDerivatives() ) {
          unsigned ndmax = (values[i]->getPntrToAction())->getNumberOfDerivatives();
          for(unsigned k=0; k<myvals.getNumberActive(sind); ++k) {
            unsigned kindex = myvals.getActiveIndex(sind,k);
            plumed_dbg_massert( bufstart+1+kindex<buffer.size(), "problem in " + getLabel()  );
            buffer[bufstart + 1 + kindex] += myvals.getDerivative(sind,kindex);
          }
        }
        // This looks after storing of data
      } else if( values[i]->storedata ) {
        gatherStoredValue( i, taskCode, myvals, bufstart, buffer ); 
      }
    }
    // Special method for dealing with centers
    // const ActionWithVirtualAtom* av = dynamic_cast<const ActionWithVirtualAtom*>( this );
    if(thisAsActionWithVatom) thisAsActionWithVatom->gatherForVirtualAtom( myvals, buffer );
  }

  if( action_to_do_after ) action_to_do_after->gatherAccumulators( taskCode, myvals, buffer );
}

void ActionWithValue::retrieveAllScalarValuesInLoop( const std::string& ulab, unsigned& nargs, std::vector<Value*>& myvals ) {
  for(unsigned i=0; i<values.size(); ++i) {
    if( values[i]->getRank()==0 ) {
      bool found=false;
      for(unsigned j=0; j<myvals.size(); ++j) {
        if( values[i]->getName()==myvals[j]->getName() ) { found=true; break; }
      }
      if( !found ) values[i]->interpretDataRequest( ulab, nargs, myvals, "" );
    }
  }
  if( action_to_do_after ) action_to_do_after->retrieveAllScalarValuesInLoop( ulab, nargs, myvals );
}

void ActionWithValue::finishComputations( const std::vector<double>& buffer ) {
  if( isActive() ) {
    for(unsigned i=0; i<values.size(); ++i) {
      unsigned bufstart = values[i]->bufstart;
      if( values[i]->reset ) values[i]->data.assign( values[i]->data.size(), 0 );
      if( (values[i]->getRank()>0 && values[i]->hasDerivatives()) || values[i]->storedata ) {
        unsigned sz_v = values[i]->data.size();
        for(unsigned j=0; j<sz_v; ++j) values[i]->add( j, buffer[bufstart+j] );
        if( values[i]->getRank()==2 && !values[i]->hasDerivatives() ) comm.Sum( values[i]->matindexes );
      }
      if( !doNotCalculateDerivatives() && values[i]->hasDeriv && values[i]->getRank()==0 ) {
        for(unsigned j=0; j<values[i]->getNumberOfDerivatives(); ++j) values[i]->setDerivative( j, buffer[bufstart+1+j] );
      }
    }
    transformFinalValueAndDerivatives( buffer );
  }
  if( action_to_do_after ) action_to_do_after->finishComputations( buffer );
}

bool ActionWithValue::getForcesFromValues( std::vector<double>& forces ) {
  unsigned type=0;
  if( values[0]->shape.size()==0 ) type=1;
  else if( values[0]->hasDeriv ) type=2;
  else plumed_assert( values[0]->shape.size()>0 );

#ifndef DNDEBUG
  if( type==0 ) {
    for(unsigned i=0; i<values.size(); ++i) plumed_dbg_assert( values[i]->shape.size()>0 && !values[i]->hasDeriv );
  } else if( type==1 && getName()!="DIAGONALIZE" ) {
    for(unsigned i=0; i<values.size(); ++i) plumed_dbg_assert( values[i]->shape.size()==0 );
  } else if( type==2 ) {
    for(unsigned i=0; i<values.size(); ++i) plumed_dbg_assert( values[i]->shape.size()>0 && values[i]->hasDeriv );
  } else if( getName()!="DIAGONALIZE" ) {
    plumed_merror("value type not defined");
  }
#endif
  bool at_least_one_forced=false;
  if( type==1 ) {
    for(unsigned i=0; i<values.size(); ++i) {
      if( values[i]->applyForce( forces ) ) at_least_one_forced=true;
    }
  } else {
    // Check if there are any forces
    for(unsigned i=0; i<values.size(); ++i) {
      if( values[i]->hasForce ) at_least_one_forced=true;
    }
    if( !at_least_one_forced ) return false;

    // Get the action that calculates these quantitites
    ActionWithValue* av = getActionThatCalculates();
    nactive_tasks = av->nactive_tasks;
    // Setup MPI parallel loop
    unsigned stride=comm.Get_size();
    unsigned rank=comm.Get_rank();
    if(serial) { stride=1; rank=0; }

    // Get number of threads for OpenMP
    unsigned nt=OpenMP::getNumThreads();
    if( nt*stride*10>nactive_tasks ) nt=nactive_tasks/stride/10;
    if( nt==0 || no_openmp ) nt=1;

    // Now determine how big the multivalue needs to be
    unsigned nderiv=0; av->getNumberOfStreamedDerivatives( nderiv );
    unsigned nquants=0, ncols=0, nmatrices=0; 
    av->getNumberOfStreamedQuantities( nquants, ncols, nmatrices );
    av->setupVirtualAtomStashes( nquants );
    #pragma omp parallel num_threads(nt)
    {
      std::vector<double> omp_forces;
      if( nt>1 ) omp_forces.resize( forces.size(), 0.0 );
      MultiValue myvals( nquants, nderiv, ncols, nmatrices, forces.size() );
      myvals.clearAll();

      #pragma omp for nowait
      for(unsigned i=rank; i<nactive_tasks; i+=stride) {
        unsigned itask = av->indexOfTaskInFullList[i];
        av->runTask( itask, av->partialTaskList[i], myvals );

        // Now get the forces
        if( nt>1 ) av->gatherForces( av->partialTaskList[i], myvals, omp_forces );
        else av->gatherForces( av->partialTaskList[i], myvals, forces );

        myvals.clearAll();
        myvals.clearStoredForces();
      }
      #pragma omp critical
      if(nt>1) for(unsigned i=0; i<forces.size(); ++i) forces[i]+=omp_forces[i];
    }
    // MPI Gather on forces
    if( !serial ) comm.Sum( forces );
    // And clear all the forces that have been added in one shot
    av->clearAllForcesInChain();
  }

  return at_least_one_forced;
}

void ActionWithValue::gatherForces( const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const {
  bool hasforce=false;
  for(unsigned i=0; i<values.size(); ++i) {
    if( values[i]->hasForce ) { hasforce=true; break; }
  }
  if( hasforce && isActive() ) {
    for(unsigned k=0; k<values.size(); ++k) {
      if( values[k]->hasForce && values[k]->getRank()==2 && !values[k]->hasDeriv ) {
        unsigned matind = values[k]->getPositionInMatrixStash();
        for(unsigned j=0; j<forces.size(); ++j) forces[j] += myvals.getStashedMatrixForce( matind, j );
      } else if( values[k]->getRank()>0 && values[k]->hasForce ) {
        unsigned sspos = values[k]->streampos; double fforce = values[k]->getForce(itask);
        for(unsigned j=0; j<myvals.getNumberActive(sspos); ++j) {
          unsigned jder=myvals.getActiveIndex(sspos, j); forces[jder] += fforce*myvals.getDerivative( sspos, jder );
        }
      }
    }
  }
  if( action_to_do_after ) action_to_do_after->gatherForces( itask, myvals, forces );
}

void ActionWithValue::clearAllForcesInChain() {
  clearInputForces(); if( action_to_do_after ) action_to_do_after->clearAllForcesInChain();
}

std::string ActionWithValue::getCleanGraphLabel( const std::string& gilab ) {
  std::string glab = gilab; 
  if( glab.find("@")!=std::string::npos ){ std::size_t at=glab.find_first_of("@"); glab=glab.substr(at+1); }
  for(unsigned j=0;;++j) {
      std::size_t pos = glab.find_first_of("-"); if( pos==std::string::npos ) break;
      glab = glab.substr(0,pos) + "h" + glab.substr(pos+1);
  }
  for(unsigned j=0;;++j) {
     std::size_t pos = glab.find_first_of("["); if( pos==std::string::npos ) break;
     glab = glab.substr(0,pos) + glab.substr(pos+1);
  }
  for(unsigned j=0;;++j) {
     std::size_t pos = glab.find_first_of("]"); if( pos==std::string::npos ) break;
     glab = glab.substr(0,pos) + glab.substr(pos+1);
  }
  return glab;
}

void ActionWithValue::generateGraphNodes( OFile& ofile, std::vector<std::string>& graph_actions ) const {
  if( action_to_do_before ) return ;
  const ActionWithVirtualAtom* ava=dynamic_cast<const ActionWithVirtualAtom*>(this);
  if( ava ) return ;

  // Check we have not dealt with this node already
  for(unsigned i=0;i<graph_actions.size();++i) {
      if( graph_actions[i]==getLabel() ) return;
  }

  // Now retrieve all the labels in this chain
  std::vector<std::string> chain; getAllActionLabelsInChain( chain );
  if( chain.size()>1 ) { 
      ofile.printf("   subgraph cluster%d { \n", graph_actions.size() );
      ofile.printf("      node [style=filled,fillcolor=lightgrey];\n");
      ofile.printf("      penwidth=3;\n");
      ofile.printf("      color=black;\n");
      // Now create all the nodes in this chain
      unsigned gstart = graph_actions.size(); 
      for(unsigned i=0;i<chain.size();++i){
          ActionWithValue* av=plumed.getActionSet().selectWithLabel<ActionWithValue*>( chain[i] ); std::string exline, num; 
          std::string label=getCleanGraphLabel( av->getLabel() );
          if( av->writeInGraph(exline) ) {
              ofile.printf("     %s [label=\"%d \\n %s: \\n %s \\n %s\"] \n", label.c_str(), graph_actions.size() - gstart +1, av->getLabel().c_str(), av->getName().c_str(), exline.c_str() );
          } else {
              ofile.printf("     %s [label=\"%d \\n %s: \\n %s\"] \n", label.c_str(), graph_actions.size() - gstart +1, av->getLabel().c_str(), av->getName().c_str() );   
          }
          for(unsigned j=0;j<av->values.size();++j) {
              for(const auto & p : (av->values[j])->userdata) {
                  // Check if the action is only being used within the chain
                  bool inchain=false;
                  for(unsigned k=0;k<chain.size();++k) {
                      if( p.first==chain[k] ){ inchain=true; break; }
                  }
                  if( inchain ) {
                      ActionWithValue* av2=plumed.getActionSet().selectWithLabel<ActionWithValue*>( p.first ); plumed_assert( av2 );
                      std::string label2=getCleanGraphLabel( av2->getLabel() ); 
                      std::string color="orange";
                      if( (av->values[j])->getRank()>0 && (av->values[j])->hasDerivatives() ) color="green";
                      else if( (av->values[j])->getRank()==2 ) color="red";
                      else if( (av->values[j])->getRank()==1 ) color="blue"; 
                      ofile.printf("     %s -> %s [label=\"%s\", color=%s, fontcolor=%s]; \n", label.c_str(), label2.c_str(), 
                                                               (av->values[j])->getName().c_str(), color.c_str(), color.c_str() ); 
                  }
              }
          }
          graph_actions.push_back( av->getLabel() );
      }
      ofile.printf("   }\n");
  } else {
      std::string label=getCleanGraphLabel( getLabel() ); 
      const ActionWithValue* av = dynamic_cast<const ActionWithValue*>( this ); std::string exline; graph_actions.push_back( getLabel() );
      if( av ) {  
          if( av->writeInGraph(exline) ) {
              ofile.printf("     %s [label=\"%s: \\n %s \\n %s\"] \n", label.c_str(), getLabel().c_str(), getName().c_str(), exline.c_str() );
          } else {
              ofile.printf("     %s [label=\"%s: \\n %s\"] \n", label.c_str(), getLabel().c_str(), getName().c_str() );
          }
      } else ofile.printf("     %s [label=\"%s: \\n %s\"] \n", label.c_str(), getLabel().c_str(), getName().c_str() ); 
  }
  // Now create the links to the nodes outside of this chain
  for(unsigned i=0;i<chain.size();++i) {
      ActionWithValue* av=plumed.getActionSet().selectWithLabel<ActionWithValue*>( chain[i] ); std::string label=getCleanGraphLabel( av->getLabel() ); 
      for(unsigned j=0;j<av->values.size();++j) {
          for(const auto & p : (av->values[j])->userdata) {
              // Check if the action is only being used within the chain
              bool inchain=false;
              for(unsigned k=0;k<chain.size();++k) {
                  if( p.first==chain[k] ){ inchain=true; break; }
              }
              if( !inchain ) {
                  Action* av2=plumed.getActionSet().selectWithLabel<Action*>( p.first ); plumed_assert( av2 );
                  ActionWithValue* avv2 = dynamic_cast<ActionWithValue*>( av2 );
                  if( !avv2 ) {
                      bool found=false;
                      for(unsigned i=0;i<graph_actions.size();++i) {
                          if( av2->getLabel()==graph_actions[i] ){ found=true; break; }
                      } 
                      if( !found ) {
                          std::string label=getCleanGraphLabel( av2->getLabel() ); 
                          ofile.printf("     %s [label=\"%s: \\n %s\"] \n", label.c_str(), av2->getLabel().c_str(), av2->getName().c_str() );
                          graph_actions.push_back( av2->getLabel() );
                      }
                  }
                  std::string label2=getCleanGraphLabel( av2->getLabel() ); std::string color="orange";
                  if( (av->values[j])->getRank()>0 && (av->values[j])->hasDerivatives() ) color="green";
                  else if( (av->values[j])->getRank()==2 ) color="red";
                  else if( (av->values[j])->getRank()==1 ) color="blue";
                  ofile.printf("     %s -> %s [label=\"%s\", color=%s, fontcolor=%s]; \n", label.c_str(), label2.c_str(), 
                                                               (av->values[j])->getName().c_str(), color.c_str(), color.c_str() ); 
              }
          }
      }
  }
}

}
