/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "ActionWithVector.h"
#include "ActionWithMatrix.h"
#include "PlumedMain.h"
#include "ActionSet.h"
#include "tools/OpenMP.h"
#include "tools/Communicator.h"

namespace PLMD {

enum class Option { no, yes };

Option interpretEnvString(const char* env,const char* str) {
  if(!str) return Option::yes;
  if(!std::strcmp(str,"yes"))return Option::yes;
  if(!std::strcmp(str,"no"))return Option::no;
  plumed_error()<<"Cannot understand env var "<<env<<"\nPossible values: yes/no\nActual value: "<<str;
}

/// Switch on/off chains of actions using PLUMED environment variable
/// export PLUMED_FORBID_CHAINS=yes  # forbid the use of chains in this run
/// export PLUMED_FORBID_CHAINS=no   # allow chains to be used in the run
/// default: yes
Option getenvChainForbidden() {
  static const char* name="PLUMED_FORBID_CHAINS";
  static const auto opt = interpretEnvString(name,std::getenv(name));
  return opt;
}

void ActionWithVector::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys ); keys.remove("NUMERICAL_DERIVATIVES");
  ActionWithArguments::registerKeywords( keys );
  keys.addFlag("SERIAL",false,"do the calculation in serial.  Do not parallelize");
  keys.reserve("optional","MASK","the label for a sparse matrix that should be used to determine which elements of the matrix should be computed");
}

ActionWithVector::ActionWithVector(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  nmask(-1),
  serial(false),
  forwardPass(false),
  action_to_do_before(NULL),
  action_to_do_after(NULL),
  never_reduce_tasks(false),
  reduce_tasks(false),
  atomsWereRetrieved(false)
{
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    ActionWithVector* av = dynamic_cast<ActionWithVector*>( getPntrToArgument(i)->getPntrToAction() );
    if( av && av->getNumberOfMasks()>=0 ) nmask=0;
  }

  if( keywords.exists("SERIAL") ) parseFlag("SERIAL",serial);
  if( keywords.exists("MASK") ) {
    std::vector<Value*> mask; parseArgumentList("MASK",mask);
    if( mask.size()>0 ) {
      if( nmask>=0 && getNumberOfArguments()==1 ) {
        ActionWithVector* av=dynamic_cast<ActionWithVector*>( getPntrToArgument(0)->getPntrToAction() );
        plumed_massert( av, "input should be a vector from ActionWithVector" ); unsigned j=0, nargs = av->getNumberOfArguments();
        for(unsigned i=nargs-av->nmask; i<nargs; ++i) {
          if( av->getPntrToArgument(i)!=mask[j] ) error("the masks in subsequent actions do not match");
          j++;
        }
      }
      if( getNumberOfArguments()>0 && getPntrToArgument(0)->hasDerivatives() ) error("input for mask should be vector or matrix");
      else if( mask[0]->getRank()==2 ) {
        if( mask.size()>1 ) error("MASK should only have one argument");
        log.printf("  only computing elements of matrix that correspond to non-zero elements of matrix %s \n", mask[0]->getName().c_str() );
      } else if( mask[0]->getRank()==1 ) {
        log.printf("  only computing elements of vector that correspond to non-zero elements of vectors %s", mask[0]->getName().c_str() );
        for(unsigned i=1; i<mask.size(); ++i) {
          if( mask[i]->getRank()!=1 ) { log.printf("\n"); error("input to mask should be vector"); }
          log.printf(", %s", mask[i]->getName().c_str() );
        }
        log.printf("\n");
      }
      std::vector<Value*> allargs( getArguments() ); nmask=mask.size();
      for(unsigned i=0; i<mask.size(); ++i) allargs.push_back( mask[i] );
      requestArguments( allargs );
    }
  }
}

ActionWithVector::~ActionWithVector() {
  if( action_to_do_before ) action_to_do_before->action_to_do_after=NULL;
}

void ActionWithVector::lockRequests() {
  ActionAtomistic::lockRequests();
  ActionWithArguments::lockRequests();
}

void ActionWithVector::unlockRequests() {
  ActionAtomistic::unlockRequests();
  ActionWithArguments::unlockRequests();
}

void ActionWithVector::calculateNumericalDerivatives(ActionWithValue* av) {
  plumed_merror("cannot calculate numerical derivative for action " + getName() + " with label " + getLabel() );
}

void ActionWithVector::clearDerivatives( const bool& force ) {
  if( !force && actionInChain() ) return;
  ActionWithValue::clearDerivatives();
  if( action_to_do_after ) action_to_do_after->clearDerivatives( true );
}

void ActionWithVector::clearInputForces( const bool& force ) {
  if( !force && actionInChain() ) return;
  ActionWithValue::clearInputForces();
  if( action_to_do_after ) action_to_do_after->clearInputForces( true );
}

const ActionWithVector* ActionWithVector::getFirstActionInChain() const {
  if( !actionInChain() ) return this;
  return action_to_do_before->getFirstActionInChain();
}

ActionWithVector* ActionWithVector::getFirstActionInChain() {
  if( !actionInChain() ) return this;
  return action_to_do_before->getFirstActionInChain();
}

void ActionWithVector::retrieveAtoms( const bool& force ) {
  if( !force && actionInChain() || atomsWereRetrieved ) return;
  ActionAtomistic::retrieveAtoms(); atomsWereRetrieved = !actionInChain();
  if( action_to_do_after ) action_to_do_after->retrieveAtoms( true );
}

bool ActionWithVector::hasStoredArguments() const {
  std::string headstr=getFirstActionInChain()->getLabel();
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( !getPntrToArgument(i)->ignoreStoredValue(headstr) ) return true;
  }
  return false;
}

bool ActionWithVector::argumentDependsOn( const std::string& headstr, ActionWithVector* faction, Value* thearg ) {
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( this!=faction && thearg==getPntrToArgument(i) ) return true;
    ActionWithVector* av = dynamic_cast<ActionWithVector*>( getPntrToArgument(i)->getPntrToAction() );
    if( av && (av->getFirstActionInChain())->getLabel()==headstr ) {
      if( av->argumentDependsOn( headstr, faction, thearg ) ) return true;;
    }
  }
  return false;
}

unsigned ActionWithVector::buildArgumentStore( const unsigned& argstart ) {
  return reallyBuildArgumentStore( argstart );
}

unsigned ActionWithVector::reallyBuildArgumentStore( const unsigned& argstart ) {
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) { if( getPntrToArgument(i)->getRank()>0 ) getPntrToArgument(i)->buildDataStore(); }
  unsigned nder=0; arg_deriv_starts.resize( getNumberOfArguments() );
  for(unsigned i=0; i<getNumberOfArguments(); ++i) { arg_deriv_starts[i] = nder; nder += getPntrToArgument(i)->getNumberOfValues(); }
  return nder;
}

ActionWithVector* ActionWithVector::getActionWithDerivatives( ActionWithVector* depaction ) {
  if( depaction==this || depaction->checkForDependency(this) ) {
    if( getNumberOfAtoms()>0 ) return this;
    std::string c=getFirstActionInChain()->getLabel();
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      if( !getPntrToArgument(i)->ignoreStoredValue(c) && !getPntrToArgument(i)->isConstant() ) return this;
    }
  }
  plumed_assert( action_to_do_before );
  return action_to_do_before->getActionWithDerivatives(depaction);
}

bool ActionWithVector::addActionToChain( const std::vector<std::string>& alabels, ActionWithVector* act ) {
  if( action_to_do_after ) { bool state=action_to_do_after->addActionToChain( alabels, act ); return state; }
  // Check action is not already in chain
  std::vector<std::string> mylabels; getFirstActionInChain()->getAllActionLabelsInChain( mylabels );
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
      plumed_massert( av, "could not cast " + alabels[i] ); bool storingall=true;
      for(int j=0; j<av->getNumberOfComponents(); ++j) {
        if( !(av->getPntrToComponent(j))->storedata ) storingall=false;
      }
      if( !storingall ) return false;
    }
  }
  // This checks that there is nothing that will cause problems in the chain
  mylabels.resize(0); getFirstActionInChain()->getAllActionLabelsInChain( mylabels );
  for(unsigned i=0; i<mylabels.size(); ++i) {
    ActionWithVector* av1=plumed.getActionSet().selectWithLabel<ActionWithVector*>( mylabels[i] );
    for(unsigned j=0; j<i; ++j) {
      ActionWithVector* av2=plumed.getActionSet().selectWithLabel<ActionWithVector*>( mylabels[j] );
      if( !av1->canBeAfterInChain( av2 ) ) error("must calculate " + mylabels[j] + " before " + mylabels[i] );
    }
  }
  action_to_do_after=act; act->action_to_do_before=this; updateTaskListReductionStatus();
  ActionWithVector* head = getFirstActionInChain();
  head->broadcastThatTasksAreReduced( head ); head->finishChainBuild( act );
  return true;
}

void ActionWithVector::updateTaskListReductionStatus() {
  ActionWithVector* head = getFirstActionInChain();
  std::vector<ActionWithVector*> task_reducing_actions; head->canReduceTasks( task_reducing_actions );
  if( task_reducing_actions.size()>0 ) head->reduce_tasks=true;
}

void ActionWithVector::broadcastThatTasksAreReduced( ActionWithVector* aselect ) {
  std::string c=getFirstActionInChain()->getLabel();
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( !getPntrToArgument(i)->ignoreStoredValue(c) ) {
      ActionWithVector* av = dynamic_cast<ActionWithVector*>( getPntrToArgument(i)->getPntrToAction() );
      if( av ) {
        bool found=false;
        ActionWithVector* av_head = av->getFirstActionInChain();
        for(unsigned i=0; i<av_head->task_control_list.size(); ++i) {
          if( aselect==av_head->task_control_list[i] ) { found=true; break; }
        }
        if( !found ) av_head->task_control_list.insert( av_head->task_control_list.begin(), aselect );

        av_head->reduce_tasks=true; av_head->updateTaskReductionFlag( av_head->reduce_tasks );
      }
    }
  }
  if( action_to_do_after ) action_to_do_after->broadcastThatTasksAreReduced( aselect );
}

void ActionWithVector::updateTaskReductionFlag( bool& head_reduce_tasks ) {
  if( actionInChain() ) {
    plumed_assert( task_control_list.size()==0 );
  } else {
    for(unsigned i=0; i<task_control_list.size(); ++i) {
      if( !(task_control_list[i]->getFirstActionInChain())->reduce_tasks ) head_reduce_tasks=false;
    }
  }
  broadcastThatTasksAreReduced( getFirstActionInChain() );
  if( action_to_do_after ) action_to_do_after->updateTaskReductionFlag( head_reduce_tasks );
}

void ActionWithVector::canReduceTasks( std::vector<ActionWithVector*>& task_reducing_actions ) {
  areAllTasksRequired( task_reducing_actions );
  if( action_to_do_after ) action_to_do_after->canReduceTasks( task_reducing_actions );
}

void ActionWithVector::finishChainBuild( ActionWithVector* act ) {
  if( action_to_do_after ) action_to_do_after->finishChainBuild( act );
}

void ActionWithVector::getAllActionLabelsInChain( std::vector<std::string>& mylabels ) const {
  bool found = false ;
  for(unsigned i=0; i<mylabels.size(); ++i) {
    if( getLabel()==mylabels[i] ) { found=true; }
  }
  if( !found ) mylabels.push_back( getLabel() );
  if( action_to_do_after ) action_to_do_after->getAllActionLabelsInChain( mylabels );
}

void ActionWithVector::taskIsActive( const unsigned& current, int& flag ) const {
  flag = checkTaskStatus( current, flag );
  if( flag<=0 && action_to_do_after ) action_to_do_after->taskIsActive( current, flag );
}

void ActionWithVector::getAdditionalTasksRequired( ActionWithVector* action, std::vector<unsigned>& atasks ) {
  for(unsigned i=0; i<task_control_list.size(); ++i ) task_control_list[i]->getAdditionalTasksRequired( action, atasks );
}

void ActionWithVector::prepare() {
  active_tasks.resize(0); atomsWereRetrieved=false;
}

int ActionWithVector::checkTaskIsActive( const unsigned& itask ) const {
  unsigned nargs = getNumberOfArguments();
  if( nargs==0 ) {
    return 1;
  } else if( nmask>0 ) {
    for(unsigned j=nargs-nmask; j<nargs; ++j) {
      Value* myarg = getPntrToArgument(j);
      if( myarg->getRank()==1 && !myarg->hasDerivatives() ) {
        if( fabs(myarg->get(itask))>0 ) return 1;
      } else if( myarg->getRank()==2 && !myarg->hasDerivatives() ) {
        unsigned ncol = myarg->getRowLength(itask);
        unsigned base = itask*myarg->getNumberOfColumns();
        for(unsigned k=0; k<ncol; ++k) {
          if( fabs(myarg->get(base+k,false))>0 ) return 1;
        }
      } else plumed_merror("only matrices and vectors should be used as masks");
    }
  } else {
    for(unsigned i=0; i<nargs; ++i) {
      if( getName()=="OUTER_PRODUCT" && i>0 ) return -1;

      Value* myarg = getPntrToArgument(i);
      if( !myarg->isDerivativeZeroWhenValueIsZero() ) return 1;

      if( myarg->getRank()==0 ) {
        return 1;
      } else if( myarg->getRank()==1 && !myarg->hasDerivatives() ) {
        if( fabs(myarg->get(itask))>0 ) return 1;
      } else if( myarg->getRank()==2 && !myarg->hasDerivatives() ) {
        unsigned ncol = myarg->getRowLength(itask);
        unsigned base = itask*myarg->getNumberOfColumns();
        for(unsigned k=0; k<ncol; ++k) {
          if( fabs(myarg->get(base+k,false))>0 ) return 1;
        }
      } else if( myarg->getRank()>0 ) {
        return 1;
      } else plumed_merror("should not be in action " + getName() );
    }
  }
  return -1;
}

std::vector<unsigned>& ActionWithVector::getListOfActiveTasks( ActionWithVector* action ) {
  if( active_tasks.size()>0 ) return active_tasks;
  unsigned ntasks=0; getNumberOfTasks( ntasks );

  if( getenvChainForbidden()==Option::yes ) {
    std::vector<int> taskFlags( ntasks, -1 );
    for(unsigned i=0; i<ntasks; ++i) taskFlags[i] = checkTaskIsActive(i);
    unsigned nt=0;
    for(unsigned i=0; i<ntasks; ++i) {
      if( taskFlags[i]>0 ) nt++;
    }
    active_tasks.resize(nt); nt=0;
    for(unsigned i=0; i<ntasks; ++i) {
      if( taskFlags[i]>0 ) { active_tasks[nt]=i; nt++; }
    }
    return active_tasks;
  }

  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial) { stride=1; rank=0; }

  // Get number of threads for OpenMP
  unsigned nt=OpenMP::getNumThreads();
  if( nt*stride*10>ntasks ) nt=ntasks/stride/10;
  if( nt==0 ) nt=1;

  if( !never_reduce_tasks && reduce_tasks ) {
    if( task_control_list.size()>0 ) {
      // Get the list of tasks that are active in the action that uses the output of this action
      for(unsigned i=0; i<task_control_list.size(); ++i) {
        task_control_list[i]->retrieveAtoms();
        active_tasks = task_control_list[i]->getListOfActiveTasks( action );
      }
      // Now work out else we need from here to calculate the later action
      getAdditionalTasksRequired( action, active_tasks );
    } else {
      std::vector<int> taskFlags( ntasks, -1 );

      #pragma omp parallel num_threads(nt)
      {
        #pragma omp for nowait
        for(unsigned i=rank; i<ntasks; i+=stride ) {
          taskIsActive( i, taskFlags[i] );
        }
      }
      for(unsigned i=0; i<ntasks; ++i) taskFlags[i] = std::abs( taskFlags[i] );
      if( !serial ) comm.Sum( taskFlags );

      unsigned nt=0;
      for(unsigned i=0; i<ntasks; ++i) {
        if( taskFlags[i]>=stride ) nt++;
      }
      active_tasks.resize(nt); nt=0;
      for(unsigned i=0; i<ntasks; ++i) {
        if( taskFlags[i]>=stride ) { active_tasks[nt]=i; nt++; }
      }
      getAdditionalTasksRequired( this, active_tasks );
    }
  } else {
    active_tasks.resize( ntasks );
    for(unsigned i=0; i<ntasks; ++i) active_tasks[i]=i;
  }
  return active_tasks;
}

bool ActionWithVector::doNotCalculateDerivatives() const {
  if( forwardPass ) return true;
  return ActionWithValue::doNotCalculateDerivatives();
}

void ActionWithVector::clearMatrixBookeeping() {
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    Value* myval = getPntrToComponent(i); if( !myval->storedata ) continue;
    if( myval->getRank()==2 && myval->getNumberOfColumns()<myval->getShape()[1] ) {
      std::fill(myval->matrix_bookeeping.begin(), myval->matrix_bookeeping.end(), 0);
    }
    myval->set(0);
  }
}

void ActionWithVector::runAllTasks() {
// Skip this if this is done elsewhere
  if( action_to_do_before ) return;

  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial) { stride=1; rank=0; }

  // Clear matrix bookeeping arrays
  ActionWithMatrix* am=dynamic_cast<ActionWithMatrix*>(this);
  if( am && stride>1 ) clearMatrixBookeeping();

  // Get the list of active tasks
  std::vector<unsigned> & partialTaskList( getListOfActiveTasks( this ) );
  unsigned nactive_tasks=partialTaskList.size();

  // Get number of threads for OpenMP
  unsigned nt=OpenMP::getNumThreads();
  if( nt*stride*10>nactive_tasks ) nt=nactive_tasks/stride/10;
  if( nt==0 ) nt=1;

  // Now do all preparations required to run all the tasks
  // prepareForTaskLoop();

  if( !action_to_do_after ) {
    forwardPass=true;
    for(unsigned i=0; i<getNumberOfComponents(); ++i) {
      if( getConstPntrToComponent(i)->getRank()==0 ) { forwardPass=false; break; }
    }
  }
  // Get the total number of streamed quantities that we need
  unsigned nquants=0, nmatrices=0, maxcol=0;
  setupStreamedComponents( getLabel(), nquants, nmatrices, maxcol );
  // Get size for buffer
  unsigned bufsize=0; getSizeOfBuffer( nactive_tasks, bufsize );
  if( buffer.size()!=bufsize ) buffer.resize( bufsize );
  // Clear buffer
  buffer.assign( buffer.size(), 0.0 );

  // Recover the number of derivatives we require
  unsigned nderivatives = 0; bool gridsInStream=checkForGrids(nderivatives);
  if( !doNotCalculateDerivatives() && !gridsInStream ) getNumberOfStreamedDerivatives( nderivatives, NULL );

  #pragma omp parallel num_threads(nt)
  {
    std::vector<double> omp_buffer;
    if( nt>1 ) omp_buffer.resize( bufsize, 0.0 );
    MultiValue myvals( nquants, nderivatives, nmatrices, maxcol );
    myvals.clearAll();

    #pragma omp for nowait
    for(unsigned i=rank; i<nactive_tasks; i+=stride) {
      // Calculate the stuff in the loop for this action
      runTask( partialTaskList[i], myvals );

      // Now transfer the data to the actions that accumulate values from the calculated quantities
      if( nt>1 ) gatherAccumulators( partialTaskList[i], myvals, omp_buffer );
      else gatherAccumulators( partialTaskList[i], myvals, buffer );

      // Clear the value
      myvals.clearAll();
    }
    #pragma omp critical
    if( nt>1 ) for(unsigned i=0; i<bufsize; ++i) buffer[i]+=omp_buffer[i];
  }

  // MPI Gather everything
  if( !serial ) {
    if( buffer.size()>0 ) comm.Sum( buffer );
    gatherProcesses();
  }
  finishComputations( buffer ); forwardPass=false;
}

void ActionWithVector::gatherProcesses() {
  ActionWithMatrix* am=dynamic_cast<ActionWithMatrix*>(this);
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    Value* myval = getPntrToComponent(i);
    if( myval->storedata && !myval->hasDeriv ) {
      comm.Sum( myval->data ); if( am && myval->getRank()==2 && myval->getNumberOfColumns()<myval->getShape()[1] ) comm.Sum( myval->matrix_bookeeping );
    }
  }
  if( action_to_do_after ) action_to_do_after->gatherProcesses();
}

bool ActionWithVector::checkForGrids( unsigned& nder ) const {
  for(int i=0; i<getNumberOfComponents(); ++i) {
    if( getConstPntrToComponent(i)->getRank()>0 && getConstPntrToComponent(i)->hasDerivatives() ) {
      nder=getConstPntrToComponent(i)->getNumberOfGridDerivatives(); return true;
    }
  }
  if( action_to_do_after ) return action_to_do_after->checkForGrids(nder);
  return false;
}

void ActionWithVector::getNumberOfTasks( unsigned& ntasks ) {
  if( ntasks==0 ) {
    if( getNumberOfArguments()==1 && getNumberOfComponents()==1 && getPntrToComponent(0)->getRank()==0 ) {
      if( !getPntrToArgument(0)->hasDerivatives() && getPntrToArgument(0)->getRank()==2 ) ntasks = getPntrToArgument(0)->getShape()[0];
      else ntasks = getPntrToArgument(0)->getNumberOfValues();
    } else {
      plumed_assert( getNumberOfComponents()>0 && getPntrToComponent(0)->getRank()>0 );
      if( getPntrToComponent(0)->hasDerivatives() ) ntasks = getPntrToComponent(0)->getNumberOfValues();
      else ntasks = getPntrToComponent(0)->getShape()[0];
    }
  }
  for(int i=0; i<getNumberOfComponents(); ++i) {
    if( getPntrToComponent(i)->getRank()==0 ) {
      if( getNumberOfArguments()!=1 ) error("mismatched numbers of tasks in streamed quantities");
      if( getPntrToArgument(0)->hasDerivatives() && ntasks!=getPntrToArgument(0)->getNumberOfValues() ) error("mismatched numbers of tasks in streamed quantities");
      else if ( !getPntrToArgument(0)->hasDerivatives() && ntasks!=getPntrToArgument(0)->getShape()[0] ) error("mismatched numbers of tasks in streamed quantities");
    } else if( getPntrToComponent(i)->hasDerivatives() && ntasks!=getPntrToComponent(i)->getNumberOfValues() ) error("mismatched numbers of tasks in streamed quantities");
    else if( !getPntrToComponent(i)->hasDerivatives() && ntasks!=getPntrToComponent(i)->getShape()[0] ) error("mismatched numbers of tasks in streamed quantities");
  }
  if( action_to_do_after ) action_to_do_after->getNumberOfTasks( ntasks );
}

void ActionWithVector::setupStreamedComponents( const std::string& headstr, unsigned& nquants, unsigned& nmat, unsigned& maxcol ) {
  for(int i=0; i<getNumberOfComponents(); ++i) { nquants++; }
}

void ActionWithVector::getSizeOfBuffer( const unsigned& nactive_tasks, unsigned& bufsize ) {
  for(int i=0; i<getNumberOfComponents(); ++i) {
    getPntrToComponent(i)->bufstart=bufsize;
    if( getPntrToComponent(i)->hasDerivatives() || getPntrToComponent(i)->getRank()==0 ) bufsize += getPntrToComponent(i)->data.size();
  }
  if( action_to_do_after ) action_to_do_after->getSizeOfBuffer( nactive_tasks, bufsize );
}

void ActionWithVector::getNumberOfStreamedDerivatives( unsigned& nderivatives, Value* stopat ) {
  std::string c=getFirstActionInChain()->getLabel();
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( !getPntrToArgument(i)->ignoreStoredValue(c) ) {
      if( getPntrToArgument(i)==stopat ) return;
      nderivatives += getPntrToArgument(i)->getNumberOfValues();
    }
  }
  if( getNumberOfAtoms()>0 ) nderivatives += 3*getNumberOfAtoms() + 9;
  // Don't do the whole chain if we have been told to stop early
  if( stopat && stopat->getPntrToAction()==this ) return;

  if( action_to_do_after ) action_to_do_after->getNumberOfStreamedDerivatives( nderivatives, stopat );
}

bool ActionWithVector::getNumberOfStoredValues( Value* startat, unsigned& nvals, const unsigned& astart, const std::vector<Value*>& stopat ) {
  for(unsigned j=astart; j<stopat.size(); ++j) {
    if( stopat[j] && (stopat[j]->getPntrToAction()==this || (stopat[j]->getPntrToAction())->checkForDependency(this)) ) return true;
  }

  std::string c=getFirstActionInChain()->getLabel();
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( !getPntrToArgument(i)->ignoreStoredValue(c) ) {
      for(unsigned j=astart; j<stopat.size(); ++j) {
        if( getPntrToArgument(i)==stopat[j] ) return true;
      }
      nvals += getPntrToArgument(i)->getNumberOfValues();
    }
  }
  if( startat->getPntrToAction()!=this && getNumberOfAtoms()>0 ) return false;

  if( action_to_do_after ) return action_to_do_after->getNumberOfStoredValues( startat, nvals, astart, stopat );
  return false;
}

void ActionWithVector::runTask( const unsigned& current, MultiValue& myvals ) const {
  const ActionWithMatrix* am = dynamic_cast<const ActionWithMatrix*>(this);
  myvals.setTaskIndex(current); myvals.vector_call=true; performTask( current, myvals );
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    const Value* myval = getConstPntrToComponent(i);
    if( am || myval->hasDerivatives() || !myval->valueIsStored() ) continue;
    Value* myv = const_cast<Value*>( myval );
    if( getName()=="RMSD_VECTOR" && myv->getRank()==2 ) continue;
    myv->set( current, myvals.get( i ) );
  }
  if( action_to_do_after ) action_to_do_after->runTask( current, myvals );
}

void ActionWithVector::gatherAccumulators( const unsigned& taskCode, const MultiValue& myvals, std::vector<double>& buffer ) const {
  for(int i=0; i<getNumberOfComponents(); ++i) {
    unsigned bufstart = getConstPntrToComponent(i)->bufstart;
    // This looks after storing of scalars that are summed from vectors/matrices
    if( getConstPntrToComponent(i)->getRank()==0 ) {
      plumed_dbg_massert( bufstart<buffer.size(), "problem in " + getLabel() );
      buffer[bufstart] += myvals.get(i);
      if( getConstPntrToComponent(i)->hasDerivatives() ) {
        for(unsigned k=0; k<myvals.getNumberActive(i); ++k) {
          unsigned kindex = myvals.getActiveIndex(i,k);
          plumed_dbg_massert( bufstart+1+kindex<buffer.size(), "problem in " + getLabel()  );
          buffer[bufstart + 1 + kindex] += myvals.getDerivative(i,kindex);
        }
      }
      // This looks after storing of vectors
    } else if( getConstPntrToComponent(i)->storedata ) gatherStoredValue( i, taskCode, myvals, bufstart, buffer );
  }
  if( action_to_do_after ) action_to_do_after->gatherAccumulators( taskCode, myvals, buffer );
}

void ActionWithVector::finishComputations( const std::vector<double>& buf ) {
  for(int i=0; i<getNumberOfComponents(); ++i) {
    // This gathers vectors and grids at the end of the calculation
    unsigned bufstart = getPntrToComponent(i)->bufstart;
    if( (getPntrToComponent(i)->getRank()>0 && getPntrToComponent(i)->hasDerivatives()) ) {
      unsigned sz_v = getPntrToComponent(i)->data.size(); getPntrToComponent(i)->data.assign( getPntrToComponent(i)->data.size(), 0 );
      for(unsigned j=0; j<sz_v; ++j) {
        plumed_dbg_assert( bufstart+j<buf.size() );
        getPntrToComponent(i)->add( j, buf[bufstart+j] );
      }
      // Make sure single values are set
    } else if( getPntrToComponent(i)->getRank()==0 ) getPntrToComponent(i)->set( buf[bufstart] );
    // This gathers derivatives of scalars
    if( !doNotCalculateDerivatives() && getPntrToComponent(i)->hasDeriv && getPntrToComponent(i)->getRank()==0 ) {
      for(unsigned j=0; j<getPntrToComponent(i)->getNumberOfDerivatives(); ++j) getPntrToComponent(i)->setDerivative( j, buf[bufstart+1+j] );
    }
  }
  if( action_to_do_after ) action_to_do_after->finishComputations( buf );
}

bool ActionWithVector::checkChainForNonScalarForces() const {
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    if( getConstPntrToComponent(i)->getRank()>0 && getConstPntrToComponent(i)->forcesWereAdded() ) return true;
  }
  if( action_to_do_after ) return action_to_do_after->checkChainForNonScalarForces();
  return false;
}

bool ActionWithVector::checkForForces() {
  if( getPntrToComponent(0)->getRank()==0 ) return ActionWithValue::checkForForces();
  else if( actionInChain() ) return false;

  // Check if there are any forces
  if( !checkChainForNonScalarForces() ) return false;

  // Setup MPI parallel loop
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial) { stride=1; rank=0; }

  // Get the number of tasks
  std::vector<unsigned> force_tasks;
  if( getenvChainForbidden()==Option::yes ) {
    std::vector<unsigned> & partialTaskList( getListOfActiveTasks( this ) );
    for(unsigned i=0; i<partialTaskList.size(); ++i) force_tasks.push_back( partialTaskList[i] );
  } else {
    getForceTasks( force_tasks ); Tools::removeDuplicates(force_tasks);
  }
  unsigned nf_tasks=force_tasks.size();
  if( nf_tasks==0 ) return false;

  // Get number of threads for OpenMP
  unsigned nt=OpenMP::getNumThreads();
  if( nt*stride*10>nf_tasks ) nt=nf_tasks/stride/10;
  if( nt==0 ) nt=1;

  // Now determine how big the multivalue needs to be
  unsigned nquants=0, nmatrices=0, maxcol=0;
  setupStreamedComponents( getLabel(), nquants, nmatrices, maxcol );
  // Recover the number of derivatives we require (this should be equal to the number of forces)
  unsigned nderiv=0; getNumberOfStreamedDerivatives( nderiv, NULL );
  if( !action_to_do_after && arg_deriv_starts.size()>0 ) {
    nderiv = 0;
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      arg_deriv_starts[i] = nderiv; nderiv += getPntrToArgument(i)->getNumberOfStoredValues();
    }
    ActionAtomistic* aa=castToActionAtomistic();
    if(aa) nderiv += 3*aa->getNumberOfAtoms() + 9;
  }
  if( forcesForApply.size()!=nderiv ) forcesForApply.resize( nderiv );
  // Clear force buffer
  forcesForApply.assign( forcesForApply.size(), 0.0 );

  #pragma omp parallel num_threads(nt)
  {
    std::vector<double> omp_forces;
    if( nt>1 ) omp_forces.resize( forcesForApply.size(), 0.0 );
    MultiValue myvals( nquants, nderiv, nmatrices, maxcol );
    myvals.clearAll();

    #pragma omp for nowait
    for(unsigned i=rank; i<nf_tasks; i+=stride) {
      runTask( force_tasks[i], myvals );

      // Now get the forces
      if( nt>1 ) gatherForces( force_tasks[i], myvals, omp_forces );
      else gatherForces( force_tasks[i], myvals, forcesForApply );

      myvals.clearAll();
    }
    #pragma omp critical
    if(nt>1) for(unsigned i=0; i<forcesForApply.size(); ++i) forcesForApply[i]+=omp_forces[i];
  }
  // MPI Gather on forces
  if( !serial ) comm.Sum( forcesForApply );
  return true;
}

bool ActionWithVector::checkComponentsForForce() const {
  for(unsigned i=0; i<values.size(); ++i) {
    if( getConstPntrToComponent(i)->forcesWereAdded() ) return true;
  }
  return false;
}

bool ActionWithVector::checkForTaskForce( const unsigned& itask, const Value* myval ) const {
  plumed_dbg_assert( !(myval->getRank()==2 && !myval->hasDerivatives()) );
  return fabs(myval->getForce(itask))>epsilon;
}

void ActionWithVector::updateForceTasksFromValue( const Value* myval, std::vector<unsigned>& force_tasks ) const {
  if( myval->getRank()>0 && myval->forcesWereAdded() ) {
    unsigned nt = myval->getNumberOfValues();
    if( !myval->hasDerivatives() ) nt = myval->getShape()[0];
    for(unsigned i=0; i<nt; ++i) {
      if( checkForTaskForce(i, myval) ) force_tasks.push_back( i );
    }
  }
}

void ActionWithVector::getForceTasks( std::vector<unsigned>& force_tasks ) const {
  if( checkComponentsForForce() ) {
    for(unsigned k=0; k<values.size(); ++k) {
      updateForceTasksFromValue( getConstPntrToComponent(k), force_tasks );
    }
  }
  if( action_to_do_after ) action_to_do_after->getForceTasks( force_tasks );
}

void ActionWithVector::gatherForcesOnStoredValue( const unsigned& ival, const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const {
  const Value* myval = getConstPntrToComponent(ival); plumed_dbg_assert( myval->storedata );
  double fforce = myval->getForce(itask);
  for(unsigned j=0; j<myvals.getNumberActive(ival); ++j) {
    unsigned jder=myvals.getActiveIndex(ival, j); plumed_dbg_assert( jder<forces.size() );
    forces[jder] += fforce*myvals.getDerivative( ival, jder );
  }
}

void ActionWithVector::gatherForces( const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const {
  if( checkComponentsForForce() ) {
    for(unsigned k=0; k<getNumberOfComponents(); ++k) {
      const Value* myval=getConstPntrToComponent(k);
      if( myval->getRank()>0 && myval->forcesWereAdded() ) gatherForcesOnStoredValue( k, itask, myvals, forces );
    }
  }
  if( action_to_do_after ) action_to_do_after->gatherForces( itask, myvals, forces );
}

void ActionWithVector::apply() {
  if( !checkForForces() ) return;
  // Find the top of the chain and add forces
  unsigned ind=0; addForcesOnArguments( 0, forcesForApply, ind ); setForcesOnAtoms( forcesForApply, ind );
}

}
