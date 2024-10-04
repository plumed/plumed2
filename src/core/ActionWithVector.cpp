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
#include "PlumedMain.h"
#include "ActionSet.h"
#include "tools/OpenMP.h"
#include "tools/Communicator.h"

namespace PLMD {

void ActionWithVector::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES");
  ActionWithArguments::registerKeywords( keys );
  keys.addFlag("SERIAL",false,"do the calculation in serial.  Do not parallelize");
}

ActionWithVector::ActionWithVector(const ActionOptions&ao):
  Action(ao),
  ActionAtomistic(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  serial(false),
  action_to_do_before(NULL),
  action_to_do_after(NULL),
  never_reduce_tasks(false),
  reduce_tasks(false),
  atomsWereRetrieved(false),
  done_in_chain(false) {
  if( keywords.exists("SERIAL") ) {
    parseFlag("SERIAL",serial);
  }
}

ActionWithVector::~ActionWithVector() {
  if( action_to_do_before ) {
    action_to_do_before->action_to_do_after=NULL;
  }
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
  if( !force && actionInChain() ) {
    return;
  }
  ActionWithValue::clearDerivatives();
  if( action_to_do_after ) {
    action_to_do_after->clearDerivatives( true );
  }
}

void ActionWithVector::clearInputForces( const bool& force ) {
  if( !force && actionInChain() ) {
    return;
  }
  ActionWithValue::clearInputForces();
  if( action_to_do_after ) {
    action_to_do_after->clearInputForces( true );
  }
}

const ActionWithVector* ActionWithVector::getFirstActionInChain() const {
  if( !actionInChain() ) {
    return this;
  }
  return action_to_do_before->getFirstActionInChain();
}

ActionWithVector* ActionWithVector::getFirstActionInChain() {
  if( !actionInChain() ) {
    return this;
  }
  return action_to_do_before->getFirstActionInChain();
}

void ActionWithVector::retrieveAtoms( const bool& force ) {
  if( !force && actionInChain() || atomsWereRetrieved ) {
    return;
  }
  ActionAtomistic::retrieveAtoms();
  atomsWereRetrieved = !actionInChain();
  if( action_to_do_after ) {
    action_to_do_after->retrieveAtoms( true );
  }
}

bool ActionWithVector::hasStoredArguments() const {
  std::string headstr=getFirstActionInChain()->getLabel();
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( !getPntrToArgument(i)->ignoreStoredValue(headstr) ) {
      return true;
    }
  }
  return false;
}

bool ActionWithVector::argumentDependsOn( const std::string& headstr, ActionWithVector* faction, Value* thearg ) {
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( this!=faction && thearg==getPntrToArgument(i) ) {
      return true;
    }
    ActionWithVector* av = dynamic_cast<ActionWithVector*>( getPntrToArgument(i)->getPntrToAction() );
    if( av && (av->getFirstActionInChain())->getLabel()==headstr ) {
      if( av->argumentDependsOn( headstr, faction, thearg ) ) {
        return true;
      };
    }
  }
  return false;
}

unsigned ActionWithVector::buildArgumentStore( const unsigned& argstart ) {
  // Don't use chains for grids
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->isConstant() ) {
      continue;
    }
    ActionWithVector* av=dynamic_cast<ActionWithVector*>(getPntrToArgument(i)->getPntrToAction());
    if( !av || getPntrToArgument(i)->getRank()>0 && getPntrToArgument(i)->hasDerivatives() ) {
      done_in_chain=false;
      break;
    }
  }
  if( done_in_chain ) {
    std::vector<std::string> alabels;
    std::vector<ActionWithVector*> f_actions;
    for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
      bool found=false;
      std::string mylab = (getPntrToArgument(i)->getPntrToAction())->getLabel();
      for(unsigned j=0; j<alabels.size(); ++j) {
        if( alabels[j]==mylab ) {
          found=true;
          break;
        }
      }
      if( !found ) {
        alabels.push_back( mylab );
      }

      // If this is calculated in setup we never need to add to chain
      if( getPntrToArgument(i)->isConstant() ) {
        continue;
      }
      // Find the chain we need to add this to from the arguments
      ActionWithVector* av=dynamic_cast<ActionWithVector*>(getPntrToArgument(i)->getPntrToAction());
      plumed_assert( av );
      found=false;
      ActionWithVector* myact = av->getFirstActionInChain();
      if( getPntrToArgument(i)->getRank()==1 && getPntrToArgument(i)->storedata ) {
        for(unsigned j=argstart; j<getNumberOfArguments(); ++j) {
          ActionWithArguments* aarg=dynamic_cast<ActionWithArguments*>( getPntrToArgument(j)->getPntrToAction() );
          if( !aarg || i==j ) {
            continue;
          }
          for(unsigned k=0; k<aarg->getNumberOfArguments(); ++k) {
            if( aarg->getPntrToArgument(k)==getPntrToArgument(i) ) {
              done_in_chain=false;
              return reallyBuildArgumentStore( argstart );
            }
          }
        }
      }
      for(unsigned j=0; j<f_actions.size(); ++j) {
        if( f_actions[j]==myact ) {
          found=true;
          break;
        }
      }
      if( !found ) {
        if( f_actions.size()>0 ) {
          if( f_actions[0]->checkForDependency(myact) ) {
            getPntrToArgument(i)->buildDataStore();
          }
          if( myact->checkForDependency(f_actions[0]) ) {
            error("cannot deal with arguments in this order. Try swapping argument order");
          }
        }
        if( !getPntrToArgument(i)->storedata && getPntrToArgument(i)->getRank()>0 ) {
          f_actions.push_back( myact );
        }
      }
    }
    // Now make sure that everything we need is in the chain
    if( f_actions.size()>0 ) {
      // Check everything for later f_actions is done before f_actions[0]
      for(unsigned i=1; i<f_actions.size(); ++i) {
        ActionWithArguments* aarg = dynamic_cast<ActionWithArguments*>( f_actions[i] );
        if( !aarg || aarg->getNumberOfArguments()==0 ) {
          continue;
        }
        for(unsigned j=0; j<aarg->getNumberOfArguments(); ++j) {
          if( (aarg->getPntrToArgument(j))->isConstant() ) {
            continue ;
          }
          bool found=false;
          std::string dep_argname = (aarg->getPntrToArgument(j))->getPntrToAction()->getLabel();
          for(const auto & pp : plumed.getActionSet()) {
            Action* p(pp.get());
            // Check if this is the dependency
            if( p->getLabel()==dep_argname ) {
              found=true;
              break;
            }
            // Check if this is the first of the arguments that will appear in this chain
            else if( p->getLabel()==f_actions[0]->getLabel() ) {
              break;
            }
          }
          if( !found ) {
            done_in_chain=false;
            break;
          }
        }
        // Stop trying to add things in the chain if we cannot
        if( !done_in_chain ) {
          return reallyBuildArgumentStore( argstart );
        }
      }
      std::vector<std::string> empty(1);
      empty[0] = f_actions[0]->getLabel();
      for(unsigned i=1; i<f_actions.size(); ++i) {
        f_actions[0]->addActionToChain( empty, f_actions[i] );
      }
    }
    // Now add this argument to the chain
    bool added=false;
    for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
      // Add this function to jobs to do in recursive loop in previous action
      if( getPntrToArgument(i)->getRank()>0 && !getPntrToArgument(i)->isConstant() ) {
        ActionWithVector* av=dynamic_cast<ActionWithVector*>( getPntrToArgument(i)->getPntrToAction() );
        if( av && av->addActionToChain( alabels, this ) ) {
          added=true;
          break;
        }
      }
    }
    plumed_massert(added, "could not add action " + getLabel() + " to chain of any of its arguments");
    // And get the number of derivatives
    ActionWithVector* head=getFirstActionInChain();
    unsigned nder=0;
    arg_deriv_starts.resize( getNumberOfArguments() );
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      // Ensures that we ignore the grid (first argument) if we are evaluating a function on a grid
      if( i==0 && getName().find("EVALUATE_FUNCTION_FROM_GRID")!=std::string::npos ) {
        continue ;
      }
      ActionWithVector* iaction=dynamic_cast<ActionWithVector*>(getPntrToArgument(i)->getPntrToAction());
      if( actionInChain() && !getPntrToArgument(i)->ignoreStoredValue(head->getLabel()) ) {
        arg_deriv_starts[i] = 0;
        head->getNumberOfStreamedDerivatives( arg_deriv_starts[i], getPntrToArgument(i) );
      } else if( iaction && iaction->isInSubChain(nder) ) {
        arg_deriv_starts[i] = nder;
        // Add the total number of derivatives that we have by this point in the chain to nder
        if( iaction ) {
          nder=0;
          head->getNumberOfStreamedDerivatives( nder, getPntrToArgument(i) );
        }
      } else {
        // Check if we have already found this action
        int k=-1;
        if( iaction ) {
          ActionWithVector* ider_action=iaction->getActionWithDerivatives( iaction );
          for(unsigned j=0; j<i; ++j) {
            if( j==0 && getName().find("EVALUATE_FUNCTION_FROM_GRID")!=std::string::npos ) {
              continue ;
            }
            ActionWithVector* jaction=dynamic_cast<ActionWithVector*>(getPntrToArgument(j)->getPntrToAction());
            if( jaction->getActionWithDerivatives(jaction)==ider_action || jaction->checkForDependency(ider_action) ) {
              k=j;
              break;
            }
          }
          if( k>=0 ) {
            arg_deriv_starts[i] = arg_deriv_starts[k];
            continue;
          }
        }

        if( i>0 ) {
          // This is a fudge so that inputs like this work:
          // c: CONTACT_MATRIX ATOMS=1-100
          // d: MATRIX_PRODUCT ARG=mat1,mat2
          // e: CUSTOM ARG=c,d
          // f: MATRIX_PRODUCT ARG=mat3,mat4
          // g: CUSTOM ARG=c,f
          // See symfunc rt-nbonds-q6 for an example
          // In this example when we set arg_deriv_starts[1] for f in g nder=number of derivatives of c
          // mder is equal to the number of derivatives by the time you get to f minus the number of derivatives for c
          unsigned mder=0;
          ActionWithVector* jaction=dynamic_cast<ActionWithVector*>(getPntrToArgument(i-1)->getPntrToAction());
          if( jaction->action_to_do_after && !(jaction->action_to_do_after)->getNumberOfStoredValues( getPntrToArgument(i-1), mder, i, getArguments() ) ) {
            mder=0;
          }
          if( mder>0 ) {
            nder = nder + mder;
          }
        }

        arg_deriv_starts[i] = nder;
        // Add the total number of derivatives that we have by this point in the chain to nder
        if( iaction ) {
          nder=0;
          if( (getPntrToArgument(i)->getPntrToAction())->getName().find("DIFFERENCE")!=std::string::npos ) {
            ActionWithArguments* aarg=dynamic_cast<ActionWithArguments*>( getPntrToArgument(i)->getPntrToAction() );
            plumed_assert( aarg && aarg->getNumberOfArguments()==2 );
            head->getNumberOfStreamedDerivatives( nder, aarg->getPntrToArgument(0) );
            nder += (aarg->getPntrToArgument(1))->getNumberOfValues();
          } else {
            head->getNumberOfStreamedDerivatives( nder, getPntrToArgument(i) );
          }
        }
      }
    }
    nder=0;
    head->getNumberOfStreamedDerivatives( nder, NULL );
    return nder;
  }
  return reallyBuildArgumentStore( argstart );
}

unsigned ActionWithVector::reallyBuildArgumentStore( const unsigned& argstart ) {
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    if( getPntrToArgument(i)->getRank()>0 ) {
      getPntrToArgument(i)->buildDataStore();
    }
  }
  unsigned nder=0;
  arg_deriv_starts.resize( getNumberOfArguments() );
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    arg_deriv_starts[i] = nder;
    nder += getPntrToArgument(i)->getNumberOfValues();
  }
  return nder;
}

ActionWithVector* ActionWithVector::getActionWithDerivatives( ActionWithVector* depaction ) {
  if( depaction==this || depaction->checkForDependency(this) ) {
    if( getNumberOfAtoms()>0 ) {
      return this;
    }
    std::string c=getFirstActionInChain()->getLabel();
    for(unsigned i=0; i<getNumberOfArguments(); ++i) {
      if( !getPntrToArgument(i)->ignoreStoredValue(c) && !getPntrToArgument(i)->isConstant() ) {
        return this;
      }
    }
  }
  plumed_assert( action_to_do_before );
  return action_to_do_before->getActionWithDerivatives(depaction);
}

bool ActionWithVector::addActionToChain( const std::vector<std::string>& alabels, ActionWithVector* act ) {
  if( action_to_do_after ) {
    bool state=action_to_do_after->addActionToChain( alabels, act );
    return state;
  }
  // Check action is not already in chain
  std::vector<std::string> mylabels;
  getFirstActionInChain()->getAllActionLabelsInChain( mylabels );
  for(unsigned i=0; i<mylabels.size(); ++i) {
    if( act->getLabel()==mylabels[i] ) {
      return true;
    }
  }

  // Check that everything that is required has been calculated
  for(unsigned i=0; i<alabels.size(); ++i) {
    bool found=false;
    for(unsigned j=0; j<mylabels.size(); ++j) {
      if( alabels[i]==mylabels[j] ) {
        found=true;
        break;
      }
    }
    if( !found ) {
      ActionWithValue* av=plumed.getActionSet().selectWithLabel<ActionWithValue*>( alabels[i] );
      plumed_massert( av, "could not cast " + alabels[i] );
      bool storingall=true;
      for(int j=0; j<av->getNumberOfComponents(); ++j) {
        if( !(av->getPntrToComponent(j))->storedata ) {
          storingall=false;
        }
      }
      if( !storingall ) {
        return false;
      }
    }
  }
  // This checks that there is nothing that will cause problems in the chain
  mylabels.resize(0);
  getFirstActionInChain()->getAllActionLabelsInChain( mylabels );
  for(unsigned i=0; i<mylabels.size(); ++i) {
    ActionWithVector* av1=plumed.getActionSet().selectWithLabel<ActionWithVector*>( mylabels[i] );
    for(unsigned j=0; j<i; ++j) {
      ActionWithVector* av2=plumed.getActionSet().selectWithLabel<ActionWithVector*>( mylabels[j] );
      if( !av1->canBeAfterInChain( av2 ) ) {
        error("must calculate " + mylabels[j] + " before " + mylabels[i] );
      }
    }
  }
  action_to_do_after=act;
  act->action_to_do_before=this;
  updateTaskListReductionStatus();
  ActionWithVector* head = getFirstActionInChain();
  head->broadcastThatTasksAreReduced( head );
  head->finishChainBuild( act );
  return true;
}

void ActionWithVector::updateTaskListReductionStatus() {
  ActionWithVector* head = getFirstActionInChain();
  std::vector<ActionWithVector*> task_reducing_actions;
  head->canReduceTasks( task_reducing_actions );
  if( task_reducing_actions.size()>0 ) {
    head->reduce_tasks=true;
  }
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
          if( aselect==av_head->task_control_list[i] ) {
            found=true;
            break;
          }
        }
        if( !found ) {
          av_head->task_control_list.insert( av_head->task_control_list.begin(), aselect );
        }

        av_head->reduce_tasks=true;
        av_head->updateTaskReductionFlag( av_head->reduce_tasks );
      }
    }
  }
  if( action_to_do_after ) {
    action_to_do_after->broadcastThatTasksAreReduced( aselect );
  }
}

void ActionWithVector::updateTaskReductionFlag( bool& head_reduce_tasks ) {
  if( actionInChain() ) {
    plumed_assert( task_control_list.size()==0 );
  } else {
    for(unsigned i=0; i<task_control_list.size(); ++i) {
      if( !(task_control_list[i]->getFirstActionInChain())->reduce_tasks ) {
        head_reduce_tasks=false;
      }
    }
  }
  broadcastThatTasksAreReduced( getFirstActionInChain() );
  if( action_to_do_after ) {
    action_to_do_after->updateTaskReductionFlag( head_reduce_tasks );
  }
}

void ActionWithVector::canReduceTasks( std::vector<ActionWithVector*>& task_reducing_actions ) {
  areAllTasksRequired( task_reducing_actions );
  if( action_to_do_after ) {
    action_to_do_after->canReduceTasks( task_reducing_actions );
  }
}

void ActionWithVector::finishChainBuild( ActionWithVector* act ) {
  if( action_to_do_after ) {
    action_to_do_after->finishChainBuild( act );
  }
}

void ActionWithVector::getAllActionLabelsInChain( std::vector<std::string>& mylabels ) const {
  bool found = false ;
  for(unsigned i=0; i<mylabels.size(); ++i) {
    if( getLabel()==mylabels[i] ) {
      found=true;
    }
  }
  if( !found ) {
    mylabels.push_back( getLabel() );
  }
  if( action_to_do_after ) {
    action_to_do_after->getAllActionLabelsInChain( mylabels );
  }
}

void ActionWithVector::taskIsActive( const unsigned& current, int& flag ) const {
  if( isActive() ) {
    flag = checkTaskStatus( current, flag );
  }
  if( flag<=0 && action_to_do_after ) {
    action_to_do_after->taskIsActive( current, flag );
  }
}

void ActionWithVector::getAdditionalTasksRequired( ActionWithVector* action, std::vector<unsigned>& atasks ) {
  for(unsigned i=0; i<task_control_list.size(); ++i ) {
    task_control_list[i]->getAdditionalTasksRequired( action, atasks );
  }
}

void ActionWithVector::prepare() {
  active_tasks.resize(0);
  atomsWereRetrieved=false;
}

std::vector<unsigned>& ActionWithVector::getListOfActiveTasks( ActionWithVector* action ) {
  if( active_tasks.size()>0 ) {
    return active_tasks;
  }
  unsigned ntasks=0;
  getNumberOfTasks( ntasks );

  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial) {
    stride=1;
    rank=0;
  }

  // Get number of threads for OpenMP
  unsigned nt=OpenMP::getNumThreads();
  if( nt*stride*10>ntasks ) {
    nt=ntasks/stride/10;
  }
  if( nt==0 ) {
    nt=1;
  }

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
      for(unsigned i=0; i<ntasks; ++i) {
        taskFlags[i] = std::abs( taskFlags[i] );
      }
      if( !serial ) {
        comm.Sum( taskFlags );
      }

      unsigned nt=0;
      for(unsigned i=0; i<ntasks; ++i) {
        if( taskFlags[i]>=stride ) {
          nt++;
        }
      }
      active_tasks.resize(nt);
      nt=0;
      for(unsigned i=0; i<ntasks; ++i) {
        if( taskFlags[i]>=stride ) {
          active_tasks[nt]=i;
          nt++;
        }
      }
      getAdditionalTasksRequired( this, active_tasks );
    }
  } else {
    active_tasks.resize( ntasks );
    for(unsigned i=0; i<ntasks; ++i) {
      active_tasks[i]=i;
    }
  }
  return active_tasks;
}

void ActionWithVector::runAllTasks() {
// Skip this if this is done elsewhere
  if( action_to_do_before ) {
    return;
  }

  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial) {
    stride=1;
    rank=0;
  }

  // Get the list of active tasks
  std::vector<unsigned> & partialTaskList( getListOfActiveTasks( this ) );
  unsigned nactive_tasks=partialTaskList.size();

  // Get number of threads for OpenMP
  unsigned nt=OpenMP::getNumThreads();
  if( nt*stride*10>nactive_tasks ) {
    nt=nactive_tasks/stride/10;
  }
  if( nt==0 ) {
    nt=1;
  }

  // Now do all preparations required to run all the tasks
  // prepareForTaskLoop();

  // Get the total number of streamed quantities that we need
  unsigned nquants=0, nmatrices=0, maxcol=0, nbooks=0;
  getNumberOfStreamedQuantities( getLabel(), nquants, nmatrices, maxcol, nbooks );
  // Get size for buffer
  unsigned bufsize=0;
  getSizeOfBuffer( nactive_tasks, bufsize );
  if( buffer.size()!=bufsize ) {
    buffer.resize( bufsize );
  }
  // Clear buffer
  buffer.assign( buffer.size(), 0.0 );

  // Recover the number of derivatives we require
  unsigned nderivatives = 0;
  bool gridsInStream=checkForGrids(nderivatives);
  if( !doNotCalculateDerivatives() && !gridsInStream ) {
    getNumberOfStreamedDerivatives( nderivatives, NULL );
  }

  #pragma omp parallel num_threads(nt)
  {
    std::vector<double> omp_buffer;
    if( nt>1 ) {
      omp_buffer.resize( bufsize, 0.0 );
    }
    MultiValue myvals( nquants, nderivatives, nmatrices, maxcol, nbooks );
    myvals.clearAll();

    #pragma omp for nowait
    for(unsigned i=rank; i<nactive_tasks; i+=stride) {
      // Calculate the stuff in the loop for this action
      runTask( partialTaskList[i], myvals );

      // Now transfer the data to the actions that accumulate values from the calculated quantities
      if( nt>1 ) {
        gatherAccumulators( partialTaskList[i], myvals, omp_buffer );
      } else {
        gatherAccumulators( partialTaskList[i], myvals, buffer );
      }

      // Clear the value
      myvals.clearAll();
    }
    #pragma omp critical
    gatherThreads( nt, bufsize, omp_buffer, buffer, myvals );
  }

  // MPI Gather everything
  if( !serial && buffer.size()>0 ) {
    gatherProcesses( buffer );
  }
  finishComputations( buffer );
}

void ActionWithVector::gatherThreads( const unsigned& nt, const unsigned& bufsize, const std::vector<double>& omp_buffer, std::vector<double>& buffer, MultiValue& myvals ) {
  if( nt>1 )
    for(unsigned i=0; i<bufsize; ++i) {
      buffer[i]+=omp_buffer[i];
    }
}

void ActionWithVector::gatherProcesses( std::vector<double>& buffer ) {
  comm.Sum( buffer );
}

bool ActionWithVector::checkForGrids( unsigned& nder ) const {
  for(int i=0; i<getNumberOfComponents(); ++i) {
    if( getConstPntrToComponent(i)->getRank()>0 && getConstPntrToComponent(i)->hasDerivatives() ) {
      nder=getConstPntrToComponent(i)->getNumberOfGridDerivatives();
      return true;
    }
  }
  if( action_to_do_after ) {
    return action_to_do_after->checkForGrids(nder);
  }
  return false;
}

void ActionWithVector::getNumberOfTasks( unsigned& ntasks ) {
  if( ntasks==0 ) {
    if( getNumberOfArguments()==1 && getNumberOfComponents()==1 && getPntrToComponent(0)->getRank()==0 ) {
      if( !getPntrToArgument(0)->hasDerivatives() && getPntrToArgument(0)->getRank()==2 ) {
        ntasks = getPntrToArgument(0)->getShape()[0];
      } else {
        ntasks = getPntrToArgument(0)->getNumberOfValues();
      }
    } else {
      plumed_assert( getNumberOfComponents()>0 && getPntrToComponent(0)->getRank()>0 );
      if( getPntrToComponent(0)->hasDerivatives() ) {
        ntasks = getPntrToComponent(0)->getNumberOfValues();
      } else {
        ntasks = getPntrToComponent(0)->getShape()[0];
      }
    }
  }
  for(int i=0; i<getNumberOfComponents(); ++i) {
    if( getPntrToComponent(i)->getRank()==0 ) {
      if( getNumberOfArguments()!=1 ) {
        error("mismatched numbers of tasks in streamed quantities");
      }
      if( getPntrToArgument(0)->hasDerivatives() && ntasks!=getPntrToArgument(0)->getNumberOfValues() ) {
        error("mismatched numbers of tasks in streamed quantities");
      } else if ( !getPntrToArgument(0)->hasDerivatives() && ntasks!=getPntrToArgument(0)->getShape()[0] ) {
        error("mismatched numbers of tasks in streamed quantities");
      }
    } else if( getPntrToComponent(i)->hasDerivatives() && ntasks!=getPntrToComponent(i)->getNumberOfValues() ) {
      error("mismatched numbers of tasks in streamed quantities");
    } else if( !getPntrToComponent(i)->hasDerivatives() && ntasks!=getPntrToComponent(i)->getShape()[0] ) {
      error("mismatched numbers of tasks in streamed quantities");
    }
  }
  if( action_to_do_after ) {
    action_to_do_after->getNumberOfTasks( ntasks );
  }
}

void ActionWithVector::setupStreamedComponents( const std::string& headstr, unsigned& nquants, unsigned& nmat, unsigned& maxcol, unsigned& nbookeeping ) {
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( !getPntrToArgument(i)->ignoreStoredValue(headstr) ) {
      getPntrToArgument(i)->streampos=nquants;
      nquants++;
    }
  }
  for(int i=0; i<getNumberOfComponents(); ++i) {
    getPntrToComponent(i)->streampos=nquants;
    nquants++;
  }
}

void ActionWithVector::getNumberOfStreamedQuantities( const std::string& headstr, unsigned& nquants, unsigned& nmat, unsigned& maxcol, unsigned& nbookeeping ) {
  setupStreamedComponents( headstr, nquants, nmat, maxcol, nbookeeping );
  if( action_to_do_after ) {
    action_to_do_after->getNumberOfStreamedQuantities( headstr, nquants, nmat, maxcol, nbookeeping );
  }
}

void ActionWithVector::getSizeOfBuffer( const unsigned& nactive_tasks, unsigned& bufsize ) {
  for(int i=0; i<getNumberOfComponents(); ++i) {
    getPntrToComponent(i)->bufstart=bufsize;
    bufsize += getPntrToComponent(i)->data.size();
  }
  if( action_to_do_after ) {
    action_to_do_after->getSizeOfBuffer( nactive_tasks, bufsize );
  }
}

void ActionWithVector::getNumberOfStreamedDerivatives( unsigned& nderivatives, Value* stopat ) {
  std::string c=getFirstActionInChain()->getLabel();
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( !getPntrToArgument(i)->ignoreStoredValue(c) ) {
      if( getPntrToArgument(i)==stopat ) {
        return;
      }
      nderivatives += getPntrToArgument(i)->getNumberOfValues();
    }
  }
  if( getNumberOfAtoms()>0 ) {
    nderivatives += 3*getNumberOfAtoms() + 9;
  }
  // Don't do the whole chain if we have been told to stop early
  if( stopat && stopat->getPntrToAction()==this ) {
    return;
  }

  if( action_to_do_after ) {
    action_to_do_after->getNumberOfStreamedDerivatives( nderivatives, stopat );
  }
}

bool ActionWithVector::getNumberOfStoredValues( Value* startat, unsigned& nvals, const unsigned& astart, const std::vector<Value*>& stopat ) {
  for(unsigned j=astart; j<stopat.size(); ++j) {
    if( stopat[j] && (stopat[j]->getPntrToAction()==this || (stopat[j]->getPntrToAction())->checkForDependency(this)) ) {
      return true;
    }
  }

  std::string c=getFirstActionInChain()->getLabel();
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if( !getPntrToArgument(i)->ignoreStoredValue(c) ) {
      for(unsigned j=astart; j<stopat.size(); ++j) {
        if( getPntrToArgument(i)==stopat[j] ) {
          return true;
        }
      }
      nvals += getPntrToArgument(i)->getNumberOfValues();
    }
  }
  if( startat->getPntrToAction()!=this && getNumberOfAtoms()>0 ) {
    return false;
  }

  if( action_to_do_after ) {
    return action_to_do_after->getNumberOfStoredValues( startat, nvals, astart, stopat );
  }
  return false;
}

void ActionWithVector::runTask( const unsigned& current, MultiValue& myvals ) const {
  if( isActive() ) {
    myvals.setTaskIndex(current);
    myvals.vector_call=true;
    performTask( current, myvals );
  }
  if( action_to_do_after ) {
    action_to_do_after->runTask( current, myvals );
  }
}

void ActionWithVector::gatherAccumulators( const unsigned& taskCode, const MultiValue& myvals, std::vector<double>& buffer ) const {
  if( isActive() ) {
    for(int i=0; i<getNumberOfComponents(); ++i) {
      unsigned bufstart = getConstPntrToComponent(i)->bufstart;
      // This looks after storing of scalars that are summed from vectors/matrices
      if( getConstPntrToComponent(i)->getRank()==0 ) {
        plumed_dbg_massert( bufstart<buffer.size(), "problem in " + getLabel() );
        unsigned sind = getConstPntrToComponent(i)->streampos;
        buffer[bufstart] += myvals.get(sind);
        if( getConstPntrToComponent(i)->hasDerivatives() ) {
          for(unsigned k=0; k<myvals.getNumberActive(sind); ++k) {
            unsigned kindex = myvals.getActiveIndex(sind,k);
            plumed_dbg_massert( bufstart+1+kindex<buffer.size(), "problem in " + getLabel()  );
            buffer[bufstart + 1 + kindex] += myvals.getDerivative(sind,kindex);
          }
        }
        // This looks after storing of vectors
      } else if( getConstPntrToComponent(i)->storedata ) {
        gatherStoredValue( i, taskCode, myvals, bufstart, buffer );
      }
    }
  }
  if( action_to_do_after ) {
    action_to_do_after->gatherAccumulators( taskCode, myvals, buffer );
  }
}

void ActionWithVector::gatherStoredValue( const unsigned& valindex, const unsigned& taskCode, const MultiValue& myvals,
    const unsigned& bufstart, std::vector<double>& buffer ) const {
  plumed_dbg_assert( getConstPntrToComponent(valindex)->getRank()==1 && !getConstPntrToComponent(valindex)->hasDeriv );
  unsigned vindex = getConstPntrToComponent(valindex)->bufstart + taskCode;
  plumed_dbg_massert( vindex<buffer.size(), "failing in " + getLabel() );
  buffer[vindex] += myvals.get(getConstPntrToComponent(valindex)->streampos);
}

void ActionWithVector::finishComputations( const std::vector<double>& buf ) {
  if( isActive() ) {
    for(int i=0; i<getNumberOfComponents(); ++i) {
      // This gathers vectors and grids at the end of the calculation
      unsigned bufstart = getPntrToComponent(i)->bufstart;
      getPntrToComponent(i)->data.assign( getPntrToComponent(i)->data.size(), 0 );
      if( (getPntrToComponent(i)->getRank()>0 && getPntrToComponent(i)->hasDerivatives()) || getPntrToComponent(i)->storedata ) {
        unsigned sz_v = getPntrToComponent(i)->data.size();
        for(unsigned j=0; j<sz_v; ++j) {
          plumed_dbg_assert( bufstart+j<buf.size() );
          getPntrToComponent(i)->add( j, buf[bufstart+j] );
        }
        // Make sure single values are set
      } else if( getPntrToComponent(i)->getRank()==0 ) {
        getPntrToComponent(i)->set( buf[bufstart] );
      }
      // This gathers derivatives of scalars
      if( !doNotCalculateDerivatives() && getPntrToComponent(i)->hasDeriv && getPntrToComponent(i)->getRank()==0 ) {
        for(unsigned j=0; j<getPntrToComponent(i)->getNumberOfDerivatives(); ++j) {
          getPntrToComponent(i)->setDerivative( j, buf[bufstart+1+j] );
        }
      }
    }
  }
  if( action_to_do_after ) {
    action_to_do_after->finishComputations( buf );
  }
}

bool ActionWithVector::checkChainForNonScalarForces() const {
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    if( getConstPntrToComponent(i)->getRank()>0 && getConstPntrToComponent(i)->forcesWereAdded() ) {
      return true;
    }
  }
  if( action_to_do_after ) {
    return action_to_do_after->checkChainForNonScalarForces();
  }
  return false;
}

bool ActionWithVector::checkForForces() {
  if( getPntrToComponent(0)->getRank()==0 ) {
    return ActionWithValue::checkForForces();
  } else if( actionInChain() ) {
    return false;
  }

  // Check if there are any forces
  if( !checkChainForNonScalarForces() ) {
    return false;
  }

  // Setup MPI parallel loop
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial) {
    stride=1;
    rank=0;
  }

  // Get the number of tasks
  std::vector<unsigned> force_tasks;
  getForceTasks( force_tasks );
  Tools::removeDuplicates(force_tasks);
  unsigned nf_tasks=force_tasks.size();

  // Get number of threads for OpenMP
  unsigned nt=OpenMP::getNumThreads();
  if( nt*stride*10>nf_tasks ) {
    nt=nf_tasks/stride/10;
  }
  if( nt==0 ) {
    nt=1;
  }

  // Now determine how big the multivalue needs to be
  unsigned nquants=0, nmatrices=0, maxcol=0, nbooks=0;
  getNumberOfStreamedQuantities( getLabel(), nquants, nmatrices, maxcol, nbooks );
  // Recover the number of derivatives we require (this should be equal to the number of forces)
  unsigned nderiv=0;
  getNumberOfStreamedDerivatives( nderiv, NULL );
  if( forcesForApply.size()!=nderiv ) {
    forcesForApply.resize( nderiv );
  }
  // Clear force buffer
  forcesForApply.assign( forcesForApply.size(), 0.0 );

  #pragma omp parallel num_threads(nt)
  {
    std::vector<double> omp_forces;
    if( nt>1 ) {
      omp_forces.resize( forcesForApply.size(), 0.0 );
    }
    MultiValue myvals( nquants, nderiv, nmatrices, maxcol, nbooks );
    myvals.clearAll();

    #pragma omp for nowait
    for(unsigned i=rank; i<nf_tasks; i+=stride) {
      runTask( force_tasks[i], myvals );

      // Now get the forces
      if( nt>1 ) {
        gatherForces( force_tasks[i], myvals, omp_forces );
      } else {
        gatherForces( force_tasks[i], myvals, forcesForApply );
      }

      myvals.clearAll();
    }
    #pragma omp critical
    if(nt>1)
      for(unsigned i=0; i<forcesForApply.size(); ++i) {
        forcesForApply[i]+=omp_forces[i];
      }
  }
  // MPI Gather on forces
  if( !serial ) {
    comm.Sum( forcesForApply );
  }
  return true;
}

bool ActionWithVector::checkComponentsForForce() const {
  for(unsigned i=0; i<values.size(); ++i) {
    if( getConstPntrToComponent(i)->forcesWereAdded() ) {
      return true;
    }
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
    if( !myval->hasDerivatives() ) {
      nt = myval->getShape()[0];
    }
    for(unsigned i=0; i<nt; ++i) {
      if( checkForTaskForce(i, myval) ) {
        force_tasks.push_back( i );
      }
    }
  }
}

void ActionWithVector::getForceTasks( std::vector<unsigned>& force_tasks ) const {
  if( isActive() && checkComponentsForForce() ) {
    for(unsigned k=0; k<values.size(); ++k) {
      updateForceTasksFromValue( getConstPntrToComponent(k), force_tasks );
    }
  }
  if( action_to_do_after ) {
    action_to_do_after->getForceTasks( force_tasks );
  }
}

void ActionWithVector::gatherForcesOnStoredValue( const Value* myval, const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const {
  plumed_dbg_assert( myval->storedata && !(myval->getRank()==2 && !myval->hasDerivatives()) );
  double fforce = myval->getForce(itask);
  unsigned sspos = myval->getPositionInStream();
  for(unsigned j=0; j<myvals.getNumberActive(sspos); ++j) {
    unsigned jder=myvals.getActiveIndex(sspos, j);
    plumed_dbg_assert( jder<forces.size() );
    forces[jder] += fforce*myvals.getDerivative( sspos, jder );
  }
}

void ActionWithVector::gatherForces( const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const {
  if( isActive() && checkComponentsForForce() ) {
    for(unsigned k=0; k<getNumberOfComponents(); ++k) {
      const Value* myval=getConstPntrToComponent(k);
      if( myval->getRank()>0 && myval->forcesWereAdded() ) {
        gatherForcesOnStoredValue( myval, itask, myvals, forces );
      }
    }
  }
  if( action_to_do_after ) {
    action_to_do_after->gatherForces( itask, myvals, forces );
  }
}

void ActionWithVector::apply() {
  if( !checkForForces() ) {
    return;
  }
  // Find the top of the chain and add forces
  unsigned ind=0;
  getFirstActionInChain()->addForcesToInput( getForcesToApply(), ind );
}

void ActionWithVector::addForcesToInput( const std::vector<double>& forcesToApply, unsigned& ind ) {
  if( ind>=forcesToApply.size() ) {
    return;
  }
  addForcesOnArguments( 0, forcesToApply, ind, getFirstActionInChain()->getLabel() );
  setForcesOnAtoms( forcesToApply, ind );
  if( action_to_do_after ) {
    action_to_do_after->addForcesToInput( forcesToApply, ind );
  }
}

}
