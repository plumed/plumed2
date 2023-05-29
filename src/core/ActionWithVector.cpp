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
  ActionWithValue::registerKeywords( keys ); keys.remove("NUMERICAL_DERIVATIVES");
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
  done_in_chain(false)
{
  if( keywords.exists("SERIAL") ) parseFlag("SERIAL",serial);
}

void ActionWithVector::lockRequests() {
  ActionAtomistic::lockRequests();
  ActionWithArguments::lockRequests();
} 
  
void ActionWithVector::unlockRequests() {
  ActionAtomistic::unlockRequests();
  ActionWithArguments::unlockRequests();
}

void ActionWithVector::clearDerivatives( const bool& force ) {
  if( !force && actionInChain() ) return;
  ActionWithValue::clearDerivatives();
  if( action_to_do_after ) action_to_do_after->clearDerivatives( true );
}

unsigned ActionWithVector::buildArgumentStore( const unsigned& argstart ) {
  if( done_in_chain ) {
      std::vector<std::string> alabels; std::vector<ActionWithVector*> f_actions;
      for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
          bool found=false; std::string mylab = (getPntrToArgument(i)->getPntrToAction())->getLabel();
          for(unsigned j=0; j<alabels.size(); ++j) {
            if( alabels[j]==mylab ) { found=true; break; }
          }          
          if( !found ) alabels.push_back( mylab );
                     
          // If this is calculated in setup we never need to add to chain 
          if( getPntrToArgument(i)->isConstant() ) continue;
          // Check if 
          ActionWithVector* av=dynamic_cast<ActionWithVector*>(getPntrToArgument(i)->getPntrToAction()); plumed_assert( av );
          found=false; ActionWithVector* myact = av->getActionThatCalculates(); 
          for(unsigned j=0; j<f_actions.size(); ++j) {
            if( f_actions[j]==myact ) { found=true; break; }
          }
          if( !found ) {
              if( !getPntrToArgument(i)->storedata && getPntrToArgument(i)->getRank()>0 ) f_actions.push_back( myact );
          }
      }
      // Now make sure that everything we need is in the chain 
      if( f_actions.size()>0 ) {
          std::vector<std::string> empty(1); empty[0] = f_actions[0]->getLabel();
          for(unsigned i=1; i<f_actions.size(); ++i) f_actions[0]->addActionToChain( empty, f_actions[i] );
      } 
      // Now add this argument to the chain
      bool added=false;
      for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
        // Add this function to jobs to do in recursive loop in previous action
        if( getPntrToArgument(i)->getRank()>0 && !getPntrToArgument(i)->isConstant() ) {
            ActionWithVector* av=dynamic_cast<ActionWithVector*>( getPntrToArgument(i)->getPntrToAction() );
            if( av && av->addActionToChain( alabels, this ) ) { added=true; break; }
        }
      }
      plumed_massert(added, "could not add action " + getLabel() + " to chain of any of its arguments"); 

      // Setup the distinct arguments and the number of derivatives
      ActionWithVector* av=dynamic_cast<ActionWithVector*>( getPntrToArgument(argstart)->getPntrToAction() ); plumed_assert( av );
      if( av->getNumberOfArguments()==0 || av->mustBeTreatedAsDistinctArguments(0) ) distinct_arguments.push_back( av ); 
      else distinct_arguments.push_back( av->getFirstNonStream() );

      arg_deriv_starts.push_back(0); unsigned nder;
      if( getPntrToArgument(argstart)->getRank()==0 || getPntrToArgument(0)->storedata ) nder = getPntrToArgument(0)->getNumberOfValues();
      else nder = distinct_arguments[0]->getNumberOfDerivatives();

      // Now deal with the rest of the arguments
      for(unsigned i=argstart+1; i<getNumberOfArguments(); ++i) {
          ActionWithVector* myval; ActionWithVector* av=dynamic_cast<ActionWithVector*>( getPntrToArgument(i)->getPntrToAction() ); plumed_assert( av );
          if( av->getNumberOfArguments()==0 || av->mustBeTreatedAsDistinctArguments(0) ) myval = av;
          else myval = av->getFirstNonStream();

          // Check we haven't already dealt with this argument
          int argno=-1;
          for(unsigned j=0; j<distinct_arguments.size(); ++j) {
            if( myval==distinct_arguments[j] ) { argno=j; break; }
          }
          if( argno>=0 ) { arg_deriv_starts.push_back( arg_deriv_starts[argno] ); continue; }
          // Add if if we haven't dealt with it
          distinct_arguments.push_back( myval ); arg_deriv_starts.push_back( nder );
          if( getPntrToArgument(i)->getRank()==0 || getPntrToArgument(i)->storedata ) nder += getPntrToArgument(i)->getNumberOfValues();
          else nder += myval->getNumberOfDerivatives();
      }
      arg_deriv_starts.push_back( nder ); return nder;
  } 
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) { if( getPntrToArgument(i)->getRank()>0 ) getPntrToArgument(i)->buildDataStore(); }
  unsigned nder=0; for(unsigned i=0; i<getNumberOfArguments(); ++i) nder += getPntrToArgument(i)->getNumberOfValues();
  return nder;
}

bool ActionWithVector::mustBeTreatedAsDistinctArguments( const unsigned& argstart ) {
  if( !done_in_chain ) return true;

  if( getNumberOfArguments()-argstart==1 ) {
      ActionWithVector* cal = getActionThatCalculates();
      if( cal->getNumberOfAtoms()>0 ) return true;
      return false;
  } 
    
  std::vector<const ActionWithValue*> allvals; 
  ActionWithVector* av=dynamic_cast<ActionWithVector*>(getPntrToArgument(argstart)->getPntrToAction());
  if( av ) av->getAllActionsRequired( allvals );
  for(unsigned j=argstart+1;j<getNumberOfArguments();++j) {
      std::vector<const ActionWithValue*> tvals; av = dynamic_cast<ActionWithVector*>(getPntrToArgument(j)->getPntrToAction());
      if( av ) av->getAllActionsRequired( tvals );
      if( allvals.size()!=tvals.size() ) return true;
    
      for(unsigned k=0;k<tvals.size();++k) {
          if( tvals[k]!=allvals[k] ) return true;
      } 
  }       
  return false;
} 

void ActionWithVector::getAllActionsRequired( std::vector<const ActionWithValue*>& allvals ) const {
  if( getNumberOfArguments()==0 ) {
      bool found = false;
      for(unsigned k=0;k<allvals.size();++k) {
          if( allvals[k]==this ) { found=true; break; }
      }                          
      if( !found ) allvals.push_back( this );
  } else {                       
     for(unsigned j=0;j<getNumberOfArguments();++j) {
        ActionWithVector* av=dynamic_cast<ActionWithVector*>(getPntrToArgument(j)->getPntrToAction());
        if(av) av->getAllActionsRequired( allvals );
        bool found = false;
        for(unsigned k=0;k<allvals.size();++k) {
            if( allvals[k]==getPntrToArgument(j)->getPntrToAction() ) { found=true; break; }
        }
        if( !found ) allvals.push_back( getPntrToArgument(j)->getPntrToAction() );
     }
  }
}    

ActionWithVector* ActionWithVector::getFirstNonStream() {
  ActionWithVector* av=dynamic_cast<ActionWithVector*>( getPntrToArgument(0)->getPntrToAction() );
  if( av->getNumberOfArguments()==0 || av->mustBeTreatedAsDistinctArguments(0) ) return this;

  for(unsigned i=1;i<getNumberOfArguments();++i) {
      ActionWithVector* aav=dynamic_cast<ActionWithVector*>( getPntrToArgument(i)->getPntrToAction() ); 
      if( av!=aav ) return this; 
  } 
  return av->getFirstNonStream();
}

bool ActionWithVector::addActionToChain( const std::vector<std::string>& alabels, ActionWithVector* act ) {
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
      ActionWithVector* av=plumed.getActionSet().selectWithLabel<ActionWithVector*>( alabels[i] );
      bool storingall=true;
      for(int j=0; j<av->getNumberOfComponents(); ++j) {
        if( !(av->getPntrToComponent(j))->storedata ) storingall=false;
      }
      if( !storingall ) return false;
    }                            
  }
  // This checks that there is nothing that will cause problems in the chain
  mylabels.resize(0); getActionThatCalculates()->getAllActionLabelsInChain( mylabels );
  for(unsigned i=0;i<mylabels.size();++i) {
      ActionWithVector* av1=plumed.getActionSet().selectWithLabel<ActionWithVector*>( mylabels[i] ); 
      for(unsigned j=0;j<i;++j) {
          ActionWithVector* av2=plumed.getActionSet().selectWithLabel<ActionWithVector*>( mylabels[j] );
          if( !av1->canBeAfterInChain( av2 ) ) error("must calculate " + mylabels[j] + " before " + mylabels[i] );
      }
  }
  action_to_do_after=act;
  act->action_to_do_before=this;
  return true;
}

ActionWithVector* ActionWithVector::getActionThatCalculates() {
  // Return this if we have no dependencies
  if( !action_to_do_before ) return this;
  // Recursively go through actiosn before
  return action_to_do_before->getActionThatCalculates();
}

void ActionWithVector::getAllActionLabelsInChain( std::vector<std::string>& mylabels ) const {
  bool found = false ;
  for(unsigned i=0; i<mylabels.size(); ++i) {
    if( getLabel()==mylabels[i] ) { found=true; }
  }
  if( !found ) mylabels.push_back( getLabel() );
  if( action_to_do_after ) action_to_do_after->getAllActionLabelsInChain( mylabels );
}

void ActionWithVector::runAllTasks() {
// Skip this if this is done elsewhere
  if( action_to_do_before ) return;
      
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(serial) { stride=1; rank=0; }
        
  // Get the number of tasks
  unsigned ntasks=0; getNumberOfTasks( ntasks );
  // Determine if some tasks can be switched off
  // setTaskFlags( ntasks, taskSet );
  // Get the number of tasks
  // nactive_tasks = taskSet.size();
  nactive_tasks = ntasks;
  // Create a vector from the task set 
  // std::vector<AtomNumber> partialTaskList( taskSet.begin(), taskSet.end() );
  // Get number of threads for OpenMP
  unsigned nt=OpenMP::getNumThreads();
  if( nt*stride*10>nactive_tasks ) nt=nactive_tasks/stride/10;
  if( nt==0 ) nt=1;
  
  // Now do all preparations required to run all the tasks
  // prepareForTaskLoop();
  
  // Get the total number of streamed quantities that we need
  unsigned nquantities = 0, ncols=0, nmatrices=0;
  getNumberOfStreamedQuantities( nquantities, ncols, nmatrices );
  // Get size for buffer
  unsigned bufsize=0; getSizeOfBuffer( nactive_tasks, bufsize );
  if( buffer.size()!=bufsize ) buffer.resize( bufsize );
  // Clear buffer
  buffer.assign( buffer.size(), 0.0 ); 
  
  // Recover the number of derivatives we require
  unsigned nderivatives = 0; bool gridsInStream=checkForGrids();
  if( !doNotCalculateDerivatives() || gridsInStream ) getNumberOfStreamedDerivatives( nderivatives );

  #pragma omp parallel num_threads(nt)
  {
    std::vector<double> omp_buffer;
    if( nt>1 ) omp_buffer.resize( bufsize, 0.0 );
    MultiValue myvals( nquantities, nderivatives, ncols, nmatrices );
    myvals.clearAll();

    #pragma omp for nowait
    for(unsigned i=rank; i<nactive_tasks; i+=stride) {
      // Calculate the stuff in the loop for this action
      // runTask( partialTaskList[i].index(), myvals );
      runTask( i, myvals );

      // Now transfer the data to the actions that accumulate values from the calculated quantities
      if( nt>1 ) {
        // gatherAccumulators( partialTaskList[i].index(), myvals, omp_buffer );
        gatherAccumulators( i, myvals, omp_buffer );
      } else {
        // gatherAccumulators( partialTaskList[i].index(), myvals, buffer );
        gatherAccumulators( i, myvals, buffer );
      }

      // Clear the value
      myvals.clearAll();
    }
    #pragma omp critical
    if(nt>1) for(unsigned i=0; i<bufsize; ++i) buffer[i]+=omp_buffer[i];
  }
  // MPI Gather everything
  if( !serial && buffer.size()>0 ) comm.Sum( buffer );
  finishComputations( buffer );
}

unsigned ActionWithVector::getArgumentPositionInStream( const unsigned& jder, MultiValue& myvals ) const { 
  if( getPntrToArgument(jder)->storedata ) plumed_merror("You still have to implement setting these values Gareth");
  return getPntrToArgument(jder)->getPositionInStream();
}

bool ActionWithVector::checkForGrids() const {
  for(int i=0; i<getNumberOfComponents(); ++i) {
    if( getConstPntrToComponent(i)->getRank()>0 && getConstPntrToComponent(i)->hasDerivatives() ) return true;
  }
  if( action_to_do_after ) return action_to_do_after->checkForGrids();
  return false;
}

void ActionWithVector::getNumberOfTasks( unsigned& ntasks ) { 
  if( ntasks==0 ) { 
      plumed_assert( getNumberOfComponents()>0 && getPntrToComponent(0)->getRank()>0 ); 
      if( getPntrToComponent(0)->hasDerivatives() ) ntasks = getPntrToComponent(0)->getNumberOfValues();
      else ntasks = getPntrToComponent(0)->getShape()[0]; 
  }
  for(int i=0; i<getNumberOfComponents(); ++i) {
      if( getPntrToComponent(i)->getRank()==0 ) { 
          if( getNumberOfArguments()>1  || ntasks!=getPntrToArgument(0)->getNumberOfValues() ) error("mismatched numbers of tasks in streamed quantities");
      } else if( getPntrToComponent(i)->hasDerivatives() && ntasks!=getPntrToComponent(i)->getNumberOfValues() ) error("mismatched numbers of tasks in streamed quantities");
      else if( ntasks!=getPntrToComponent(i)->getShape()[0] ) error("mismatched numbers of tasks in streamed quantities");
  }
  if( action_to_do_after ) action_to_do_after->getNumberOfTasks( ntasks );
}

void ActionWithVector::getNumberOfStreamedQuantities( unsigned& nquants, unsigned& ncols, unsigned& nmat ) {
  for(unsigned i=0; i<getNumberOfArguments(); ++i) { getPntrToArgument(i)->streampos=nquants; nquants++; }
  for(int i=0; i<getNumberOfComponents(); ++i) { getPntrToComponent(i)->streampos=nquants; nquants++; }
  if( action_to_do_after ) action_to_do_after->getNumberOfStreamedQuantities( nquants, ncols, nmat );
}

void ActionWithVector::getSizeOfBuffer( const unsigned& nactive_tasks, unsigned& bufsize ) {
  for(int i=0; i<getNumberOfComponents(); ++i) { getPntrToComponent(i)->bufstart=bufsize; bufsize += getPntrToComponent(i)->data.size(); }
  if( action_to_do_after ) action_to_do_after->getSizeOfBuffer( nactive_tasks, bufsize );
}

void ActionWithVector::getNumberOfStreamedDerivatives( unsigned& nderivatives ) {
  unsigned nnd=0;
  if( !doNotCalculateDerivatives() ) {
    nnd = getNumberOfDerivatives();
  } else {
    for(int i=0; i<getNumberOfComponents(); ++i) {
      if( getConstPntrToComponent(i)->getRank()>0 && getConstPntrToComponent(i)->hasDerivatives() ) { nnd = getNumberOfDerivatives(); break; }
    }
  }
  if( nnd>nderivatives ) nderivatives = nnd;
  if( action_to_do_after ) action_to_do_after->getNumberOfStreamedDerivatives( nderivatives );
} 

void ActionWithVector::runTask( const unsigned& current, MultiValue& myvals ) const {
  if( isActive() ) {
    myvals.setTaskIndex(current); myvals.vector_call=true; performTask( current, myvals );
  }
  if( action_to_do_after ) action_to_do_after->runTask( current, myvals );
}

void ActionWithVector::gatherAccumulators( const unsigned& taskCode, const MultiValue& myvals, std::vector<double>& buffer ) const {
  if( isActive() ) {
    for(int i=0; i<getNumberOfComponents(); ++i) {
      unsigned bufstart = getConstPntrToComponent(i)->bufstart;
      // This looks after storing of scalars that are summed from vectors/matrices
      if( getConstPntrToComponent(i)->getRank()==0 ) {
        plumed_dbg_massert( bufstart<buffer.size(), "problem in " + getLabel() );
        unsigned sind = getConstPntrToComponent(i)->streampos; buffer[bufstart] += myvals.get(sind);
        if( getConstPntrToComponent(i)->hasDerivatives() ) {
          for(unsigned k=0; k<myvals.getNumberActive(sind); ++k) {
            unsigned kindex = myvals.getActiveIndex(sind,k);
            plumed_dbg_massert( bufstart+1+kindex<buffer.size(), "problem in " + getLabel()  );
            buffer[bufstart + 1 + kindex] += myvals.getDerivative(sind,kindex);
          }
        }
      // This looks after storing of vectors
      } else if( getConstPntrToComponent(i)->storedata ) {
        plumed_dbg_assert( getConstPntrToComponent(i)->getRank()==1 && !getConstPntrToComponent(i)->hasDeriv );
        unsigned vindex = getConstPntrToComponent(i)->bufstart + taskCode; plumed_dbg_massert( vindex<buffer.size(), "failing in " + getLabel() );
        buffer[vindex] += myvals.get(getConstPntrToComponent(i)->streampos);
      }
    }
  }
  if( action_to_do_after ) action_to_do_after->gatherAccumulators( taskCode, myvals, buffer ); 
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
      } else if( getPntrToComponent(i)->getRank()==0 ) getPntrToComponent(i)->set( buf[bufstart] );
      // This gathers derivatives of scalars
      if( !doNotCalculateDerivatives() && getPntrToComponent(i)->hasDeriv && getPntrToComponent(i)->getRank()==0 ) {
        for(unsigned j=0; j<getPntrToComponent(i)->getNumberOfDerivatives(); ++j) getPntrToComponent(i)->setDerivative( j, buf[bufstart+1+j] );
      }
    }
    transformFinalValueAndDerivatives( buf ); 
  } 
  if( action_to_do_after ) action_to_do_after->finishComputations( buf );
}

}
