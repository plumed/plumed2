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
  action_to_do_after(NULL)
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

bool ActionWithVector::checkForGrids() const {
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
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
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
      if( getPntrToComponent(i)->hasDerivatives() && ntasks!=getPntrToComponent(i)->getNumberOfValues() ) error("mismatched numbers of tasks in streamed quantities");
      else if( ntasks!=getPntrToComponent(i)->getShape()[0] ) error("mismatched numbers of tasks in streamed quantities");
  }
  if( action_to_do_after ) action_to_do_after->getNumberOfTasks( ntasks );
}

void ActionWithVector::getNumberOfStreamedQuantities( unsigned& nquants, unsigned& ncols, unsigned& nmat ) {
  for(unsigned i=0; i<getNumberOfArguments(); ++i) { getPntrToArgument(i)->streampos=nquants; nquants++; }
  for(unsigned i=0; i<getNumberOfComponents(); ++i) { getPntrToComponent(i)->streampos=nquants; nquants++; }
  if( action_to_do_after ) action_to_do_after->getNumberOfStreamedQuantities( nquants, ncols, nmat );
}

void ActionWithVector::getSizeOfBuffer( const unsigned& nactive_tasks, unsigned& bufsize ) {
  for(unsigned i=0; i<getNumberOfComponents(); ++i) { getPntrToComponent(i)->bufstart=bufsize; bufsize += getPntrToComponent(i)->data.size(); }
  if( action_to_do_after ) action_to_do_after->getSizeOfBuffer( nactive_tasks, bufsize );
}

void ActionWithVector::getNumberOfStreamedDerivatives( unsigned& nderivatives ) {
  unsigned nnd=0;
  if( !doNotCalculateDerivatives() ) {
    nnd = getNumberOfDerivatives();
  } else {
    for(unsigned i=0; i<getNumberOfComponents(); ++i) {
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
    for(unsigned i=0; i<getNumberOfComponents(); ++i) {
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
    for(unsigned i=0; i<getNumberOfComponents(); ++i) {
      // This gathers vectors and grids at the end of the calculation
      unsigned bufstart = getPntrToComponent(i)->bufstart;
      getPntrToComponent(i)->data.assign( getPntrToComponent(i)->data.size(), 0 );
      if( (getPntrToComponent(i)->getRank()>0 && getPntrToComponent(i)->hasDerivatives()) || getPntrToComponent(i)->storedata ) {
        unsigned sz_v = getPntrToComponent(i)->data.size();
        for(unsigned j=0; j<sz_v; ++j) getPntrToComponent(i)->add( j, buffer[bufstart+j] );
      }
      // This gathers derivatives of scalars
      if( !doNotCalculateDerivatives() && getPntrToComponent(i)->hasDeriv && getPntrToComponent(i)->getRank()==0 ) {
        for(unsigned j=0; j<getPntrToComponent(i)->getNumberOfDerivatives(); ++j) getPntrToComponent(i)->setDerivative( j, buffer[bufstart+1+j] );
      }
    }
    transformFinalValueAndDerivatives( buffer );
  } 
  if( action_to_do_after ) action_to_do_after->finishComputations( buffer );
}

}
