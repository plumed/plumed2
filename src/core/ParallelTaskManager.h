/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2023 The plumed team
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
#ifndef __PLUMED_core_ParallelTaskManager_h
#define __PLUMED_core_ParallelTaskManager_h

#include "tools/Communicator.h"
#include "ActionWithVector.h"
#include "ActionWithMatrix.h"
#include "tools/OpenMP.h"

namespace PLMD {

template <class T>
class ParallelTaskManager {
private:
/// The underlying action for which we are managing parallel tasks
  ActionWithVector* action;
/// The MPI communicator
  Communicator& comm;
/// Is this an action with matrix
  bool ismatrix;
/// The buffer that we use (we keep a copy here to avoid resizing)
  std::vector<double> buffer;
/// A tempory vector of MultiValue so we can avoid doing lots of resizes
  std::vector<MultiValue> myvals;
/// This holds indices for creating derivatives
  std::vector<std::size_t> indices;
/// This holds all the input data that is required to calculate all values for all tasks
  std::vector<double> inputdata;
public:
  ParallelTaskManager(ActionWithVector* av);
/// Setup an array to hold all the indices that are used for derivatives
  void setupIndexList( const std::vector<std::size_t>& ind );
/// This runs all the tasks
  void runAllTasks( const unsigned& natoms=0 );
/// This runs each of the tasks
  void runTask( const unsigned& current, MultiValue& myvals ) const ;
};

template <class T>
ParallelTaskManager<T>::ParallelTaskManager(ActionWithVector* av):
  action(av),
  comm(av->comm),
  ismatrix(false)
{
  ActionWithMatrix* am=dynamic_cast<ActionWithMatrix*>(av);
  if(am) ismatrix=true;
}

template <class T>
void ParallelTaskManager<T>::setupIndexList( const std::vector<std::size_t>& ind ) {
  indices.resize( ind.size() ); for(unsigned i=0; i<ind.size(); ++i) indices[i] = ind[i];
}

template <class T>
void ParallelTaskManager<T>::runAllTasks( const unsigned& natoms ) {
  unsigned stride=comm.Get_size();
  unsigned rank=comm.Get_rank();
  if(action->runInSerial()) { stride=1; rank=0; }

  // Clear matrix bookeeping arrays
  // if( ismatrix && stride>1 ) clearMatrixBookeeping();

  // Get the list of active tasks
  std::vector<unsigned> & partialTaskList( action->getListOfActiveTasks( action ) );
  unsigned nactive_tasks=partialTaskList.size();

  // Get number of threads for OpenMP
  unsigned nt=OpenMP::getNumThreads();
  if( nt*stride*10>nactive_tasks ) nt=nactive_tasks/stride/10;
  if( nt==0 ) nt=1;
  if( myvals.size()!=nt ) myvals.resize(nt);

  // Get the total number of streamed quantities that we need
  // Get size for buffer
  unsigned bufsize=0, nderivatives = 0; bool gridsInStream=false;
  if( buffer.size()!=bufsize ) buffer.resize( bufsize );
  // Clear buffer
  buffer.assign( buffer.size(), 0.0 );

  // Get all the input data so we can broadcast it to the GPU
  action->getInputData( inputdata );

  #pragma omp parallel num_threads(nt)
  {
    std::vector<double> omp_buffer;
    const unsigned t=OpenMP::getThreadNum();
    if( nt>1 ) omp_buffer.resize( bufsize, 0.0 );
    if( myvals[t].getNumberOfValues()!=action->getNumberOfComponents() || myvals[t].getNumberOfDerivatives()!=nderivatives || myvals[t].getAtomVector().size()!=natoms ) {
      myvals[t].resize( action->getNumberOfComponents(), nderivatives, natoms );
    }
    myvals[t].clearAll();

    #pragma omp for nowait
    for(unsigned i=rank; i<nactive_tasks; i+=stride) {
      // Calculate the stuff in the loop for this action
      runTask( partialTaskList[i], myvals[t] );

      // Clear the value
      myvals[t].clearAll();
    }
    #pragma omp critical
    if( nt>1 ) for(unsigned i=0; i<bufsize; ++i) buffer[i]+=omp_buffer[i];
  }
  // MPI Gather everything
  if( !action->runInSerial() ) {
    if( buffer.size()>0 ) comm.Sum( buffer );
    for(unsigned i=0; i<action->getNumberOfComponents(); ++i) (action->copyOutput(i))->MPIGatherTasks( !ismatrix, comm );
  }
}

template <class T>
void ParallelTaskManager<T>::runTask( const unsigned& current, MultiValue& myvals ) const {
  const ActionWithMatrix* am = dynamic_cast<const ActionWithMatrix*>(action);
  if( am ) return ;
  myvals.setTaskIndex(current); action->performTask( current, myvals );
  for(unsigned i=0; i<action->getNumberOfComponents(); ++i) {
    const Value* myval = action->getConstPntrToComponent(i);
    if( myval->hasDerivatives() || (action->getName()=="RMSD_VECTOR" && myval->getRank()==2) ) continue;
    Value* myv = const_cast<Value*>( myval ); myv->set( current, myvals.get( i ) );
  }
}

}
#endif
