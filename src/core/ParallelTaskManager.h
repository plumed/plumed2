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

template <class D>
class ParallelActionsInput {
public:
/// Do we need to calculate the derivatives
  bool noderiv;
/// Periodic boundary conditions
  const Pbc& pbc;
/// The number of components the underlying action is computing
  unsigned ncomponents;
//// The number of indexes that you use per task
  unsigned nindices_per_task;
/// This holds all the input data that is required to calculate all values for all tasks
  std::vector<double> inputdata;
/// This holds data for that the underlying action needs to do the calculation
  D actiondata;
/// Default constructor
  ParallelActionsInput( const Pbc& box ) : noderiv(false), pbc(box), ncomponents(0), nindices_per_task(0) {}
};

template <class T, class D>
class ParallelTaskManager {
private:
/// The underlying action for which we are managing parallel tasks
  ActionWithVector* action;
/// The MPI communicator
  Communicator& comm;
/// Is this an action with matrix
  bool ismatrix;
/// Are we using acc for parallisation
  bool useacc;
/// This holds the values before we pass them to the value
  Matrix<double> value_mat;
/// A tempory set of vectors for holding forces over threads
  std::vector<std::vector<double> > omp_forces;
/// An action to hold data that we pass to and from the static function
  ParallelActionsInput<D> myinput;
public:
  static void registerKeywords( Keywords& keys );
  ParallelTaskManager(ActionWithVector* av);
/// Setup an array to hold all the indices that are used for derivatives
  void setNumberOfIndicesPerTask( const std::size_t& nind );
/// Copy the data from the underlying colvar into this parallel action
  void setActionInput( const D& adata );
/// This runs all the tasks
  void runAllTasks();
/// Apply the forces on the parallel object
  void applyForces( std::vector<double>& forcesForApply );
};

template <class T, class D>
void ParallelTaskManager<T, D>::registerKeywords( Keywords& keys ) {
  keys.addFlag("USEGPU",false,"run this calculation on the GPU");
}

template <class T, class D>
ParallelTaskManager<T, D>::ParallelTaskManager(ActionWithVector* av):
  action(av),
  comm(av->comm),
  ismatrix(false),
  useacc(false),
  myinput(av->getPbc())
{
  ActionWithMatrix* am=dynamic_cast<ActionWithMatrix*>(av);
  if(am) ismatrix=true;
  action->parseFlag("USEGPU",useacc);
#ifdef __PLUMED_HAS_OPENACC
  if( useacc ) action->log.printf("  using GPU to calculate this action\n");
#else
  if( useacc ) action->error("cannot use USEGPU flag as PLUMED has not been compiled with openacc");
#endif
}

template <class T, class D>
void ParallelTaskManager<T, D>::setNumberOfIndicesPerTask( const std::size_t& nind ) {
  plumed_massert( action->getNumberOfComponents()>0, "there should be some components wen you setup the index list" );
  std::size_t valuesize=(action->getConstPntrToComponent(0))->getNumberOfStoredValues();
  for(unsigned i=1; i<action->getNumberOfComponents(); ++i) plumed_assert( valuesize==(action->getConstPntrToComponent(i))->getNumberOfStoredValues() );
  myinput.ncomponents = action->getNumberOfComponents();
  value_mat.resize( valuesize, action->getNumberOfComponents() ); myinput.nindices_per_task = nind;
}

template <class T, class D>
void ParallelTaskManager<T, D>::setActionInput( const D& adata ) {
  myinput.actiondata=adata;
}

template <class T, class D>
void ParallelTaskManager<T, D>::runAllTasks() {
  // Get the list of active tasks
  std::vector<unsigned> & partialTaskList( action->getListOfActiveTasks( action ) );
  unsigned nactive_tasks=partialTaskList.size();
  // Clear the value matrix
  value_mat = 0;

  // Get all the input data so we can broadcast it to the GPU
  myinput.noderiv = true;
  action->getInputData( myinput.inputdata );

  if( useacc ) {
#ifdef __PLUMED_HAS_OPENACC
     #pragma acc data copyin(nactive_tasks) copyin(partialTaskList) copyin(myinput) copy(value_mat)
       {
     #pragma acc parallel loop
         for(unsigned i=0; i<nactive_tasks; ++i) {
           // Calculate the stuff in the loop for this action
           const auto [values, derivs] = T::performTask( partialTaskList[i], myinput );
     
           // Transfer the data to the values
           T::transferToValue( partialTaskList[i], values, value_mat );
         }
       }
#else
       plumed_merror("cannot use USEGPU flag if PLUMED has not been compiled with openACC");
#endif
  } else {
     // Get the MPI details
     unsigned stride=comm.Get_size();
     unsigned rank=comm.Get_rank();
     if(action->runInSerial()) { stride=1; rank=0; }

     // Get number of threads for OpenMP
     unsigned nt=OpenMP::getNumThreads();
     if( nt*stride*10>nactive_tasks ) nt=nactive_tasks/stride/10;
     if( nt==0 ) nt=1;

     #pragma omp parallel num_threads(nt)
     {
       #pragma omp for nowait
       for(unsigned i=rank; i<nactive_tasks; i+=stride) {
         // Calculate the stuff in the loop for this action
         const auto [values, derivs] = T::performTask( partialTaskList[i], myinput );

         // Transfer the data to the values
         T::transferToValue( partialTaskList[i], values, value_mat );
       }
     }
     // MPI Gather everything
     if( !action->runInSerial() ) comm.Sum( value_mat );
  }

  // And transfer the value to the output values
  for(unsigned i=0; i<action->getNumberOfComponents(); ++i) {
    Value* myval = action->copyOutput(i);
    for(unsigned j=0; j<myval->getNumberOfStoredValues(); ++j) myval->set( j, value_mat[j][i] );
  }
}

template <class T, class D>
void ParallelTaskManager<T, D>::applyForces( std::vector<double>& forcesForApply ) {
  // Get the list of active tasks
  std::vector<unsigned> & partialTaskList( action->getListOfActiveTasks( action ) );
  unsigned nactive_tasks=partialTaskList.size();
  // Clear force buffer
  forcesForApply.assign( forcesForApply.size(), 0.0 );

  // Get all the input data so we can broadcast it to the GPU
  myinput.noderiv = false;
  // Retrieve the forces from the values
  for(unsigned i=0; i<action->getNumberOfComponents(); ++i) {
    Value* myval = action->copyOutput(i);
    for(unsigned j=0; j<myval->getNumberOfStoredValues(); ++j) value_mat[j][i] = myval->getForce( j );
  }

  if( useacc ) {
#ifdef __PLUMED_HAS_OPENACC
     #pragma acc data copyin(nactive_tasks) copyin(partialTaskList) copyin(myinput) copyin(value_mat) copy(forcesForApply)
       {
     #pragma acc parallel loop reduction(forcesForApply)
         for(unsigned i=0; i<nactive_tasks; ++i) {
           // Calculate the stuff in the loop for this action
           const auto [values, derivs] = T::performTask( partialTaskList[i], myinput );
     
           // Gather the forces from the values
           T::gatherForces( partialTaskList[i], myinput, value_mat, derivs, forcesForApply );
         }
       }
#else
       plumed_merror("cannot use USEGPU flag if PLUMED has not been compiled with openACC");
#endif
  } else {
     // Get the MPI details
     unsigned stride=comm.Get_size();
     unsigned rank=comm.Get_rank();
     if(action->runInSerial()) { stride=1; rank=0; }

     // Get number of threads for OpenMP
     unsigned nt=OpenMP::getNumThreads();
     if( nt*stride*10>nactive_tasks ) nt=nactive_tasks/stride/10;
     if( nt==0 ) nt=1;

     #pragma omp parallel num_threads(nt)
     {
       const unsigned t=OpenMP::getThreadNum();
       if( nt>1 ) {
         if( omp_forces[t].size()!=forcesForApply.size() ) omp_forces[t].resize( forcesForApply.size(), 0.0 );
         else omp_forces[t].assign( forcesForApply.size(), 0.0 );
       }
       #pragma omp for nowait
       for(unsigned i=rank; i<nactive_tasks; i+=stride) {
         // Calculate the stuff in the loop for this action
         const auto [values, derivs] = T::performTask( partialTaskList[i], myinput );

         // Gather the forces from the values
         if( nt>1 ) T::gatherForces( partialTaskList[i], myinput, value_mat, derivs, omp_forces[t] );
         else T::gatherForces( partialTaskList[i], myinput, value_mat, derivs, forcesForApply );
       }
       #pragma omp critical
       if(nt>1) for(unsigned i=0; i<forcesForApply.size(); ++i) forcesForApply[i]+=omp_forces[t][i];
     }
     // MPI Gather everything
     if( !action->runInSerial() ) comm.Sum( forcesForApply );
  }

}

}
#endif
