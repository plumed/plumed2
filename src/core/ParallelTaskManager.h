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

#include "ActionWithVector.h"
#include "ActionWithMatrix.h"
#include "tools/Communicator.h"
#include "tools/OpenMP.h"
#include "tools/View.h"
#include "tools/View2D.h"

namespace PLMD {
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
/// Default constructor
  ParallelActionsInput( const Pbc& box ) : noderiv(false), pbc(box), ncomponents(0), nindices_per_task(0) {}
};

class ParallelActionsOutput {
public:
  View<double,helpers::dynamic_extent> values;
  std::vector<double>& derivatives;
  ParallelActionsOutput( std::size_t ncomp, double* v, std::vector<double>& d ) : values(v,ncomp), derivatives(d) {}
};

class ForceInput {
private:
  std::size_t ncomp;
  class DerivHelper {
  private:
    std::size_t nderiv;
    std::vector<double>& derivatives;
  public:
    DerivHelper( std::size_t n, std::vector<double>& d ) : nderiv(n), derivatives(d) {}
    View<double, helpers::dynamic_extent> operator[](std::size_t i) const {
      return View<double, helpers::dynamic_extent>( derivatives.data() + i*nderiv, nderiv );
    }
  };
public:
  View<double,helpers::dynamic_extent> force;
  DerivHelper deriv;
  ForceInput( std::size_t n, double* f, std::size_t m, std::vector<double>& d ) : ncomp(n), force(f,n), deriv(m,d) {}
};

struct ForceOutput {
  std::vector<double>& thread_safe;
  std::vector<double>& thread_unsafe;
  ForceOutput( std::vector<double>& s, std::vector<double>& u ) : thread_safe(s), thread_unsafe(u) {}
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
/// Number of derivatives for each task
  std::size_t nderivatives_per_task;
/// The number of forces on each thread
  std::size_t nthreaded_forces;
/// This holds the values before we pass them to the value
  std::vector<double> value_stash;
/// A tempory set of vectors for holding forces over threads
  std::vector<std::vector<double> > omp_forces;
/// An action to hold data that we pass to and from the static function
  ParallelActionsInput myinput;
/// This holds data for that the underlying action needs to do the calculation
  D actiondata;
public:
  static void registerKeywords( Keywords& keys );
  ParallelTaskManager(ActionWithVector* av);
/// Setup the derivatives
  void setNumberOfIndicesAndDerivativesPerTask( const std::size_t& nind, const std::size_t& nder );
/// Set the number of forces that are gathered over threads
  void setNumberOfThreadedForces( std::size_t nt );
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
  nderivatives_per_task(0),
  nthreaded_forces(0),
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
void ParallelTaskManager<T, D>::setNumberOfIndicesAndDerivativesPerTask( const std::size_t& nind, const std::size_t& nder ) {
  plumed_massert( action->getNumberOfComponents()>0, "there should be some components wen you setup the index list" );
  std::size_t valuesize=(action->getConstPntrToComponent(0))->getNumberOfStoredValues();
  for(unsigned i=1; i<action->getNumberOfComponents(); ++i) plumed_assert( valuesize==(action->getConstPntrToComponent(i))->getNumberOfStoredValues() );
  myinput.ncomponents = action->getNumberOfComponents(); nderivatives_per_task = nder;
  value_stash.resize( valuesize*action->getNumberOfComponents() ); myinput.nindices_per_task = nind;
}

template <class T, class D>
void ParallelTaskManager<T, D>::setNumberOfThreadedForces( std::size_t nt ) {
  nthreaded_forces = nt; unsigned t=OpenMP::getNumThreads(); if( useacc ) t = 1;
  omp_forces.resize(t); for(unsigned i=0; i<t; ++i) omp_forces.resize(nt);
}

template <class T, class D>
void ParallelTaskManager<T, D>::setActionInput( const D& adata ) {
  actiondata=adata;
}

template <class T, class D>
void ParallelTaskManager<T, D>::runAllTasks() {
  // Get the list of active tasks
  std::vector<unsigned> & partialTaskList( action->getListOfActiveTasks( action ) );
  unsigned nactive_tasks=partialTaskList.size();

  // Get all the input data so we can broadcast it to the GPU
  myinput.noderiv = true;
  action->getInputData( myinput.inputdata );
  value_stash.assign( value_stash.size(), 0.0 );

  if( useacc ) {
    std::size_t indatasize = myinput.inputdata.size();
#ifdef __PLUMED_HAS_OPENACC
#pragma acc data copyin(nactive_tasks) copyin(partialTaskList) copyin(myinput) copyin(myinput.inputdata[0:indatasize]) copy(value_stash)
    {
      std::vector<double> derivatives( myinput.ncomponents*nderivatives_per_task );
#pragma acc parallel loop
      for(unsigned i=0; i<nactive_tasks; ++i) {
        std::size_t task_number = partialTaskList[i];
        std::size_t val_pos = task_number*myinput.ncomponents;
        ParallelActionsOutput myout( myinput.ncomponents, value_stash.data()+val_pos, derivatives );
        // Calculate the stuff in the loop for this action
        T::performTask( task_number, actiondata, myinput, myout );
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
      std::vector<double> derivatives( myinput.ncomponents*nderivatives_per_task );
      #pragma omp for nowait
      for(unsigned i=rank; i<nactive_tasks; i+=stride) {
        std::size_t task_number = partialTaskList[i];
        std::size_t val_pos = task_number*myinput.ncomponents;
        ParallelActionsOutput myout( myinput.ncomponents, value_stash.data()+val_pos, derivatives );
        // Calculate the stuff in the loop for this action
        T::performTask( task_number, actiondata, myinput, myout );
      }
    }
    // MPI Gather everything
    if( !action->runInSerial() ) comm.Sum( value_stash );
  }

  // And transfer the value to the output values
  for(unsigned i=0; i<action->getNumberOfComponents(); ++i) {
    Value* myval = action->copyOutput(i);
    for(unsigned j=0; j<myval->getNumberOfStoredValues(); ++j) myval->set( j, value_stash[j*myinput.ncomponents+i] );
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
    for(unsigned j=0; j<myval->getNumberOfStoredValues(); ++j) value_stash[j*myinput.ncomponents+i] = myval->getForce( j );
  }

  if( useacc ) {
#ifdef __PLUMED_HAS_OPENACC
    omp_forces[0].assign( omp_forces[0].size(), 0.0 );

#pragma acc data copyin(nactive_tasks) copyin(partialTaskList) copyin(myinput) copyin(value_stash) copy(omp_forces[0]) copy(forcesForApply)
    {
#pragma acc parallel loop reduction(omp_forces)
      ForceOutput forces( omp_forces[0], forcesForApply );
      std::vector<double> fake_vals( myinput.ncomponents );
      std::vector<double> derivatives( myinput.ncomponents*nderivatives_per_task );
      for(unsigned i=0; i<nactive_tasks; ++i) {
        std::size_t task_index = partialTaskList[i];
        ParallelActionsOutput myout( myinput.ncomponents, fake_vals.data(), derivatives );
        // Calculate the stuff in the loop for this action
        T::performTask( task_index, actiondata, myinput, myout );

        // Gather the forces from the values
        T::gatherForces( task_index, actiondata, myinput, ForceInput( myinput.ncomponents, value_stash.data()+myinput.ncomponents*task_index, nderivatives_per_task, derivatives), forces );
      }
      T::gatherThreads( omp_forces[0], forcesForApply );
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
      omp_forces[t].assign( nthreaded_forces, 0.0 );
      ForceOutput forces( omp_forces[t], forcesForApply );
      std::vector<double> fake_vals( myinput.ncomponents );
      std::vector<double> derivatives( myinput.ncomponents*nderivatives_per_task );
      #pragma omp for nowait
      for(unsigned i=rank; i<nactive_tasks; i+=stride) {
        std::size_t task_index = partialTaskList[i];
        ParallelActionsOutput myout( myinput.ncomponents, fake_vals.data(), derivatives );
        // Calculate the stuff in the loop for this action
        T::performTask( task_index, actiondata, myinput, myout );

        // Gather the forces from the values
        T::gatherForces( task_index, actiondata, myinput, ForceInput( myinput.ncomponents, value_stash.data()+myinput.ncomponents*task_index, nderivatives_per_task, derivatives), forces );
      }
      #pragma omp critical
      T::gatherThreads( forces );
    }
    // MPI Gather everything
    if( !action->runInSerial() ) comm.Sum( forcesForApply );
  }

}

}
#endif
