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

#ifdef __PLUMED_USE_OPENACC
#define __PLUMED_HAS_OPENACC 1
#endif //__PLUMED_USE_OPENACC

#include "ActionWithVector.h"
#include "ActionWithMatrix.h"
#include "tools/Communicator.h"
#include "tools/OpenMP.h"
#include "tools/View.h"
#include "tools/View2D.h"

#include<iostream>
#define vdbg(...) std::cerr << __LINE__ << ":" << #__VA_ARGS__ << " " << (__VA_ARGS__) << '\n'
namespace PLMD {
//this  little tool is a RAII helper to put and remove data on the device with less risks of foregging something
template<typename accData>
class fromToDataHelper {
  accData &m;
public:
  fromToDataHelper(accData &d) : m(d) { d.toACCDevice(); }
  ~fromToDataHelper() { m.removeFromACCDevice(); }
};

struct ParallelActionsInput {
/// Do we need to calculate the derivatives
  bool noderiv{false};
/// Periodic boundary conditions
  const Pbc* pbc;
/// The number of components the underlying action is computing
  unsigned ncomponents{0};
//// The number of indexes that you use per task
  unsigned nindices_per_task{0};
/// This holds all the input data that is required to calculate all values for all tasks
  unsigned dataSize{0};
  double *inputdata{nullptr};
  // std::vector<double> inputdata{};
/// Default constructor
  ParallelActionsInput( const Pbc& box ) : pbc(&box) {}
  /// the copy is needed due to some openACC resistance to pass the inputdata as a std::vector (and some "this" problems, please do not ask, I'll get grumpy)
  ParallelActionsInput( const ParallelActionsInput& other ) :
    noderiv {other.noderiv},
    //just copies the ptr
    pbc  {other.pbc},
    ncomponents {other.ncomponents},
    nindices_per_task {other.nindices_per_task},
    dataSize {other.dataSize},
    //just copies the value of the pointer
    inputdata  {other.inputdata}
  {}

  ParallelActionsInput& operator=( const ParallelActionsInput& other )  = delete;

  //helper function to bring data to the device in a controlled way
  void toACCDevice()const {
#pragma acc enter data copyin(this[0:1], noderiv, pbc[0:1],ncomponents, nindices_per_task, dataSize, inputdata[0:dataSize])

    pbc->toACCDevice();
  }
  //helper function to remove data from the device in a controlled way
  void removeFromACCDevice() const  {
    pbc->removeFromACCDevice();
    // assuming dataSize is not changed
#pragma acc exit data delete(inputdata[0:dataSize],dataSize,nindices_per_task,ncomponents, pbc[0:1],noderiv,this[0:1])
  }
};

struct ParallelActionsOutput {
  View<double> values;
  View<double> derivatives;
  ParallelActionsOutput( std::size_t ncomp, double* v, std::size_t ndev, double* d ) : values(v,ncomp), derivatives(d,ndev) {}
};

class ForceInput {
private:
  std::size_t ncomp;
public:
  View<double> force;
  View2D<double> deriv;
  ForceInput( std::size_t nc, double* f, std::size_t nd, double* d ) : ncomp(nc), force(f,nc), deriv(d,nc,nd) {}
};

//There is no need to pass this as reference:
struct ForceOutput {
  //I would suggest to invert the name or to be clearer
// like something that recalls that "thread_safe" will be reducted
  double * thread_safe{nullptr};
  size_t const reducedSize{0};
  //these are the forces that we promise will not provoke races
  double *const thread_unsafe{nullptr};
  //const T* is a ptr to const T
  //T* const is a conts ptr to a modifiable T
  size_t  const other_size{0};
  ForceOutput(std::vector<double>& reduced, std::vector<double>& notReduced) :
    thread_safe{reduced.data()},
    reducedSize {reduced.size()},
    thread_unsafe{notReduced.data()},
    other_size{notReduced.size()}
  {}
  ForceOutput(double* reduced, size_t rs, double* notReduce, size_t osz):
    thread_safe(reduced),
    reducedSize(rs),
    thread_unsafe(notReduce),
    other_size(osz)
  {}
  ForceOutput(const ForceOutput &) = default;
  ForceOutput(ForceOutput &&) = default;
  ForceOutput& operator=(const ForceOutput &) = delete;
  ForceOutput& operator=(ForceOutput &&) = delete;
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
/// This structs is used to pass data between the parallel interface and the function caller
  ParallelActionsInput myinput;
//this holds the data for myinput that will be passed though myinput
  std::vector<double> input_buffer;
/// This holds data for that the underlying action needs to do the calculation
  D actiondata;
public:
  static void registerKeywords( Keywords& keys );
  ParallelTaskManager(ActionWithVector* av);
/// Setup the derivatives
  void setNumberOfIndicesAndDerivativesPerTask(std::size_t nind, std::size_t nder );
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
void ParallelTaskManager<T, D>::setNumberOfIndicesAndDerivativesPerTask( const std::size_t nind, const std::size_t nder ) {
  plumed_massert( action->getNumberOfComponents()>0, "there should be some components wen you setup the index list" );
  std::size_t valuesize=(action->getConstPntrToComponent(0))->getNumberOfStoredValues();
  for(unsigned i=1; i<action->getNumberOfComponents(); ++i) {
    plumed_assert( valuesize==(action->getConstPntrToComponent(i))->getNumberOfStoredValues() );
  }
  myinput.ncomponents = action->getNumberOfComponents();
  nderivatives_per_task = nder;
  value_stash.resize( valuesize*action->getNumberOfComponents() );
  myinput.nindices_per_task = nind;
}

template <class T, class D>
void ParallelTaskManager<T, D>::setNumberOfThreadedForces( std::size_t nt ) {
  //please change the name of this function to something more clear
  // because these are the "thread unsafe" forces that needs a reduction
  nthreaded_forces = nt;
  unsigned t=OpenMP::getNumThreads();
  if( useacc ) {
    t = 1;
  }
  omp_forces.resize(t);
  for(unsigned i=0; i<t; ++i) {
    omp_forces[i].resize(nt);
  }
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
  //shadowing myinput

  myinput.noderiv = true;
  action->getInputData( input_buffer );
  myinput.dataSize = input_buffer.size();
  myinput.inputdata = input_buffer.data();
  value_stash.assign( value_stash.size(), 0.0 );
  if( useacc ) {
#ifdef __PLUMED_HAS_OPENACC

    //I have a few problem with "this" <- meaning "this" pointer-  being copyed,
    // so I workarounded it with few copies
    //This is a RAII guard for moving data
    ParallelActionsInput input = myinput;
    auto myinput_acc = fromToDataHelper(input);
    D t_actiondata = actiondata;
    auto actiondata_acc = fromToDataHelper(t_actiondata);

    auto value_stash_data = value_stash.data();
    auto partialTaskList_data = partialTaskList.data();
    auto vs_size = value_stash.size();
    auto devSize = input.ncomponents*nderivatives_per_task;
#pragma acc parallel loop present(input, t_actiondata) \
                          copyin(nactive_tasks, devSize, \
                          partialTaskList_data[0:nactive_tasks])\
                          copy(value_stash_data[0:vs_size]) \
                          default(none)
    for(unsigned i=0; i<nactive_tasks; ++i) {
      std::vector<double> derivatives(devSize);
      std::size_t task_number = partialTaskList_data[i];
      std::size_t val_pos = task_number*input.ncomponents;
      ParallelActionsOutput myout( input.ncomponents,
                                   value_stash_data+val_pos,
                                   devSize,
                                   derivatives.data() );
      // Calculate the stuff in the loop for this action
      T::performTask( task_number, t_actiondata, input, myout );
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
        ParallelActionsOutput myout( myinput.ncomponents,
                                     value_stash.data()+val_pos,
                                     myinput.ncomponents*nderivatives_per_task,
                                     derivatives.data() );
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
    for(unsigned j=0; j<myval->getNumberOfStoredValues(); ++j) {
      myval->set( j, value_stash[j*myinput.ncomponents+i] );
    }
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
    auto myval = action->copyOutput(i);
    //will be this okeier?
    // auto myval = action->getConstPntrToComponent(i);
    for(unsigned j=0; j<myval->getNumberOfStoredValues(); ++j) {
      value_stash[j*myinput.ncomponents+i] = myval->getForce( j );
    }
  }

  if( useacc ) {
#ifdef __PLUMED_HAS_OPENACC
    omp_forces[0].assign( omp_forces[0].size(), 0.0 );
    auto value_stash_data = value_stash.data();
    const auto value_stash_size = value_stash.size();
    auto partialTaskList_data = partialTaskList.data();
    auto forcesForApply_data = forcesForApply.data();
    const auto ffa_size = forcesForApply.size();
    auto omp_forces_data = omp_forces[0].data();
    auto of_size = omp_forces[0].size();
    //on the cpu we only need the declaration, see the create statement below
    double* fake_vals;
    const auto nderivs =  myinput.ncomponents*nderivatives_per_task;
    const auto ndev_per_task = nderivatives_per_task;
    const auto ncomponents = myinput.ncomponents;

    //To future me/you:
    // I need to allocate this on the host to create a bigger temporay data array
    // on the device
    // by trying with double* x=nullptr, you will get a failure
    // another solution is acc_malloc and then device_ptr in the pragma
    // (but you have to remember the acc_free)
    std::vector<double> derivative(1);
    double * derivatives = derivative.data();

    ParallelActionsInput input = myinput;
    auto myinput_acc = fromToDataHelper(input);
    D t_actiondata = actiondata;
    auto actiondata_acc = fromToDataHelper(t_actiondata);


#pragma acc parallel loop reduction(+:omp_forces_data[0:of_size]) \
                            present(input,t_actiondata) \
                             copyin(nactive_tasks, \
                                    nderivs, \
                                    ndev_per_task, \
                                    partialTaskList_data[0:nactive_tasks], \
                                    value_stash_data[0:value_stash_size]) \
                               copy(omp_forces_data[0:of_size], \
                                    forcesForApply_data[0:ffa_size]) \
                             create(derivatives[0:nderivs]) \
                            default(none)
    for(unsigned i=0; i<nactive_tasks; ++i) {
      //This may be changed to a shared array
      std::vector<double> valstmp( input.ncomponents );
      std::size_t task_number = partialTaskList_data[i];
      std::size_t val_pos = task_number*input.ncomponents;
      // std::vector<double> derivative( nderivs );
      ParallelActionsOutput myout( input.ncomponents,
                                   valstmp.data(),
                                   nderivs,
                                   derivatives);
      // Calculate the stuff in the loop for this action
      T::performTask( task_number, t_actiondata, input, myout );

      // Gather the forces from the values
      T::gatherForces( task_number, t_actiondata, input,
                       ForceInput( input.ncomponents,
                                   value_stash_data+input.ncomponents*task_number,
                                   ndev_per_task,
                                   derivatives),
                       ForceOutput { omp_forces_data,of_size, forcesForApply_data,ffa_size }
                     );



    }
    T::gatherThreads( { omp_forces[0], forcesForApply});
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
      ForceOutput forces{ omp_forces[t], forcesForApply };
      std::vector<double> fake_vals( myinput.ncomponents );
      std::vector<double> derivatives( myinput.ncomponents*nderivatives_per_task );
      #pragma omp for nowait
      for(unsigned i=rank; i<nactive_tasks; i+=stride) {
        std::size_t task_index = partialTaskList[i];
        ParallelActionsOutput myout( myinput.ncomponents,
                                     fake_vals.data(),
                                     derivatives.size(),
                                     derivatives.data() );
        // Calculate the stuff in the loop for this action
        T::performTask( task_index, actiondata, myinput, myout );

        // Gather the forces from the values
        T::gatherForces( task_index, actiondata, myinput,
                         ForceInput( myinput.ncomponents,
                                     value_stash.data()+myinput.ncomponents*task_index,
                                     nderivatives_per_task,
                                     derivatives.data()),
                         forces
                       );
      }
      #pragma omp critical
      T::gatherThreads( forces );
    }
    // MPI Gather everything (this may be extended to the gpu thing, after makning it mpi-aware)
    if( !action->runInSerial() ) comm.Sum( forcesForApply );
  }

}

} // namespace PLMD
#endif
