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

#include "tools/ColvarOutput.h"
#include "tools/OpenACC.h"

namespace PLMD {

struct ParallelActionsInput {
/// Do we need to calculate the derivatives
  bool noderiv{false};
/// Periodic boundary conditions
  const Pbc* pbc;
/// The number of components the underlying action is computing
  unsigned ncomponents{0};
//// The number of indexes that you use per task
  unsigned nindices_per_task{0};
/// The end of the derivatives that are only affected by one task
  unsigned threadsafe_derivatives_end{0};
/// The start of the thread unsafe forces
  unsigned threadunsafe_forces_start{0};
/// This holds all the input data that is required to calculate all values for all tasks
  unsigned dataSize{0};
  double *inputdata{nullptr};
  // std::vector<double> inputdata{};
/// Default constructor
  ParallelActionsInput( const Pbc& box )
    : pbc(&box) {}
  /// the copy is needed due to some openACC resistance to pass the inputdata as a std::vector (and some "this" problems, please do not ask, I'll get grumpy)
  ParallelActionsInput( const ParallelActionsInput& other )
    : noderiv {other.noderiv},
      //just copies the ptr
      pbc  {other.pbc},
      ncomponents {other.ncomponents},
      nindices_per_task {other.nindices_per_task},
      threadsafe_derivatives_end {other.threadsafe_derivatives_end},
      threadunsafe_forces_start {other.threadunsafe_forces_start},
      dataSize {other.dataSize},
      //just copies the value of the pointer
      inputdata  {other.inputdata}
  {}

  ParallelActionsInput& operator=( const ParallelActionsInput& other )  = delete;

  //helper function to bring data to the device in a controlled way
  void toACCDevice()const {
#pragma acc enter data copyin(this[0:1], noderiv, pbc[0:1],ncomponents, \
  nindices_per_task, threadsafe_derivatives_end, threadunsafe_forces_start,\
   dataSize, inputdata[0:dataSize])

    pbc->toACCDevice();
  }
  //helper function to remove data from the device in a controlled way
  void removeFromACCDevice() const  {
    pbc->removeFromACCDevice();
    // assuming dataSize is not changed
#pragma acc exit data delete(inputdata[0:dataSize],dataSize,nindices_per_task, \
  ncomponents, threadsafe_derivatives_end, threadunsafe_forces_start, \
  pbc[0:1],noderiv,this[0:1])
  }
};

struct ParallelActionsOutput {
  View<double> values;
  View<double> derivatives;
  ParallelActionsOutput( std::size_t ncomp, double* v, std::size_t ndev, double* d )
    : values(v,ncomp),
      derivatives(d,ndev) {}
};

class ForceInput {
private:
  std::size_t ncomp;
public:
  View<double> force;
  View2D<double> deriv;
  ForceInput( std::size_t nc, double* f, std::size_t nd, double* d ) :
    ncomp(nc),
    force(f,nc),
    deriv(d,nc,nd) {}
};

//There is no need to pass this as reference:
struct ForceOutput {
  //I would suggest to invert the name or to be clearer
// like something that recalls that "thread_safe" will be need to be reducted (hence it is NOT thread safe)
  View<double> thread_safe;
  //these are the forces that we promise will not provoke races
  View<double> thread_unsafe;
  //const T* is a ptr to const T
  //T* const is a conts ptr to a modifiable T
  ForceOutput(std::vector<double>& reduced, std::vector<double>& notReduced) :
    thread_safe{reduced.data(),
                reduced.size()},
    thread_unsafe{notReduced.data()
                  ,notReduced.size()}
  {}
  ForceOutput(double* reduced, size_t rs, double* notReduce, size_t nrsz):
    thread_safe{reduced,rs},
    thread_unsafe{notReduce,nrsz}
  {}
  ForceOutput(const ForceOutput &) = default;
  ForceOutput(ForceOutput &&) = default;
  ForceOutput& operator=(const ForceOutput &) = delete;
  ForceOutput& operator=(ForceOutput &&) = delete;
};

namespace PTMUtils {
template<class, class = void>
constexpr bool has_gatherForces_custom = false;

//this verifies that T has a method gatherForces_custom that can be called with this signature
template<class T>
constexpr bool has_gatherForces_custom <
T,
std::void_t<
decltype(T::gatherForces_custom(
           std::declval<unsigned >(),
           std::declval<size_t >(),
           std::declval<size_t >(),
           std::declval<const typename T::input_type & >(),
           std::declval<const ParallelActionsInput& >(),
           std::declval<View<unsigned> >(),
           std::declval<double *>(),
           std::declval<double *>(),
           std::declval<View<double> >()
         ))
>
> = true;
} //namespace PTMUtils

template <class T>
class ParallelTaskManager {
public:
  using input_type= typename T::input_type;
  static constexpr bool has_custom_gather=PTMUtils::has_gatherForces_custom<T>;
private:
/// The underlying action for which we are managing parallel tasks
  ActionWithVector* action;
/// The MPI communicator
  Communicator& comm;
/// Is this an action with matrix
  bool ismatrix;
/// Are we using acc for parallisation
  bool useacc;
/// Number of derivatives for quantity being calculated
  std::size_t nderivatives_per_component;
/// This holds the values before we pass them to the value
  std::vector<double> value_stash;
/// A tempory set of vectors for holding forces over threads
  std::vector<std::vector<double> > omp_forces;
/// This structs is used to pass data between the parallel interface and the function caller
  ParallelActionsInput myinput;
//this holds the data for myinput that will be passed though myinput
  std::vector<double> input_buffer;
/// This holds data for that the underlying action needs to do the calculation
  input_type actiondata;
/// This is used internally to gather the forces on the threads
  void gatherThreads( ForceOutput forces );
public:
  static void registerKeywords( Keywords& keys );
  ParallelTaskManager(ActionWithVector* av);
/// Setup the parallel task manager the three arguments are
/// nind = the number of indices that are used per task
/// nder = number of derivatives per component that is being calculated
/// nder_tus = number of derivatives that can be directly added to the applyForces array
/// nforce_ts = size of part of forces array that need to be handled with thread safe procedure
  void setupParallelTaskManager( std::size_t nind, std::size_t nder, std::size_t nder_tus, std::size_t nforce_ts );
/// Copy the data from the underlying colvar into this parallel action
  void setActionInput( const input_type& adata );
/// Get the action input so we can use it
  input_type& getActionInput();
/// This runs all the tasks
  void runAllTasks();
/// Apply the forces on the parallel object
  void applyForces( std::vector<double>& forcesForApply );
/// This is used to gather forces that are thread safe
  static void gatherThreadSafeForces( const ParallelActionsInput& input,
                                      View<std::size_t,helpers::dynamic_extent> force_indices,
                                      const ForceInput& fdata,
                                      View<double,helpers::dynamic_extent> forces );
/// This is used to gather forces that are not thread safe
  static void gatherThreadUnsafeForces( const ParallelActionsInput& input,
                                        View<std::size_t,helpers::dynamic_extent> force_indices,
                                        const ForceInput& fdata,
                                        View<double,helpers::dynamic_extent> forces );
};

template <class T>
void ParallelTaskManager<T>::registerKeywords( Keywords& keys ) {
  keys.addFlag("USEGPU",false,"run this calculation on the GPU");
}

template <class T>
ParallelTaskManager<T>::ParallelTaskManager(ActionWithVector* av):
  action(av),
  comm(av->comm),
  ismatrix(false),
  useacc(false),
  nderivatives_per_component(0),
  myinput(av->getPbc()) {
  ActionWithMatrix* am=dynamic_cast<ActionWithMatrix*>(av);
  if(am) {
    ismatrix=true;
  }
  action->parseFlag("USEGPU",useacc);
#ifdef __PLUMED_USE_OPENACC
  if( useacc ) {
    action->log.printf("  using GPU to calculate this action\n");
  }
#else
  if( useacc ) {
    action->error("cannot use USEGPU flag as PLUMED has not been compiled with openacc");
  }
#endif
}

template <class T>
void ParallelTaskManager<T>::setupParallelTaskManager(
  std::size_t nind,
  std::size_t nder,
  std::size_t nder_tus,
  std::size_t nforce_ts ) {
  //int nt ) {
  plumed_massert( action->getNumberOfComponents()>0, "there should be some components wen you setup the index list" );
  std::size_t valuesize=(action->getConstPntrToComponent(0))->getNumberOfStoredValues();
  for(unsigned i=1; i<action->getNumberOfComponents(); ++i) {
    plumed_assert( valuesize==(action->getConstPntrToComponent(i))->getNumberOfStoredValues() );
  }
  myinput.ncomponents = action->getNumberOfComponents();
  nderivatives_per_component = nder;
  value_stash.resize( valuesize*action->getNumberOfComponents() );
  myinput.nindices_per_task = nind;
  plumed_massert( nder_tus<nder, "number of derivatives that can be handled in thread safe manner must be less than total number of derivatives");
  myinput.threadsafe_derivatives_end = nder_tus;
  myinput.threadunsafe_forces_start = action->getNumberOfDerivatives() - nforce_ts;
  unsigned t=OpenMP::getNumThreads();
  if( useacc ) {
    t = 1;
  }
  omp_forces.resize(t);
  for(unsigned i=0; i<t; ++i) {
    omp_forces[i].resize(nforce_ts);
  }
}

template <class T>
void ParallelTaskManager<T>::setActionInput( const input_type& adata ) {
  actiondata=adata;
}

template <class T>
typename ParallelTaskManager<T>::input_type& ParallelTaskManager<T>::getActionInput() {
  return actiondata;
}

template <class T>
void ParallelTaskManager<T>::runAllTasks() {
  // Get the list of active tasks
  std::vector<unsigned> & partialTaskList( action->getListOfActiveTasks( action ) );
  unsigned nactive_tasks=partialTaskList.size();
  // Get all the input data so we can broadcast it to the GPU
  myinput.noderiv = true;
  action->getInputData( input_buffer );
  myinput.dataSize = input_buffer.size();
  myinput.inputdata = input_buffer.data();
  std::fill (value_stash.begin(),value_stash.end(), 0.0);
  if( useacc ) {
#ifdef __PLUMED_USE_OPENACC
    //I have a few problem with "this" <- meaning "this" pointer-  being copyed,
    // so I workarounded it with few copies
    ParallelActionsInput input = myinput;
    auto myinput_acc = OpenACC::fromToDataHelper(input);
    input_type t_actiondata = actiondata;
    auto actiondata_acc = OpenACC::fromToDataHelper(t_actiondata);

    //template type is deduced
    OpenACC::memoryManager vs{value_stash};
    auto value_stash_data = vs.devicePtr();

    OpenACC::memoryManager ptl{partialTaskList};
    auto partialTaskList_data = ptl.devicePtr();

    const auto nderivPerComponent = nderivatives_per_component;
    const auto ndev_per_task = input.ncomponents*nderivPerComponent;

    OpenACC::memoryManager<double>dev(ndev_per_task*nactive_tasks);
    auto derivatives = dev.devicePtr();

#pragma acc parallel loop present(input, t_actiondata) \
                           copyin(nactive_tasks, \
                                 ndev_per_task, \
                                 nderivPerComponent)\
                        deviceptr(derivatives, \
                                  partialTaskList_data, \
                                  value_stash_data) \
                          default(none)
    for(unsigned i=0; i<nactive_tasks; ++i) {
      std::size_t task_index = partialTaskList_data[i];
      std::size_t val_pos = task_index*input.ncomponents;
      ParallelActionsOutput myout {input.ncomponents,
                                   value_stash_data+val_pos,
                                   ndev_per_task,
                                   derivatives+ndev_per_task*i};
      // Calculate the stuff in the loop for this action
      T::performTask( task_index, t_actiondata, input, myout );
    }
    vs.copyFromDevice(value_stash.data());
#else
    plumed_merror("cannot use USEGPU flag if PLUMED has not been compiled with openACC");
#endif
  } else {
    // Get the MPI details
    unsigned stride=comm.Get_size();
    unsigned rank=comm.Get_rank();
    if(action->runInSerial()) {
      stride=1;
      rank=0;
    }

    // Get number of threads for OpenMP
    unsigned nt=OpenMP::getNumThreads();
    if( nt*stride*10>nactive_tasks ) {
      nt=nactive_tasks/stride/10;
    }
    if( nt==0 ) {
      nt=1;
    }

    #pragma omp parallel num_threads(nt)
    {
      std::vector<double> derivatives( myinput.ncomponents*nderivatives_per_component );
      #pragma omp for nowait
      for(unsigned i=rank; i<nactive_tasks; i+=stride) {
        std::size_t task_index = partialTaskList[i];
        std::size_t val_pos = task_index*myinput.ncomponents;
        ParallelActionsOutput myout( myinput.ncomponents,
                                     value_stash.data()+val_pos,
                                     nderivatives_per_component,
                                     derivatives.data() );
        // Calculate the stuff in the loop for this action
        T::performTask( task_index, actiondata, myinput, myout );
      }
    }
    // MPI Gather everything
    if( !action->runInSerial() ) {
      comm.Sum( value_stash );
    }
  }

  // And transfer the value to the output values
  for(unsigned i=0; i<action->getNumberOfComponents(); ++i) {
    Value* myval = action->copyOutput(i);
    for(unsigned j=0; j<myval->getNumberOfStoredValues(); ++j) {
      myval->set( j, value_stash[j*myinput.ncomponents+i] );
    }
  }
}

template <class T>
void ParallelTaskManager<T>::applyForces( std::vector<double>& forcesForApply ) {
  // Get the list of active tasks
  std::vector<unsigned> & partialTaskList( action->getListOfActiveTasks( action ) );
  unsigned nactive_tasks=partialTaskList.size();
  // Clear force buffer
  forcesForApply.assign( forcesForApply.size(), 0.0 );
  //TODO: check if std::fill is faster (i get conflicting answers on the net)
  //std::fill (forcesForApply.begin(),forcesForApply.end(), 0.0);
  // Get all the input data so we can broadcast it to the GPU
  myinput.noderiv = false;
  // Retrieve the forces from the values
  for(unsigned i=0; i<action->getNumberOfComponents(); ++i) {
    auto myval = action->getConstPntrToComponent(i);
    for(unsigned j=0; j<myval->getNumberOfStoredValues(); ++j) {
      value_stash[j*myinput.ncomponents+i] = myval->getForce( j );
    }
  }

  if( useacc ) {
#ifdef __PLUMED_USE_OPENACC
    omp_forces[0].assign( omp_forces[0].size(), 0.0 );
    ParallelActionsInput input = myinput;
    auto myinput_acc = OpenACC::fromToDataHelper(input);
    input_type t_actiondata = actiondata;
    auto actiondata_acc = OpenACC::fromToDataHelper(t_actiondata);

    //template type is deduced
    OpenACC::memoryManager vs{value_stash};
    auto value_stash_data = vs.devicePtr();

    OpenACC::memoryManager ptl{partialTaskList};
    auto partialTaskList_data = ptl.devicePtr();

    OpenACC::memoryManager ffa {forcesForApply};
    auto forcesForApply_data = ffa.devicePtr();
    const auto forcesForApply_size = ffa.size();

//    OpenACC::memoryManager ofd {omp_forces[0]};
//    auto omp_forces_data = ofd.devicePtr();
//    const auto omp_forces_size = ofd.size();
//    const auto reduction_size = (has_custom_gather)
//                                ? 0:omp_forces_size;
    const auto nderivPerComponent = nderivatives_per_component;
    const auto ndev_per_task = input.ncomponents*nderivPerComponent;

    OpenACC::memoryManager<double> dev{ndev_per_task*nactive_tasks};
    auto derivatives = dev.devicePtr();
    OpenACC::memoryManager<std::size_t> ind{ndev_per_task*nactive_tasks};
    auto indices = ind.devicePtr();
    OpenACC::memoryManager<double> vtmp{input.ncomponents*nactive_tasks};
    auto valstmp = vtmp.devicePtr();

#pragma acc data present(input,t_actiondata) \
                  copyin(nactive_tasks, \
                         nderivPerComponent, \
                         forcesForApply_size, \
                         ndev_per_task) \
               deviceptr(derivatives, \
                         indices, \
                         value_stash_data, \
                         partialTaskList_data, \
                         forcesForApply_data, \
                         valstmp) \
                 default(none)
    {
#pragma acc parallel loop
      for(unsigned i=0; i<nactive_tasks; ++i) {
        std::size_t task_index = partialTaskList_data[i];
        ParallelActionsOutput myout { input.ncomponents,
                                      valstmp+input.ncomponents*i,
                                      ndev_per_task,
                                      derivatives+ndev_per_task*i};
        // Calculate the stuff in the loop for this action
        T::performTask( task_index, t_actiondata, input, myout );
        // Get the indices that are used here
        View<std::size_t> force_indices( indices + ndev_per_task*i,
                                         ndev_per_task );
        T::getForceIndices( task_index,
                            forcesForApply_size,
                            t_actiondata,
                            input,
                            force_indices );

        ForceInput finput( input.ncomponents,
                           value_stash_data+input.ncomponents*task_index,
                           nderivPerComponent,
                           derivatives+ndev_per_task*i);

        // Gather forces that can be gathered locally
        gatherThreadSafeForces( input,
                                force_indices,
                                finput,
                                View<double>(forcesForApply_data,forcesForApply_size));

      }

#pragma acc parallel loop
      for(unsigned v=input.threadunsafe_forces_start; v<forcesForApply_size; ++v) {
        double tmp = 0.0;
#pragma acc loop reduction(+:tmp)
        for(unsigned t=0; t<nactive_tasks; ++t) {
          int m=-1;
          std::size_t task_index = partialTaskList_data[t];
          View<const std::size_t> force_indices( indices + ndev_per_task*t,
                                                 ndev_per_task );
// #pragma acc loop seq
          for(unsigned d=input.threadsafe_derivatives_end; d<ndev_per_task; ++d) {
            if( force_indices[d]==v ) {
              m=d;
              break;
            }
          }
          if( m<0 ) {
            continue ;
          }
          auto fdata = ForceInput { input.ncomponents,
                                    value_stash_data+input.ncomponents*task_index,
                                    nderivPerComponent,
                                    derivatives+ndev_per_task*t};
          for(unsigned i=0; i<input.ncomponents; ++i) {
            tmp += fdata.force[i]*fdata.deriv[i][m];
          }
        }
        forcesForApply_data[v] = tmp;
      }
    }

    ffa.copyFromDevice(forcesForApply.data());
//    ofd.copyFromDevice(omp_forces[0].data());
    // for(unsigned v=0; v<9; ++v) {
    //   std::cout << omp_forces[0][v] << " ";
    // }
    // std::cout << "\n";
//    gatherThreads(ForceOutput{ omp_forces[0], forcesForApply });
#else
    plumed_merror("cannot use USEGPU flag if PLUMED has not been compiled with openACC");
#endif
  } else {
    // Get the MPI details
    unsigned stride=comm.Get_size();
    unsigned rank=comm.Get_rank();
    if(action->runInSerial()) {
      stride=1;
      rank=0;
    }

    // Get number of threads for OpenMP
    unsigned nt=OpenMP::getNumThreads();
    if( nt*stride*10>nactive_tasks ) {
      nt=nactive_tasks/stride/10;
    }
    if( nt==0 ) {
      nt=1;
    }

    #pragma omp parallel num_threads(nt)
    {
      const unsigned t=OpenMP::getThreadNum();
      omp_forces[t].assign( forcesForApply.size()-myinput.threadunsafe_forces_start, 0.0 );
      ForceOutput forces{ omp_forces[t], forcesForApply };
      std::vector<double> fake_vals( myinput.ncomponents );
      std::vector<double> derivatives( myinput.ncomponents*nderivatives_per_component );
      std::vector<std::size_t> indices( nderivatives_per_component );
      View<std::size_t,helpers::dynamic_extent> force_indices( indices.data(), nderivatives_per_component );
      #pragma omp for nowait
      for(unsigned i=rank; i<nactive_tasks; i+=stride) {
        std::size_t task_index = partialTaskList[i];
        ParallelActionsOutput myout( myinput.ncomponents,
                                     fake_vals.data(),
                                     derivatives.size(),
                                     derivatives.data() );
        // Calculate the stuff in the loop for this action
        T::performTask( task_index, actiondata, myinput, myout );

        // If this is a matrix this returns a number that isn't one as we have to loop over the columns
        std::size_t nfpt = T::getNumberOfValuesPerTask( task_index, actiondata );
        for(unsigned j=0; j<nfpt; ++j) {
            // Get the force indices
            T::getForceIndices( task_index, 
                                j,
                                forcesForApply.size(),
                                actiondata,
                                myinput,
                                force_indices );

            // Create a force input object
            ForceInput finput( myinput.ncomponents,
                               value_stash.data()+myinput.ncomponents*task_index,
                               nderivatives_per_component,
                               derivatives.data());
 
            // Gather forces that are thread safe
            gatherThreadSafeForces( myinput,
                                    force_indices,
                                    finput,
                                    View<double,helpers::dynamic_extent>(forcesForApply.data(),forcesForApply.size()) );

            // Gather forces that are not thread safe
            gatherThreadUnsafeForces( myinput,
                                      force_indices,
                                      finput,
                                      View<double,helpers::dynamic_extent>(omp_forces[t].data(),omp_forces[t].size()) );
        }
      }

      #pragma omp critical
      gatherThreads( forces );
    }
    // MPI Gather everything (this must be extended to the gpu thing, after makning it mpi-aware)
    if( !action->runInSerial() ) {
      comm.Sum( forcesForApply );
    }
  }

}

template <class T>
void ParallelTaskManager<T>::gatherThreadSafeForces( const ParallelActionsInput& input,
    const View<std::size_t,helpers::dynamic_extent> force_indices,
    const ForceInput& fdata,
    View<double,helpers::dynamic_extent> forces ) {
  for(unsigned i=0; i<input.ncomponents; ++i) {
    double ff = fdata.force[i];
    for(unsigned j=0; j<input.threadsafe_derivatives_end; ++j) {
      forces[ force_indices[j] ] += ff*fdata.deriv[i][j];
    }
  }
}

template <class T>
void ParallelTaskManager<T>::gatherThreadUnsafeForces( const ParallelActionsInput& input,
    View<std::size_t,helpers::dynamic_extent> force_indices,
    const ForceInput& fdata,
    View<double,helpers::dynamic_extent> forces ) {
  for(unsigned i=0; i<input.ncomponents; ++i) {
    double ff = fdata.force[i];
    for(unsigned j=input.threadsafe_derivatives_end; j<force_indices.size(); ++j) {
      forces[ force_indices[j] - input.threadunsafe_forces_start ] += ff*fdata.deriv[i][j];
    }
  }
}

template <class T>
void ParallelTaskManager<T>::gatherThreads( ForceOutput forces ) {
  //Forceoutput is basically two spans, so it is ok to pass it by value
  unsigned k=0;
  for(unsigned n=forces.thread_unsafe.size()-forces.thread_safe.size(); n<forces.thread_unsafe.size(); ++n) {
    forces.thread_unsafe[n] += forces.thread_safe[k];
    ++k;
  }
}

} // namespace PLMD
#endif
