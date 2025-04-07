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
/// The number of scalars we are calculating for each task
  unsigned nscalars{0};
//// The number of indexes that you use per task
  unsigned nindices_per_task{0};
/// This holds all the input data that is required to calculate all values for all tasks
  unsigned dataSize{0};
  double *inputdata{nullptr};
/// Bookeeping stuff for arguments
  std::size_t nargs{0};
  std::vector<std::size_t> ranks;
  std::vector<std::size_t> shapestarts;
  std::vector<std::size_t> shapedata;
  std::vector<std::size_t> ncols;
  std::vector<std::size_t> bookstarts;
  std::vector<std::size_t> booksizes;
  std::vector<std::size_t> bookeeping;
  std::vector<std::size_t> argstarts;
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
  /// Setup the arguments
  void setupArguments( const ActionWithArguments* action );
};

inline
void ParallelActionsInput::setupArguments( const ActionWithArguments* action ) {
  nargs = action->getNumberOfArguments();
  ranks.resize( nargs );
  shapestarts.resize( nargs );
  argstarts.resize( nargs );
  std::size_t s = 0, ts = 0;
  for(unsigned i=0; i<nargs; ++i) {
    Value* myarg = action->getPntrToArgument(i);
    shapestarts[i] = ts;
    ranks[i] = myarg->getRank();
    ts += ranks[i];
    argstarts[i] = s;
    s += myarg->getNumberOfStoredValues();
  }
  shapedata.resize( ts );
  ts = 0;
  ncols.resize( nargs );
  bookstarts.resize( nargs );
  booksizes.resize( nargs );
  std::size_t nbook = 0;
  for(unsigned i=0; i<nargs; ++i) {
    Value* myarg = action->getPntrToArgument(i);
    for(unsigned j=0; j<ranks[i]; ++j) {
      shapedata[ts] = myarg->getShape()[j];
      ++ts;
    }
    bookstarts[i] = nbook;
    if( ranks[i]==1 ) {
      ncols[i] = 1;
      booksizes[i] = 2*myarg->getShape()[0];
    } else if( ranks[i]==2 ) {
      ncols[i] = myarg->getNumberOfColumns();
      booksizes[i] = myarg->matrix_bookeeping.size();
    }
    nbook += booksizes[i];
  }
  bookeeping.resize( nbook );
  ts = 0;
  for(unsigned i=0; i<nargs; ++i) {
    Value* myarg = action->getPntrToArgument(i);
    if( ranks[i]==1 ) {
      for(unsigned j=0; j<myarg->getShape()[0]; ++j) {
        bookeeping[ts] = 1;
        bookeeping[ts+1] = 0;
        ts += 2;
      }
    } else if( ranks[i]==2 ) {
      for(unsigned j=0; j<myarg->matrix_bookeeping.size(); ++j) {
        bookeeping[ts] = myarg->matrix_bookeeping[j];
        ++ts;
      }
    }
  }
}

class ArgumentBookeepingHolder {
public:
  std::size_t rank;
  std::size_t ncols;
  std::size_t start;
  View<const std::size_t,helpers::dynamic_extent> shape;
  View<const std::size_t,helpers::dynamic_extent> bookeeping;
  ArgumentBookeepingHolder( std::size_t argno, const ParallelActionsInput& inp ) :
    rank(inp.ranks[argno]),
    ncols(inp.ncols[argno]),
    start(inp.argstarts[argno]),
    shape(inp.shapedata.data() + inp.shapestarts[argno], rank ),
    bookeeping(inp.bookeeping.data() + inp.bookstarts[argno], inp.booksizes[argno] ) {
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
  std::size_t nvals;
public:
  View<double> force;
  View2D<double> deriv;
  ForceInput( std::size_t nc, std::size_t nv, double* f, std::size_t nd, double* d ) :
    ncomp(nc),
    nvals(nv),
    force(f,nv),
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

template <class T>
class ParallelTaskManager {
public:
  using input_type= typename T::input_type;
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
  input_type actiondata;
//// This is used internally to get the number of elements in the value stash
  std::size_t getValueStashSize() const ;
/// This is used internally to gather the forces on the threads
  void gatherThreads( ForceOutput forces );
public:
  static void registerKeywords( Keywords& keys );
  ParallelTaskManager(ActionWithVector* av);
/// Setup the parallel task manager the three arguments are
/// nind = the number of indices that are used per task
/// nder = number of derivatives per component that is being calculated
/// nt = number of derivatives that need to gathered over the threads (default of zero means all derivatives need to be calculated in a way that is thread safe)
/// s = stride for values - default is equal to number of components
  void setupParallelTaskManager( std::size_t nind, std::size_t nder, int nt=-1, int s=-1 );
/// Copy the data from the underlying colvar into this parallel action
  void setActionInput( const input_type& adata );
/// Get the action input so we can use it
  input_type& getActionInput();
/// This runs all the tasks
  void runAllTasks();
/// Apply the forces on the parallel object
  void applyForces( std::vector<double>& forcesForApply );
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
  nthreaded_forces(0),
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
std::size_t ParallelTaskManager<T>::getValueStashSize() const {
  std::size_t valuesize=0;
  for(unsigned i=0; i<action->getNumberOfComponents(); ++i) {
    valuesize += (action->getConstPntrToComponent(i))->getNumberOfStoredValues();
  }
  return valuesize;
}


template <class T>
void ParallelTaskManager<T>::setupParallelTaskManager(
  std::size_t nind,
  std::size_t nder,
  int nt,
  int s ) {
  plumed_massert( action->getNumberOfComponents()>0, "there should be some components wen you setup the index list" );
  myinput.ncomponents = action->getNumberOfComponents();
  myinput.nscalars = myinput.ncomponents;
  if( s>0 ) {
    myinput.nscalars = s;
  }
  nderivatives_per_component = nder;
  value_stash.resize( getValueStashSize() );
  myinput.nindices_per_task = nind;

  if( nt<0 ) {
    nthreaded_forces = action->getNumberOfDerivatives();
  } else {
    nthreaded_forces = nt;
  }
  unsigned t=OpenMP::getNumThreads();
  if( useacc ) {
    t = 1;
  }
  omp_forces.resize(t);
  for(unsigned i=0; i<t; ++i) {
    omp_forces[i].resize(nthreaded_forces);
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
  // Transfer all the bookeeping information about the arguments
  myinput.setupArguments( action );
  // Reset the values at the start of the task loop
  std::size_t totalvals=getValueStashSize();
  if( value_stash.size()!=totalvals ) {
    value_stash.resize(totalvals);
  }
  value_stash.assign( value_stash.size(), 0.0 );
  if( useacc ) {
#ifdef __PLUMED_USE_OPENACC
    //I have a few problem with "this" <- meaning "this" pointer-  being copyed,
    // so I workarounded it with few copies
    ParallelActionsInput input = myinput;
    auto myinput_acc = OpenACC::fromToDataHelper(input);
    D t_actiondata = actiondata;
    auto actiondata_acc = OpenACC::fromToDataHelper(t_actiondata);

    auto value_stash_data = value_stash.data();
    auto partialTaskList_data = partialTaskList.data();
    auto vs_size = value_stash.size();

    const auto nderivPerComponent = nderivatives_per_component;
    const auto ndev_per_task = input.ncomponents*nderivPerComponent;
    //To future me/you:
    // I need to allocate this on the host to create a bigger temporay data array
    // on the device
    // by trying with double* x=nullptr, you will get a failure
    // another solution is acc_malloc and then device_ptr in the pragma
    // (but you have to remember the acc_free)
    std::vector<double> derivative(1);
    double * derivatives = derivative.data();

#pragma acc parallel loop present(input, t_actiondata) \
                          copyin(nactive_tasks, \
                                 ndev_per_task, \
                                 nderivPerComponent, \
                                 partialTaskList_data[0:nactive_tasks])\
                          copy(value_stash_data[0:vs_size]) \
                          create(derivatives[0:ndev_per_task*nactive_tasks]) \
                          default(none)
    for(unsigned i=0; i<nactive_tasks; ++i) {
      std::size_t task_index = partialTaskList_data[i];
      std::size_t val_pos = task_index*input.nscalars;
      ParallelActionsOutput myout {input.nscalars,
                                   value_stash_data+val_pos,
                                   ndev_per_task,
                                   derivatives+ndev_per_task*i};
      // Calculate the stuff in the loop for this action
      T::performTask( task_index, t_actiondata, input, myout );
    }
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
        std::size_t val_pos = task_index*myinput.nscalars;
        ParallelActionsOutput myout( myinput.nscalars,
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
  action->transferStashToValues( value_stash );
}

template <class T>
void ParallelTaskManager<T>::applyForces( std::vector<double>& forcesForApply ) {
  // Get the list of active tasks
  std::vector<unsigned> & partialTaskList( action->getListOfActiveTasks( action ) );
  unsigned nactive_tasks=partialTaskList.size();
  // Clear force buffer
  forcesForApply.assign( forcesForApply.size(), 0.0 );

  // Get all the input data so we can broadcast it to the GPU
  myinput.noderiv = false;
  // Retrieve the forces from the values
  action->transferForcesToStash( value_stash );

  if( useacc ) {
#ifdef __PLUMED_USE_OPENACC
    omp_forces[0].assign( omp_forces[0].size(), 0.0 );
    ParallelActionsInput input = myinput;
    auto myinput_acc = OpenACC::fromToDataHelper(input);
    D t_actiondata = actiondata;
    auto actiondata_acc = OpenACC::fromToDataHelper(t_actiondata);

    //passing raw pointer makes things easier in openacc
    auto value_stash_data = value_stash.data();
    const auto value_stash_size = value_stash.size();
    auto partialTaskList_data = partialTaskList.data();
    auto forcesForApply_data = forcesForApply.data();
    const auto forcesForApply_size = forcesForApply.size();
    auto omp_forces_data = omp_forces[0].data();
    const auto omp_forces_size = omp_forces[0].size();

    const auto nderivPerComponent = nderivatives_per_component;
    const auto ndev_per_task = input.ncomponents*nderivPerComponent;

    //To future me/you:
    // I need to allocate this on the host to create a bigger temporay data array
    // on the device
    // by trying with double* x=nullptr, you will get a failure
    // another solution is acc_malloc and then device_ptr in the pragma
    // (but you have to remember the acc_free)
    std::vector<double> derivative(1);
    double * derivatives = derivative.data();


#pragma acc parallel loop reduction(+:omp_forces_data[0:omp_forces_size])\
                            present(input,t_actiondata) \
                             copyin(nactive_tasks, \
                                    ndev_per_task, \
                                    nderivPerComponent ,\
                                    partialTaskList_data[0:nactive_tasks], \
                                    value_stash_data[0:value_stash_size]) \
                               copy(omp_forces_data[0:omp_forces_size], \
                                    forcesForApply_data[0:forcesForApply_size]) \
                             create(derivatives[0:ndev_per_task*nactive_tasks]) \
                            default(none)
    for(unsigned i=0; i<nactive_tasks; ++i) {
      //This may be changed to a shared array
      std::vector<double> valstmp( input.nscalars );
      std::size_t task_index = partialTaskList_data[i];
      ParallelActionsOutput myout( input.nscalars,
                                   valstmp.data(),
                                   nderivPerComponent,
                                   derivatives+ndev_per_task*i
                                 );

      // Calculate the stuff in the loop for this action
      T::performTask( task_index, t_actiondata, input, myout );
      // Gather the forces from the values
      T::gatherForces( task_index, t_actiondata, input,
                       ForceInput( input.ncomponents,
                                   input.nscalars,
                                   value_stash_data+input.nscalars*task_index,
                                   nderivPerComponent,
                                   derivatives+ndev_per_task*i),
                       ForceOutput { omp_forces_data,omp_forces_size, forcesForApply_data,forcesForApply_size }
                     );
    }
    gatherThreads({ omp_forces_data,omp_forces_size, forcesForApply_data,forcesForApply_size });
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
      omp_forces[t].assign( nthreaded_forces, 0.0 );
      ForceOutput forces{ omp_forces[t], forcesForApply };
      std::vector<double> fake_vals( myinput.nscalars );
      std::vector<double> derivatives( myinput.ncomponents*nderivatives_per_component );
      #pragma omp for nowait
      for(unsigned i=rank; i<nactive_tasks; i+=stride) {
        std::size_t task_index = partialTaskList[i];
        ParallelActionsOutput myout( myinput.nscalars,
                                     fake_vals.data(),
                                     derivatives.size(),
                                     derivatives.data() );
        // Calculate the stuff in the loop for this action
        T::performTask( task_index, actiondata, myinput, myout );

        // Gather the forces from the values
        T::gatherForces( task_index,
                         actiondata,
                         myinput,
                         ForceInput( myinput.ncomponents,
                                     myinput.nscalars,
                                     value_stash.data()+myinput.nscalars*task_index,
                                     nderivatives_per_component,
                                     derivatives.data()),
                         forces );
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
void ParallelTaskManager<T>::gatherThreads( ForceOutput forces ) {
  //Forceoutput is basically two spans, so it is ok to pass it by value
  unsigned k=0;
  for(unsigned n=forces.thread_unsafe.size()-nthreaded_forces; n<forces.thread_unsafe.size(); ++n) {
    forces.thread_unsafe[n] += forces.thread_safe[k];
    ++k;
  }
}

} // namespace PLMD
#endif
