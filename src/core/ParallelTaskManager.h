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

struct ArgumentsBookkeeping {
  std::size_t nargs{0};
  std::vector<std::size_t> ranks;
  std::vector<std::size_t> shapestarts;
  std::vector<std::size_t> shapedata;
  std::vector<std::size_t> ncols;
  std::vector<std::size_t> bookstarts;
  std::vector<std::size_t> booksizes;
  std::vector<std::size_t> bookeeping;
  std::vector<std::size_t> argstarts;
  void setupArguments( const ActionWithArguments* action );
};

inline
void ArgumentsBookkeeping::setupArguments( const ActionWithArguments* action ) {
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

struct ParallelActionsInput {
  /// Do we need to calculate the derivatives
  bool noderiv{false};
  /// Periodic boundary conditions
  const Pbc* pbc;
  /// The number of components the underlying action is computing
  unsigned ncomponents{0};
  /// The number of scalars we are calculating for each task
  unsigned nscalars{0};
  /// Number of derivatives for each scalar being calculated
  unsigned nderivatives_per_scalar{0};
  /// The start of the thread unsafe forces
  unsigned threadunsafe_forces_start{0};
  /// This holds all the input data that is required to calculate all values for all tasks
  unsigned dataSize{0};
  double *inputdata{nullptr};
  /// Bookeeping stuff for arguments
  std::size_t nargs{0};
  std::size_t ranks_size{0};
  const std::size_t* ranks{nullptr};
  std::size_t shapestarts_size{0};
  const std::size_t* shapestarts{nullptr};
  std::size_t shapedata_size{0};
  const std::size_t* shapedata{nullptr};
  std::size_t ncols_size{0};
  const std::size_t* ncols{nullptr};
  std::size_t bookstarts_size{0};
  const std::size_t* bookstarts{nullptr};
  std::size_t booksizes_size{0};
  const std::size_t* booksizes{nullptr};
  std::size_t bookeeping_size{0};
  const std::size_t* bookeeping{nullptr};
  std::size_t argstarts_size{0};
  const std::size_t* argstarts{nullptr};
  /// Default constructor
  ParallelActionsInput( const Pbc& box )
    : pbc(&box) {}
  /// the copy is needed due to some openACC resistance to pass the inputdata as a std::vector (and some "this" problems, please do not ask, I'll get grumpy)
  ParallelActionsInput( const ParallelActionsInput& other ) = default;
  ParallelActionsInput(ParallelActionsInput&& other ) = default;

  ParallelActionsInput& operator=( const ParallelActionsInput& other )  = delete;

  //helper function to bring data to the device in a controlled way
  void toACCDevice()const {
#pragma acc enter data copyin(this[0:1], noderiv, pbc[0:1],ncomponents, \
    nscalars, nderivatives_per_scalar, threadunsafe_forces_start, \
     dataSize, inputdata[0:dataSize])
    if (nargs>0) {
#pragma acc enter data copyin( nargs, \
     ranks_size, ranks[0:ranks_size], \
     shapestarts_size, shapestarts[0:shapestarts_size], \
     shapedata_size, shapedata[0:shapedata_size], \
     ncols_size, ncols[0:ncols_size], \
     bookstarts_size, bookstarts[0:bookstarts_size], \
     booksizes_size, booksizes[0:booksizes_size], \
     bookeeping_size, bookeeping[0:bookeeping_size], \
     argstarts_size, argstarts[0:argstarts_size] \
    )
    }
    pbc->toACCDevice();
  }
  //helper function to remove data from the device in a controlled way
  void removeFromACCDevice() const  {
    pbc->removeFromACCDevice();
    // assuming dataSize is not changed
    if (nargs>0) {
#pragma acc exit data delete( \
    shapestarts[0:shapestarts_size], shapestarts_size, \
    shapedata[0:shapedata_size], shapedata_size, \
    ncols[0:ncols_size], ncols_size, \
    bookstarts[0:bookstarts_size], bookstarts_size, \
    booksizes[0:booksizes_size], booksizes_size, \
    bookeeping[0:bookeeping_size], bookeeping_size, \
    argstarts[0:argstarts_size], argstarts_size, \
    ranks[0:ranks_size], ranks_size, \
    nargs )
    }

#pragma acc exit data delete( \
    inputdata[0:dataSize], dataSize, \
    threadunsafe_forces_start, nderivatives_per_scalar, \
    nscalars, ncomponents, pbc[0:1],noderiv,this[0:1])
  }
  /// Setup the arguments
  void setupArguments( const ArgumentsBookkeeping& ab );
};

inline void ParallelActionsInput::setupArguments( const ArgumentsBookkeeping& ab ) {
  nargs = ab.nargs;
  ranks = ab.ranks.data();
  ranks_size = ab.ranks.size();
  shapestarts = ab.shapestarts.data();
  shapestarts_size = ab.shapestarts.size();
  shapedata = ab.shapedata.data();
  shapedata_size = ab.shapedata.size();
  ncols = ab.ncols.data();
  ncols_size = ab.ncols.size();
  bookstarts = ab.bookstarts.data();
  bookstarts_size = ab.bookstarts.size();
  booksizes = ab.booksizes.data();
  booksizes_size = ab.booksizes.size();
  bookeeping = ab.bookeeping.data();
  bookeeping_size = ab.bookeeping.size();
  argstarts = ab.argstarts.data();
  argstarts_size = ab.argstarts.size();
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
    shape(inp.shapedata + inp.shapestarts[argno], rank ),
    bookeeping(inp.bookeeping + inp.bookstarts[argno], inp.booksizes[argno] ) {
  }
};

struct ParallelActionsOutput {
  View<double> values;
  View<double> derivatives;
  View<double> buffer;
  ParallelActionsOutput( std::size_t ncomp, double* v, std::size_t ndev, double* d, std::size_t nb, double* b )
    : values(v,ncomp),
      derivatives(d,ndev),
      buffer(b,nb) {}
};

class ForceIndexHolder {
public:
  View<std::size_t,helpers::dynamic_extent> threadsafe_derivatives_end;
  View<std::size_t,helpers::dynamic_extent> tot_indices;
  View2D<std::size_t,helpers::dynamic_extent,helpers::dynamic_extent> indices;
  ForceIndexHolder( std::size_t nc, std::size_t nd, std::size_t* ind ) :
    threadsafe_derivatives_end(ind,nc),
    tot_indices(ind+nc,nc),
    indices(ind+2*nc,nc,nd) {}
};

class ForceInput {
private:
  std::size_t nvals;
public:
  View<double> force;
  View2D<double> deriv;
  ForceInput( std::size_t nv, double* f, std::size_t nd, double* d ) :
    nvals(nv),
    force(f,nv),
    deriv(d,nv,nd) {}
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

// namespace PTMUtils {
// template<class, class = void>
// constexpr bool has_gatherForces_custom = false;
//
// //this verifies that T has a method gatherForces_custom that can be called with this signature
// template<class T>
// constexpr bool has_gatherForces_custom <
// T,
// std::void_t<
// decltype(T::gatherForces_custom(
//            std::declval<unsigned >(),
//            std::declval<size_t >(),
//            std::declval<size_t >(),
//            std::declval<const typename T::input_type & >(),
//            std::declval<const ParallelActionsInput& >(),
//            std::declval<View<unsigned> >(),
//            std::declval<double *>(),
//            std::declval<double *>(),
//            std::declval<View<double> >()
//          ))
// >
// > = true;
//
// template<class, class = void>
// constexpr bool has_gatherForces_GPU = false;
//
// //this verifies that T has a method gatherForces_custom that can be called with this signature
// template<class T>
// constexpr bool has_gatherForces_GPU <
// T,
// std::void_t<
// decltype(T::gatherForcesGPU(
//            std::declval<unsigned >(),
//            std::declval<const typename T::input_type & >(),
//            std::declval<const ParallelActionsInput& >(),
//            std::declval<const ForceInput& >(),
//            std::declval<ForceOutput >()
//          ))
// >
// > = true;
//
// /// If the template class has virialSize, otherwise is 9
// template<class, class=void>
// constexpr size_t virialSize = 9;
//
// template<class T>
// constexpr size_t virialSize<T, std::void_t<decltype(T::virialSize),
// //this ensures that T::virialSize is a static member
//           std::enable_if_t<!std::is_member_pointer_v<decltype(&T::virialSize)>>
//           >>
//           = T::virialSize;
// } //namespace PTMUtils

template <class T>
class ParallelTaskManager {
public:
  using input_type= typename T::input_type;
//  static constexpr bool has_custom_gather=PTMUtils::has_gatherForces_custom<T>;
//  static constexpr bool has_GPU_gather=PTMUtils::has_gatherForces_GPU<T>;
//  static constexpr size_t virialSize = PTMUtils::virialSize<T>;
private:
/// The underlying action for which we are managing parallel tasks
  ActionWithVector* action;
/// The MPI communicator
  Communicator& comm;
/// Is this an action with matrix
  bool ismatrix;
/// Are we using acc for parallisation
  bool useacc;
/// Number of derivatives calculated for each task
  std::size_t nderivatives_per_task;
/// The number of forces on each thread
/// The number of forces on each thread
  std::size_t nthreaded_forces;
/// This holds the values before we pass them to the value
  std::vector<double> value_stash;
/// A tempory set of vectors for holding forces over threads
  std::vector<std::vector<double> > omp_forces;
/// This structs is used to pass data between the parallel interface and the function caller
  ParallelActionsInput myinput;
  ArgumentsBookkeeping argumentsMap;
//this holds the data for myinput that will be passed though myinput
  std::vector<double> input_buffer;
/// This holds tempory data that we use in performTask
  std::size_t workspace_size;
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
/// nder = number of derivatives per scalar
/// nforce_ts = number of forces that are modified by multiple tasks
  void setupParallelTaskManager( std::size_t nder, std::size_t nforce_ts );
/// Copy the data from the underlying colvar into this parallel action
  void setActionInput( const input_type& adata );
/// Creating the size of the workspace
  void setWorkspaceSize( std::size_t size );
/// Get the action input so we can use it
  input_type& getActionInput();
/// This runs all the tasks
  void runAllTasks();
/// Apply the forces on the parallel object
  void applyForces( std::vector<double>& forcesForApply );
/// This is used to gather forces that are thread safe
  static void gatherThreadSafeForces( const ParallelActionsInput& input,
                                      ForceIndexHolder force_indices,
                                      const ForceInput& fdata,
                                      View<double,helpers::dynamic_extent> forces );
/// This is used to gather forces that are not thread safe
  static void gatherThreadUnsafeForces( const ParallelActionsInput& input,
                                        ForceIndexHolder force_indices,
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
  nderivatives_per_task(0),
  nthreaded_forces(0),
  myinput(av->getPbc()),
  workspace_size(0) {
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
void ParallelTaskManager<T>::setupParallelTaskManager( std::size_t nder,
    std::size_t nforce_ts ) {
  plumed_massert( action->getNumberOfComponents()>0, "there should be some components wen you setup the index list" );
  myinput.ncomponents = action->getNumberOfComponents();
  unsigned ntasks=0;
  action->getNumberOfTasks( ntasks );
  myinput.nscalars = 0;
  for(unsigned i=0; i<myinput.ncomponents; ++i) {
    if( (action->copyOutput(i))->getRank()==1 ) {
      myinput.nscalars += 1;
    } else if( (action->copyOutput(i))->getRank()==2 ) {
      if( ntasks==(action->copyOutput(i))->getShape()[0] ) {
        myinput.nscalars += (action->copyOutput(i))->getNumberOfColumns();
      } else {
        myinput.nscalars += 1;
      }
    }
  }
  myinput.nderivatives_per_scalar = nder;
  nderivatives_per_task = nder*myinput.nscalars;
  value_stash.resize( getValueStashSize() );
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
void ParallelTaskManager<T>::setWorkspaceSize( std::size_t size ) {
  workspace_size = size;
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
  argumentsMap.setupArguments( action );
  myinput.setupArguments( argumentsMap );
  // Reset the values at the start of the task loop
  std::size_t totalvals=getValueStashSize();
  if( value_stash.size()!=totalvals ) {
    value_stash.resize(totalvals);
  }
  std::fill (value_stash.begin(),value_stash.end(), 0.0);
  if( useacc ) {
#ifdef __PLUMED_USE_OPENACC
    if (comm.Get_rank()== 0) {// no multigpu shenanigans until this works
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

      const auto ndev_per_task = nderivatives_per_task;

      OpenACC::memoryManager<double> buff{};
      if ( workspace_size>0 ) {
        buff.resize(workspace_size*nactive_tasks);
      }
      auto workspaceSize = workspace_size;
      auto buffer = buff.devicePtr();
      OpenACC::memoryManager<double>dev(ndev_per_task*nactive_tasks);
      auto derivatives = dev.devicePtr();
#pragma acc parallel loop present(input, t_actiondata) \
                           copyin(nactive_tasks, \
                                 ndev_per_task, \
                                 workspace_size)\
                        deviceptr(derivatives, \
                                  partialTaskList_data, \
                                  value_stash_data, \
                                  buffer) \
                          default(none)
      for(unsigned i=0; i<nactive_tasks; ++i) {
        std::size_t task_index = partialTaskList_data[i];
        std::size_t val_pos = task_index*input.nscalars;
        ParallelActionsOutput myout {input.nscalars,
                                     value_stash_data+val_pos,
                                     ndev_per_task,
                                     derivatives+ndev_per_task*i,
                                     workspaceSize,
                                     (workspaceSize>0)?
                                     buffer+workspaceSize*i
                                     :nullptr  };
        // Calculate the stuff in the loop for this action
        T::performTask( task_index, t_actiondata, input, myout );
      }
      vs.copyFromDevice(value_stash.data());
    }
    comm.Bcast( value_stash.data(), value_stash.size(), 0);
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
      std::vector<double> buffer( workspace_size );
      std::vector<double> derivatives( nderivatives_per_task );
      #pragma omp for nowait
      for(unsigned i=rank; i<nactive_tasks; i+=stride) {
        std::size_t task_index = partialTaskList[i];
        std::size_t val_pos = task_index*myinput.nscalars;
        ParallelActionsOutput myout( myinput.nscalars,
                                     value_stash.data()+val_pos,
                                     nderivatives_per_task,
                                     derivatives.data(),
                                     workspace_size,
                                     buffer.data() );
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
  //TODO: check if std::fill is faster (i get conflicting answers on the net)
  //std::fill (forcesForApply.begin(),forcesForApply.end(), 0.0);
  // Get all the input data so we can broadcast it to the GPU
  myinput.noderiv = false;
  // Retrieve the forces from the values
  action->transferForcesToStash( value_stash );

  if( useacc ) {
#ifdef __PLUMED_USE_OPENACC
    if (comm.Get_rank() == 0) {
      std::fill (omp_forces[0].begin(),omp_forces[0].end(), 0.0);
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

      const auto ndev_per_task = nderivatives_per_task;
      const auto nind_per_scalar = 2 + input.nderivatives_per_scalar;
      const auto nind_per_task = input.nscalars*nind_per_scalar;

      OpenACC::memoryManager<double> dev{ndev_per_task*nactive_tasks};
      auto derivatives = dev.devicePtr();
      OpenACC::memoryManager<std::size_t> ind{nind_per_task*nactive_tasks};
      auto indices = ind.devicePtr();
      OpenACC::memoryManager<double> vtmp{input.ncomponents*nactive_tasks};
      auto valstmp = vtmp.devicePtr();
      OpenACC::memoryManager<double> buff{};
      if ( workspace_size>0 ) {
        buff.resize(workspace_size*nactive_tasks);
      }
      auto workspaceSize = workspace_size;
      auto buffer = buff.devicePtr();

#pragma acc data present(input,t_actiondata) \
                  copyin(nactive_tasks, \
                         nderivPerComponent, \
                         forcesForApply_size, \
                         ndev_per_task, \
                         workspaceSize, \
                         omp_forces_size) \
               deviceptr(derivatives, \
                         indices, \
                         value_stash_data, \
                         partialTaskList_data, \
                         forcesForApply_data, \
                         valstmp, \
                         buffer) \
                 default(none)
      {
#pragma acc parallel loop
        for(unsigned i=0; i<nactive_tasks; ++i) {
          std::size_t task_index = partialTaskList_data[i];
          ParallelActionsOutput myout { input.nscalars,
                                        valstmp+input.nscalars*i,
                                        ndev_per_task,
                                        derivatives+ndev_per_task*i,
                                        workspaceSize,
                                        (workspaceSize>0)?
                                        buffer+workspaceSize*i
                                        :nullptr  };
          // Calculate the stuff in the loop for this action
          T::performTask( task_index, t_actiondata, input, myout );

          // If this is a matrix this returns a number that isn't one as we have to loop over the columns
          std::size_t nfpt = T::getNumberOfValuesPerTask( task_index, t_actiondata );
          for(unsigned j=0; j<nfpt; ++j) {
            // Create a force index holder
            ForceIndexHolder force_indices( input.ncomponents, input.nderivatives_per_scalar, indices + i*nind_per_task + j*input.ncomponents*nind_per_scalar );
            // Get the indices for forces
            T::getForceIndices( task_index,
                                j,
                                forcesForApply_size,
                                t_actiondata,
                                input,
                                force_indices );

            // Create a force input object
            ForceInput finput( myinput.nscalars,
                               value_stash_data+input.nscalars*task_index + j*input.ncomponents,
                               input.nderivatives_per_scalar,
                               derivatives + i*ndev_per_task + j*input.ncomponents*input.nderivatives_per_scalar );

            // Gather forces that can be gathered locally
            gatherThreadSafeForces( input,
                                    force_indices,
                                    finput,
                                    View<double>(forcesForApply_data,forcesForApply_size));
          }
        }

#pragma acc parallel loop
        for(unsigned v=input.threadunsafe_forces_start; v<forcesForApply_size; ++v) {
          double tmp = 0.0;
#pragma acc loop reduction(+:tmp)
          for(unsigned t=0; t<nactive_tasks; ++t) {
            std::size_t task_index = partialTaskList_data[t];
            std::size_t nfpt = T::getNumberOfValuesPerTask( task_index, t_actiondata );
            for(unsigned j=0; j<nfpt; ++j) {
              ForceIndexHolder force_indices( input.ncomponents, input.nderivatives_per_scalar, indices + t*nind_per_task + j*input.ncomponents*nind_per_scalar );
              auto fdata = ForceInput { input.nscalars,
                                        value_stash_data+input.nscalars*task_index + j*input.ncomponents,
                                        input.nderivatives_per_scalar,
                                        derivatives + i*ndev_per_task + j*input.ncomponents*input.nderivatives_per_scalar }
              for(unsigned i=0; i<input.ncomponents; ++i) {
                for(unsigned d=force_indices.threadsafe_derivatives_end[i]; d<force_indices.tot_indices[i]; ++d) {
                  if( force_indices.indices[i][d]==v ) {
                    tmp += fdata.force[i]*fdata.deriv[i][d];
                    break;
                  }
                }
              }
            }
          }
          forcesForApply_data[v] = tmp;
        }
      }
      ffa.copyFromDevice(forcesForApply.data());
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
        omp_forces[t].assign( omp_forces[t].size(), 0.0 );
        ForceOutput forces{ omp_forces[t], forcesForApply };
        std::vector<double> buffer( workspace_size );
        std::vector<double> fake_vals( myinput.nscalars );
        std::vector<double> derivatives( nderivatives_per_task );
        std::vector<std::size_t> indices( (2+myinput.nderivatives_per_scalar)*myinput.ncomponents );
        ForceIndexHolder force_indices( myinput.ncomponents, myinput.nderivatives_per_scalar, indices.data() );
        #pragma omp for nowait
        for(unsigned i=rank; i<nactive_tasks; i+=stride) {
          std::size_t task_index = partialTaskList[i];
          ParallelActionsOutput myout( myinput.nscalars,
                                       fake_vals.data(),
                                       derivatives.size(),
                                       derivatives.data(),
                                       workspace_size,
                                       buffer.data() );
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
            ForceInput finput( myinput.nscalars,
                               value_stash.data()+myinput.nscalars*task_index + j*myinput.ncomponents,
                               myinput.nderivatives_per_scalar,
                               derivatives.data() + j*myinput.ncomponents*myinput.nderivatives_per_scalar );

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
      const ForceIndexHolder force_indices,
      const ForceInput& fdata,
      View<double,helpers::dynamic_extent> forces ) {
    for(unsigned i=0; i<input.ncomponents; ++i) {
      double ff = fdata.force[i];
      for(unsigned j=0; j<force_indices.threadsafe_derivatives_end[i]; ++j) {
        forces[ force_indices.indices[i][j] ] += ff*fdata.deriv[i][j];
      }
    }
  }

  template <class T>
  void ParallelTaskManager<T>::gatherThreadUnsafeForces( const ParallelActionsInput& input,
      const ForceIndexHolder force_indices,
      const ForceInput& fdata,
      View<double,helpers::dynamic_extent> forces ) {
    for(unsigned i=0; i<input.ncomponents; ++i) {
      double ff = fdata.force[i];
      for(unsigned j=force_indices.threadsafe_derivatives_end[i]; j<force_indices.tot_indices[i]; ++j) {
        forces[ force_indices.indices[i][j] - input.threadunsafe_forces_start ] += ff*fdata.deriv[i][j];
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
