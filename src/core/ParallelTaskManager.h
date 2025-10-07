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
  std::size_t s = 0;
  std::size_t ts = 0;
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
  std::size_t ncomponents{0};
  /// The number of scalars we are calculating for each task
  unsigned nscalars{0};
  /// The number of force scalars for each task
  unsigned nforcescalars{0};
  /// Number of derivatives for each scalar being calculated
  std::size_t nderivatives_per_scalar{0};
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
  static ParallelActionsInput create( const Pbc& box ) {
    auto toret=ParallelActionsInput();
    toret.pbc=&box;
    return toret;
  }

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
  unsigned sizeOfFakeVals() const {
    return nscalars;
  }
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

struct ArgumentBookeepingHolder {
  std::size_t rank;
  std::size_t ncols;
  std::size_t start;
  View<const std::size_t> shape;
  View<const std::size_t> bookeeping;

  static ArgumentBookeepingHolder create ( std::size_t argno, const ParallelActionsInput& inp ) {
    return ArgumentBookeepingHolder{
      inp.ranks[argno], // rank
      inp.ncols[argno], // ncols
      inp.argstarts[argno], // start
      View<const std::size_t>(inp.shapedata + inp.shapestarts[argno], inp.ranks[argno] ), // shape
      View<const std::size_t>(inp.bookeeping + inp.bookstarts[argno], inp.booksizes[argno] )  // bookeeping
    };
  }
};

struct ParallelActionsOutput {
  View<double> values;
  View<double> derivatives;
  View<double> buffer;

  static ParallelActionsOutput create( std::size_t ncomp, double* v, std::size_t ndev, double* d, std::size_t nb, double* b ) {
    return ParallelActionsOutput{
      View{v,ncomp}, //values
      View{d,ndev},  // derivatives
      View{b,nb}     // buffer
    };
  }
};

struct ForceIndexHolder {
  View<std::size_t> threadsafe_derivatives_end;
  View<std::size_t> tot_indices;
  View2D<std::size_t> indices;

/// Constructs a ForceIndexHolder object.
/// \param nc Definition here (number of components?)
/// \param nd Definition here (number of derivatives?)
/// \param ind Pointer to an array storing index data. It should have a size of at least 2*nc + nc*nd.
  static ForceIndexHolder create(const std::size_t nc,
                                 const std::size_t nd,
                                 std::size_t* ind ) {
    return ForceIndexHolder{
      View{ind,nc},        // threadsafe_derivatives_end
      View{ind+nc,nc},     // tot_indices
      View2D{ind+2*nc,nc,nd} // indices
    };
  }
  static ForceIndexHolder create(const ParallelActionsInput& inp,
                                 std::size_t* ind ) {
    return create(inp.ncomponents,
                  inp.nderivatives_per_scalar,ind);
  }

/// \brief Returns the number of indexes needed by a ForceIndexHolder
///        for each scalar.
///
/// \param nc Definition here (number of components?)
/// \param nd Definition here (number of derivatives?)
///
/// \return the number of indexes needed by a ForceIndexHolder
///         for each scalar.
  static size_t indexesPerScalar(const std::size_t nc, const std::size_t nd) {
    return nc       // threadsafe_derivatives_end
           + nc     // tot_indices
           + nc*nd; // indices
  }
  static size_t indexesPerScalar(const ParallelActionsInput& inp) {
    return indexesPerScalar(inp.ncomponents,
                            inp.nderivatives_per_scalar);
  }
};

class ForceInput {
public:
  View<double> force;
  View2D<double> deriv;
  static ForceInput create( std::size_t nv, double* f, std::size_t nd, double* d ) {
    return ForceInput{
      View{f,nv},     //force
      View2D{d,nv,nd} //deriv
    };
  }
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
  static ForceOutput create(std::vector<double>& reduced, std::vector<double>& notReduced) {
    return ForceOutput{
      View{reduced.data(),reduced.size()},      // thread_safe
      View{notReduced.data(),notReduced.size()} // thread_unsafe
    };
  }
  static ForceOutput create(double* reduced, size_t rs, double* notReduce, size_t nrsz) {
    return ForceOutput{
      View{reduced,rs},    // thread_safe
      View{notReduce,nrsz} // thread_unsafe
    };
  }
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
/// True if not using MPI for parllisation
  bool serial;
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
  const input_type& getActionInput() const ;
/// Is the calculation running in serial
  bool runInSerial() const {
    return serial;
  }
/// This runs all the tasks
  void runAllTasks();
/// Apply the forces on the parallel object
  void applyForces( std::vector<double>& forcesForApply );
/// This is used to gather forces that are thread safe
  static void gatherThreadSafeForces( const ParallelActionsInput& input,
                                      const ForceIndexHolder& force_indices,
                                      const ForceInput& fdata,
                                      View<double> forces );
/// This is used to gather forces that are not thread safe
  static void gatherThreadUnsafeForces( const ParallelActionsInput& input,
                                        const ForceIndexHolder& force_indices,
                                        const ForceInput& fdata,
                                        View<double> forces );
};

template <class T>
void ParallelTaskManager<T>::registerKeywords( Keywords& keys ) {
  keys.addFlag("SERIAL",false,"do the calculation in serial.  Do not parallelize");
  keys.addFlag("USEGPU",false,"run this calculation on the GPU");
  keys.addLinkInDocForFlag("USEGPU","gpu.md");
  keys.addLinkInDocForFlag("SERIAL", "actions.md");
}

template <class T>
ParallelTaskManager<T>::ParallelTaskManager(ActionWithVector* av):
  action(av),
  comm(av->comm),
  ismatrix(false),
  useacc(false),
  nderivatives_per_task(0),
  nthreaded_forces(0),
  myinput(ParallelActionsInput::create(av->getPbc())),
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
  action->parseFlag("SERIAL",serial);
  if( serial ) {
    action->log.printf("  not using MPI to parallelise this action\n");
  }
}

template <class T>
std::size_t ParallelTaskManager<T>::getValueStashSize() const {
  std::size_t valuesize=0;
  for(unsigned i=0; i<action->getNumberOfComponents(); ++i) {
    const Value* mycomp = action->getConstPntrToComponent(i);
    if( mycomp->hasDerivatives() ) {
      valuesize += mycomp->getNumberOfStoredValues()*(1+action->getNumberOfDerivatives());
    } else {
      valuesize += mycomp->getNumberOfStoredValues();
    }
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
  myinput.nforcescalars = 0;
  for(unsigned i=0; i<myinput.ncomponents; ++i) {
    if( (action->copyOutput(i))->hasDerivatives() ) {
      myinput.nscalars += 1 + action->getNumberOfDerivatives();
      myinput.nforcescalars += 1;
    } else if( (action->copyOutput(i))->getRank()==1 ) {
      myinput.nscalars += 1;
      myinput.nforcescalars += 1;
    } else if( (action->copyOutput(i))->getRank()==2 ) {
      if( ntasks==(action->copyOutput(i))->getShape()[0] ) {
        myinput.nscalars += (action->copyOutput(i))->getNumberOfColumns();
        myinput.nforcescalars += (action->copyOutput(i))->getNumberOfColumns();
      } else {
        myinput.nscalars += 1;
        myinput.nforcescalars += 1;
      }
    }
  }
  myinput.nderivatives_per_scalar = nder;
  nderivatives_per_task = nder*myinput.nforcescalars;
  value_stash.resize( getValueStashSize() );
  myinput.threadunsafe_forces_start = action->getNumberOfForceDerivatives() - nforce_ts;
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
const typename ParallelTaskManager<T>::input_type& ParallelTaskManager<T>::getActionInput() const {
  return actiondata;
}

template <class T>
void ParallelTaskManager<T>::setWorkspaceSize( std::size_t size ) {
  workspace_size = size;
}

#ifdef __PLUMED_USE_OPENACC
//use the __PLUMED_USE_OPENACC_TASKSMINE macro to debug the ptm ins a single file
//so that compiling witha a small modification will be faster (the ptm is included nearly everywhere)
#ifndef __PLUMED_USE_OPENACC_TASKSMINE
template <class T>
void runAllTasksACC(typename T::input_type actiondata,
                    ParallelActionsInput myinput,
                    std::vector<double>& value_stash,
                    const std::vector<unsigned> & partialTaskList,
                    const unsigned nactive_tasks,
                    const std::size_t nderivatives_per_task,
                    const std::size_t workspace_size
                   ) {
  auto myinput_acc = OpenACC::fromToDataHelper(myinput);
  auto actiondata_acc = OpenACC::fromToDataHelper(actiondata);

  //template type is deduced
  OpenACC::memoryManager vs{value_stash};
  auto value_stash_data = vs.devicePtr();

  OpenACC::memoryManager ptl{partialTaskList};
  auto partialTaskList_data = ptl.devicePtr();

  OpenACC::memoryManager<double> buff{workspace_size*nactive_tasks};

  auto buffer = buff.devicePtr();
  OpenACC::memoryManager<double> dev(nderivatives_per_task*nactive_tasks);
  auto derivatives = dev.devicePtr();
#pragma acc parallel loop present(myinput, actiondata) \
                           copyin(nactive_tasks, \
                                 nderivatives_per_task, \
                                 workspace_size)\
                        deviceptr(derivatives, \
                                  partialTaskList_data, \
                                  value_stash_data, \
                                  buffer) \
                          default(none)
  for(unsigned i=0; i<nactive_tasks; ++i) {
    std::size_t task_index = partialTaskList_data[i];
    std::size_t val_pos = task_index*myinput.nscalars;
    auto myout = ParallelActionsOutput::create (myinput.nscalars,
                 value_stash_data+val_pos,
                 nderivatives_per_task,
                 derivatives+nderivatives_per_task*i,
                 workspace_size,
                 (workspace_size>0)?
                 buffer+workspace_size*i
                 :nullptr  );
    // Calculate the stuff in the loop for this action
    T::performTask( task_index, actiondata, myinput, myout );
  }
  vs.copyFromDevice(value_stash.data());
}
#else
template <class T>
void runAllTasksACC(typename T::input_type actiondata,
                    ParallelActionsInput myinput,
                    std::vector<double>& value_stash,
                    const std::vector<unsigned> & partialTaskList,
                    const unsigned nactive_tasks,
                    const std::size_t nderivatives_per_task,
                    const std::size_t workspace_size
                   ) ;
#endif //__PLUMED_USE_OPENACC_TASKSMINE
#endif //__PLUMED_USE_OPENACC

template <class T>
void ParallelTaskManager<T>::runAllTasks() {
  // Get the list of active tasks
  std::vector<unsigned> & partialTaskList( action->getListOfActiveTasks( action ) );
  unsigned nactive_tasks=partialTaskList.size();
  // Get all the input data so we can broadcast it to the GPU
  myinput.noderiv = true;
  for(unsigned i=0; i<action->getNumberOfComponents(); ++i) {
    if( (action->getConstPntrToComponent(i))->hasDerivatives() ) {
      myinput.noderiv=false;
      break;
    }
  }
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
      runAllTasksACC<T>(
        actiondata,
        myinput,
        value_stash,
        partialTaskList,
        nactive_tasks,
        nderivatives_per_task,
        workspace_size
      );
    }
    comm.Bcast( value_stash.data(), value_stash.size(), 0);
#else
    plumed_merror("cannot use USEGPU flag if PLUMED has not been compiled with openACC");
#endif
  } else {
    // Get the MPI details
    unsigned stride=comm.Get_size();
    unsigned rank=comm.Get_rank();
    if(serial) {
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
      std::vector<double> buffer( workspace_size );
      std::vector<double> derivatives( nderivatives_per_task );
      #pragma omp for nowait
      for(unsigned i=rank; i<nactive_tasks; i+=stride) {
        std::size_t task_index = partialTaskList[i];
        std::size_t val_pos = task_index*myinput.nscalars;
        auto myout = ParallelActionsOutput::create ( myinput.nscalars,
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
    if( !serial ) {
      comm.Sum( value_stash );
    }
  }
// And transfer the value to the output values
  action->transferStashToValues( partialTaskList, value_stash );
}

#ifdef __PLUMED_USE_OPENACC
//use the __PLUMED_USE_OPENACC_FORCESMINE macro to debug the ptm ins a single file
//so that compiling witha a small modification will be faster (the ptm is included nearly everywhere)
#ifndef __PLUMED_USE_OPENACC_FORCESMINE
template <class T>
void applyForcesWithACC(PLMD::View<double> forcesForApply,
                        typename T::input_type actiondata,
                        ParallelActionsInput myinput,
                        const std::vector<double>& value_stash,
                        const std::vector<unsigned> & partialTaskList,
                        const unsigned nactive_tasks,
                        const std::size_t nderivatives_per_task,
                        const std::size_t workspace_size
                       ) {
  auto myinput_acc = OpenACC::fromToDataHelper(myinput);
  auto actiondata_acc = OpenACC::fromToDataHelper(actiondata);

  //template type is deduced
  OpenACC::memoryManager vs{value_stash};
  auto value_stash_data = vs.devicePtr();

  OpenACC::memoryManager ptl{partialTaskList};
  auto partialTaskList_data = ptl.devicePtr();

  OpenACC::memoryManager ffa {forcesForApply};
  auto forcesForApply_data = ffa.devicePtr();
  const auto forcesForApply_size = ffa.size();
  const auto nind_per_scalar = ForceIndexHolder::indexesPerScalar(myinput);
  //nscalars is >=ncomponents (see setupParallelTaskManager )
  const auto nind_per_task = nind_per_scalar*myinput.nscalars;

  OpenACC::memoryManager<double> dev{nderivatives_per_task*nactive_tasks};
  auto derivatives = dev.devicePtr();
  OpenACC::memoryManager<std::size_t> ind{nind_per_task*nactive_tasks};
  auto indices = ind.devicePtr();
  OpenACC::memoryManager<double> vtmp{myinput.sizeOfFakeVals()*nactive_tasks};
  auto valstmp = vtmp.devicePtr();
  OpenACC::memoryManager<double> buff{workspace_size*nactive_tasks};
  auto buffer = buff.devicePtr();

#define forces_indicesArg(taskID,scalarID) ForceIndexHolder::create(myinput, \
                          indices + taskID*nind_per_task + scalarID*nind_per_scalar)
#define derivativeDrift(taskID,scalarID)  taskID*nderivatives_per_task \
                           + scalarID*myinput.ncomponents*myinput.nderivatives_per_scalar
#define stashDrift(taskID,scalarID) taskID*myinput.nscalars \
                           + scalarID*myinput.ncomponents

#pragma acc data present(myinput,actiondata) \
                  copyin(nactive_tasks, \
                         forcesForApply_size, \
                         nderivatives_per_task, nind_per_task,nind_per_scalar, \
                         workspace_size) \
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
    for(unsigned t=0; t<nactive_tasks; ++t) {
      std::size_t task_index = partialTaskList_data[t];
      auto myout = ParallelActionsOutput::create( myinput.nscalars,
                   valstmp+myinput.nscalars*t,
                   nderivatives_per_task,
                   derivatives+nderivatives_per_task*t,
                   workspace_size,
                   (workspace_size>0)?buffer+workspace_size*t:nullptr);
      // Calculate the stuff in the loop for this action
      T::performTask( task_index, actiondata, myinput, myout );
      // If this is a matrix this returns a number that isn't one as we have to loop over the columns
      const std::size_t nvpt = T::getNumberOfValuesPerTask( task_index, actiondata );
#pragma acc loop seq
      for(unsigned vID=0; vID<nvpt; ++vID) {
        auto force_indices = forces_indicesArg(t,vID);
        // Create a force index holder
        // Get the indices for forces
        T::getForceIndices( task_index,
                            vID,
                            forcesForApply_size,
                            actiondata,
                            myinput,
                            force_indices );

        // Create a force input object
        auto finput = ForceInput::create ( myinput.nscalars,
                                           value_stash_data + stashDrift(task_index,vID),
                                           myinput.nderivatives_per_scalar,
                                           derivatives + derivativeDrift(t,vID));

        // Gather forces that can be gathered locally
        ParallelTaskManager<T>::gatherThreadSafeForces( myinput,
            force_indices,
            finput,
            View<double>(forcesForApply_data,
                         forcesForApply_size));
      }
    }

#pragma acc parallel loop
    for(unsigned v=myinput.threadunsafe_forces_start; v<forcesForApply_size; ++v) {
      double tmp = 0.0;
#pragma acc loop reduction(+:tmp)
      for(unsigned t=0; t<nactive_tasks; ++t) {
        const std::size_t task_index = partialTaskList_data[t];
        const std::size_t nvpt = T::getNumberOfValuesPerTask( task_index, actiondata );
        for(unsigned vID=0; vID<nvpt; ++vID) {
          auto force_indices = forces_indicesArg(t,vID);

          auto fdata = ForceInput::create( myinput.nscalars,
                                           value_stash_data + stashDrift(task_index,vID),
                                           myinput.nderivatives_per_scalar,
                                           derivatives + derivativeDrift(t,vID));
          for(unsigned i=0; i<myinput.ncomponents; ++i) {
            const double ff = fdata.force[i];
            for(unsigned d=force_indices.threadsafe_derivatives_end[i];
                d<force_indices.tot_indices[i]; ++d) {
              if( force_indices.indices[i][d]==v ) {
                tmp += ff*fdata.deriv[i][d];
                // break;
              }
            }
          }
        }
      }
      forcesForApply_data[v] = tmp;
    }
  }
#undef forces_indicesArg
#undef derivativeDrift
#undef stashDrift
  ffa.copyFromDevice(forcesForApply.data());
}
#else
template <class T>
void applyForcesWithACC(PLMD::View<double> forcesForApply,
                        typename T::input_type actiondata,
                        ParallelActionsInput myinput,
                        const std::vector<double>& value_stash,
                        const std::vector<unsigned> & partialTaskList,
                        const unsigned nactive_tasks,
                        const std::size_t nderivatives_per_task,
                        const std::size_t workspace_size
                       );
#endif //__PLUMED_USE_OPENACC_FORCESMINE
#endif //__PLUMED_USE_OPENACC
template <class T>
void ParallelTaskManager<T>::applyForces( std::vector<double>& forcesForApply ) {
  // Get the list of active tasks
  std::vector<unsigned> & partialTaskList= action->getListOfActiveTasks( action ) ;
  unsigned nactive_tasks=partialTaskList.size();
  // Clear force buffer
  forcesForApply.assign( forcesForApply.size(), 0.0 );
  //TODO: check if std::fill is faster (i get conflicting answers on the net)
  //std::fill (forcesForApply.begin(),forcesForApply.end(), 0.0);
  // Get all the input data so we can broadcast it to the GPU
  myinput.noderiv = false;
  // Retrieve the forces from the values
  action->transferForcesToStash( partialTaskList, value_stash );

  if( useacc ) {
#ifdef __PLUMED_USE_OPENACC
    std::fill (omp_forces[0].begin(),omp_forces[0].end(), 0.0);
    if (comm.Get_rank() == 0) {
      applyForcesWithACC<T>(
        PLMD::View<double> { forcesForApply.data(), forcesForApply.size() },
        actiondata,
        myinput,
        value_stash,
        partialTaskList,
        nactive_tasks,
        nderivatives_per_task,
        workspace_size
      );
    }
#else
    plumed_merror("cannot use USEGPU flag if PLUMED has not been compiled with openACC");
#endif
  } else {
    // Get the MPI details
    unsigned stride=comm.Get_size();
    unsigned rank=comm.Get_rank();
    if(serial) {
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
      std::vector<double> buffer( workspace_size );
      std::vector<double> fake_vals( myinput.sizeOfFakeVals() );
      std::vector<double> derivatives( nderivatives_per_task );
      std::vector<std::size_t> indices(ForceIndexHolder::indexesPerScalar(myinput));

      auto force_indices = ForceIndexHolder::create( myinput,indices.data() );
      #pragma omp for nowait
      for(unsigned i=rank; i<nactive_tasks; i+=stride) {
        std::size_t task_index = partialTaskList[i];
        auto myout = ParallelActionsOutput::create( myinput.nscalars,
                     fake_vals.data(),
                     derivatives.size(),
                     derivatives.data(),
                     workspace_size,
                     buffer.data() );
        // Calculate the stuff in the loop for this action
        T::performTask( task_index, actiondata, myinput, myout );

        // If this is a matrix this returns a number that isn't one as we have to loop over the columns
        const std::size_t nvpt = T::getNumberOfValuesPerTask( task_index, actiondata );
        for(unsigned j=0; j<nvpt; ++j) {
          // Get the force indices
          T::getForceIndices( task_index,
                              j,
                              forcesForApply.size(),
                              actiondata,
                              myinput,
                              force_indices );
          // Create a force input object
          auto finput=ForceInput::create( myinput.nforcescalars,
                                          value_stash.data()
                                          + myinput.nforcescalars*task_index
                                          + j*myinput.ncomponents,
                                          myinput.nderivatives_per_scalar,
                                          derivatives.data()
                                          + j*myinput.ncomponents*myinput.nderivatives_per_scalar );

          // Gather forces that are thread safe
          gatherThreadSafeForces( myinput,
                                  force_indices,
                                  finput,
                                  View<double>(forcesForApply.data(),
                                               forcesForApply.size()) );

          // Gather forces that are not thread safe
          gatherThreadUnsafeForces( myinput,
                                    force_indices,
                                    finput,
                                    View<double>(omp_forces[t].data(),
                                                 omp_forces[t].size()) );
        }
      }

      #pragma omp critical
      gatherThreads( ForceOutput::create(omp_forces[t], forcesForApply ) );
    }
    // MPI Gather everything (this must be extended to the gpu thing, after makning it mpi-aware)
    if( !serial ) {
      comm.Sum( forcesForApply );
    }
  }
}

template <class T>
void ParallelTaskManager<T>::gatherThreadSafeForces( const ParallelActionsInput& input,
    const ForceIndexHolder& force_indices,
    const ForceInput& fdata,
    View<double> forces ) {
  for(unsigned i=0; i<input.ncomponents; ++i) {
    double ff = fdata.force[i];
    for(unsigned j=0; j<force_indices.threadsafe_derivatives_end[i]; ++j) {
      forces[ force_indices.indices[i][j] ] += ff*fdata.deriv[i][j];
    }
  }
}

template <class T>
void ParallelTaskManager<T>::gatherThreadUnsafeForces(const ParallelActionsInput& input,
    const ForceIndexHolder& force_indices,
    const ForceInput& fdata,
    View<double> forces ) {
  for(unsigned i=0; i<input.ncomponents; ++i) {
    const double ff = fdata.force[i];
    for(unsigned d=force_indices.threadsafe_derivatives_end[i];
        d<force_indices.tot_indices[i]; ++d) {
      forces[ force_indices.indices[i][d] - input.threadunsafe_forces_start ]
      += ff*fdata.deriv[i][d];
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
