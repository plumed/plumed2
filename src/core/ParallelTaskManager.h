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

template<typename precision>
struct ParActionsInput {
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
  precision *inputdata{nullptr};
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
  static ParActionsInput create( const Pbc& box ) {
    auto toret=ParActionsInput();
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

template<typename precision>
inline void ParActionsInput<precision>::setupArguments( const ArgumentsBookkeeping& ab ) {
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

  template<typename precision>
  static ArgumentBookeepingHolder create ( std::size_t argno, const ParActionsInput<precision>& inp ) {
    return ArgumentBookeepingHolder{
      inp.ranks[argno], // rank
      inp.ncols[argno], // ncols
      inp.argstarts[argno], // start
      View<const std::size_t>(inp.shapedata + inp.shapestarts[argno], inp.ranks[argno] ), // shape
      View<const std::size_t>(inp.bookeeping + inp.bookstarts[argno], inp.booksizes[argno] )  // bookeeping
    };
  }
};

template<typename precision>
struct ParActionsOutput {
  View<precision> values;
  View<precision> derivatives;
  View<precision> buffer;

  static ParActionsOutput create( std::size_t ncomp,
                                  precision* v,
                                  std::size_t ndev,
                                  precision* d,
                                  std::size_t nb,
                                  precision* b ) {
    return ParActionsOutput{
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

  template<typename precision>
  static ForceIndexHolder create(const ParActionsInput<precision>& inp,
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
  template<typename precision>
  static size_t indexesPerScalar(const ParActionsInput<precision>& inp) {
    return indexesPerScalar(inp.ncomponents,
                            inp.nderivatives_per_scalar);
  }
};

template <typename precision>
class ForcesInput {
public:
  View<precision> force;
  View2D<precision> deriv;
  static ForcesInput create( std::size_t nv, precision* f, std::size_t nd, precision* d ) {
    return ForcesInput{
      View{f,nv},     //force
      View2D{d,nv,nd} //deriv
    };
  }
};

//There is no need to pass this as reference:
template <typename precision>
struct ForcesOutput {
  //I would suggest to invert the name or to be clearer
// like something that recalls that "thread_safe" will be need to be reducted (hence it is NOT thread safe)
  View<precision> thread_safe;
  //these are the forces that we promise will not provoke races
  View<precision> thread_unsafe;
  //const T* is a ptr to const T
  //T* const is a conts ptr to a modifiable T
  static ForcesOutput create(std::vector<precision>& reduced, std::vector<precision>& notReduced) {
    return ForcesOutput{
      View{reduced.data(),reduced.size()},      // thread_safe
      View{notReduced.data(),notReduced.size()} // thread_unsafe
    };
  }
  static ForcesOutput create(precision* reduced, size_t rs, precision* notReduce, size_t nrsz) {
    return ForcesOutput{
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
//            std::declval<const ParActionsInput& >(),
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
//            std::declval<const ParActionsInput& >(),
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

template <typename CV, typename=void>
struct cvprecision {
  typedef double type;
};

template <typename CV>
struct cvprecision<CV,std::void_t<typename CV::precision>> {
  typedef typename CV::precision type;
};

template <typename CV>
using cvprecision_t = typename cvprecision<CV>::type;

template <class T>
class ParallelTaskManager {
public:
  using input_type= typename T::input_type;
  using precision = cvprecision_t<T>;
  typedef ParActionsInput<precision> ParallelActionsInput;
  typedef ParActionsOutput<precision> ParallelActionsOutput;
  typedef ForcesInput<precision> ForceInput;
  typedef ForcesOutput<precision> ForceOutput;
//  static constexpr bool has_custom_gather=PTMUtils::has_gatherForces_custom<T>;
//  static constexpr bool has_GPU_gather=PTMUtils::has_gatherForces_GPU<T>;
//  static constexpr size_t virialSize = PTMUtils::virialSize<T>;
protected:
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
  std::vector<precision> value_stash;
/// A tempory set of vectors for holding forces over threads
  std::vector<std::vector<precision> > omp_forces;
/// This structs is used to pass data between the parallel interface and the function caller
  ParallelActionsInput myinput;
  ArgumentsBookkeeping argumentsMap;
//this holds the data for myinput that will be passed though myinput
  std::vector<precision> input_buffer;
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
                                      View<precision> forces );
/// This is used to gather forces that are not thread safe
  static void gatherThreadUnsafeForces( const ParallelActionsInput& input,
                                        const ForceIndexHolder& force_indices,
                                        const ForceInput& fdata,
                                        View<precision> forces );
};

struct defaultPTM {
  template <typename ACC>
  using PTM=ParallelTaskManager<ACC>;
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
  std::fill (value_stash.begin(),value_stash.end(), precision(0.0));
  if( useacc ) {
    // there should not be any avaiable path to get here
    plumed_merror("to use the USEGPU flag you need to LOAD a plugin like \"openaccPTM\"");
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
      std::vector<precision> buffer( workspace_size );
      std::vector<precision> derivatives( nderivatives_per_task );
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
  action->transferStashToValues( value_stash );
}

template<typename prec>
struct forceData {
  std::vector<prec> ffa;
  forceData( std::vector<double>& forcesForApply ):
    ffa(forcesForApply.size()),destination(forcesForApply) {}
  void update() {
    std::copy(ffa.begin(),ffa.end(),destination.begin());
  }
private:
  std::vector<double>& destination;
};

template<>
struct forceData<double> {
  PLMD::View<double> ffa;
  forceData( std::vector<double>& forcesForApply ):
    ffa(forcesForApply.data(),forcesForApply.size()) {}
  void update() {
  }
};
template <class T>
void ParallelTaskManager<T>::applyForces( std::vector<double>& forcesForApply ) {
  // Get the list of active tasks
  std::vector<unsigned> & partialTaskList= action->getListOfActiveTasks( action ) ;
  unsigned nactive_tasks=partialTaskList.size();
  forceData<precision> forces(forcesForApply);
  // Clear force buffer
  std::fill (forces.ffa.begin(),forces.ffa.end(), precision(0.0));
  // Get all the input data so we can broadcast it to the GPU
  myinput.noderiv = false;
  // Retrieve the forces from the values
  action->transferForcesToStash( value_stash );

  if( useacc ) {
    // there should not be any avaiable path to get here
    plumed_merror("to use the USEGPU flag you need to LOAD a plugin like \"openaccPTM\"");
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
      std::vector<precision> buffer( workspace_size );
      std::vector<precision> fake_vals( myinput.sizeOfFakeVals() );
      std::vector<precision> derivatives( nderivatives_per_task );
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
                              forces.ffa.size(),
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
                                  View<double>(forces.ffa.data(),
                                               forces.ffa.size()) );

          // Gather forces that are not thread safe
          gatherThreadUnsafeForces( myinput,
                                    force_indices,
                                    finput,
                                    View<precision>(omp_forces[t].data(),
                                        omp_forces[t].size()) );
        }
      }

      #pragma omp critical
      gatherThreads( ForceOutput::create(omp_forces[t], forcesForApply ) );
    }
    forces.update();
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
    View<precision> forces ) {
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
    View<precision> forces ) {
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
