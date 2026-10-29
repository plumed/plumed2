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
#ifndef __PLUMED_acc_ParallelTaskManager_h
#define __PLUMED_acc_ParallelTaskManager_h
#include "plumed/core/ActionWithVector.h"
#include "plumed/core/ActionWithMatrix.h"
#include "plumed/core/ParallelTaskManager.h"
#include "plumed/tools/Communicator.h"
#include "plumed/tools/OpenMP.h"
#include "plumed/tools/View.h"
#include "plumed/tools/View2D.h"

#include "plumed/tools/ColvarOutput.h"
#include "plumed/tools/OpenACC.h"
#include "plumed/tools/Pbc.h"

#include <algorithm>
#include <vector>

namespace PLMD {

/// Device buffers used by the OpenACC task manager.
/// They are kept between calls (one set per action) and grown when needed, so
/// that the device memory is not allocated and freed at every step.
/// Allocation happens on first use, so ranks that do not run on the GPU do not allocate.
template <typename precision>
struct AccDeviceWorkspace {
  /// device copies of the data pointed by ParActionsInput (see toDeviceInput)
  OpenACC::memoryManager<precision> inputdata;
  OpenACC::memoryManager<Pbc> pbc;
  OpenACC::memoryManager<std::size_t> argbook;
  /// host staging area for the argument bookkeeping, packed in a single array
  std::vector<std::size_t> argbook_host;
  OpenACC::memoryManager<precision> value_stash;
  OpenACC::memoryManager<unsigned> taskList;
  OpenACC::memoryManager<precision> buffer;
  OpenACC::memoryManager<precision> derivatives;
  OpenACC::memoryManager<precision> forces;
  OpenACC::memoryManager<std::size_t> indices;
  OpenACC::memoryManager<precision> valstmp;
  /// per-task partial sums of the forces on the last entries of the force array
  OpenACC::memoryManager<precision> sharedForces;
};

/// The async queue used by the task manager. All the transfers and the kernels
/// of a call are queued on it and the host waits only once, at the end.
inline constexpr int accQueue=1;

template <class T>
class AccParallelTaskManager : public ParallelTaskManager<T> {
  using ParallelTaskManager<T>::action;
  using ParallelTaskManager<T>::myinput;
  using ParallelTaskManager<T>::argumentsMap;
  using ParallelTaskManager<T>::value_stash;
  using ParallelTaskManager<T>::comm;
  using ParallelTaskManager<T>::actiondata;
  using ParallelTaskManager<T>::nderivatives_per_task;
  using ParallelTaskManager<T>::workspace_size;
  using ParallelTaskManager<T>::omp_forces;
  using ParallelTaskManager<T>::input_buffer;
  using ParallelTaskManager<T>::serial;
  using ParallelTaskManager<T>::useacc;
  using ParallelTaskManager<T>::getValueStashSize;
public:
  typedef typename ParallelTaskManager<T>::ParallelActionsInput ParallelActionsInput;
  typedef typename ParallelTaskManager<T>::ParallelActionsOutput ParallelActionsOutput;
  typedef typename ParallelTaskManager<T>::ForceInput ForceInput;
  typedef typename ParallelTaskManager<T>::ForceOutput ForceOutput;
  typedef typename ParallelTaskManager<T>::precision precision;
private:
/// Device buffers kept between steps
  AccDeviceWorkspace<precision> devws;
public:
  static void registerKeywords( Keywords& keys ) {
    ParallelTaskManager<T>::registerKeywords(keys);
  }
  AccParallelTaskManager(ActionWithVector* av):
    PLMD::ParallelTaskManager<T>(av) {
    useacc=true;
  }
/// This runs all the tasks
  void runAllTasks();
/// Apply the forces on the parallel object
  void applyForces( std::vector<double>& forcesForApply );
};

struct ACCPTM {
  template <typename T>
  using PTM=AccParallelTaskManager<T>;
};

/// Returns a copy of the input whose pointers are device addresses.
/// The data (input data, Pbc and argument bookkeeping) is copied in the
/// device buffers of the workspace, with one transfer per array: the
/// bookkeeping arrays are packed in a single buffer. The returned object can
/// be passed to the kernels as firstprivate, i.e. as a kernel argument.
/// The copies are queued on accQueue: the host data must not change until
/// the queue is waited.
/// Before this, the input was put on the device with "enter data" and removed
/// with "exit data" at every call: each pointer member costs a transfer and
/// an attach/detach (another small transfer), i.e. about 40 transfers of a
/// few microseconds each per call, which dominated the time of small actions.
template <typename precision>
ParActionsInput<precision> toDeviceInput(const ParActionsInput<precision>& in,
    AccDeviceWorkspace<precision>& ws) {
  ParActionsInput<precision> dev=in;

  ws.inputdata.resize(in.dataSize);
  ws.inputdata.copyToDeviceAsync(in.inputdata, accQueue);
  dev.inputdata=ws.inputdata.devicePtr();

  // Pbc only holds fixed-size members, so a bitwise copy is valid on the device
  ws.pbc.resize(1);
  ws.pbc.copyToDeviceAsync(in.pbc, accQueue);
  dev.pbc=ws.pbc.devicePtr();

  if(in.nargs>0) {
    const std::size_t* src[]= {in.ranks, in.shapestarts, in.shapedata, in.ncols,
                               in.bookstarts, in.booksizes, in.bookeeping, in.argstarts
                              };
    const std::size_t sizes[]= {in.ranks_size, in.shapestarts_size, in.shapedata_size, in.ncols_size,
                                in.bookstarts_size, in.booksizes_size, in.bookeeping_size, in.argstarts_size
                               };
    const std::size_t** dst[]= {&dev.ranks, &dev.shapestarts, &dev.shapedata, &dev.ncols,
                                &dev.bookstarts, &dev.booksizes, &dev.bookeeping, &dev.argstarts
                               };
    constexpr unsigned narrays=sizeof(sizes)/sizeof(sizes[0]);
    std::size_t total=0;
    for(unsigned i=0; i<narrays; ++i) {
      total+=sizes[i];
    }
    ws.argbook_host.resize(total);
    std::size_t offset=0;
    for(unsigned i=0; i<narrays; ++i) {
      std::copy(src[i], src[i]+sizes[i], ws.argbook_host.begin()+offset);
      offset+=sizes[i];
    }
    ws.argbook.resize(total);
    ws.argbook.copyToDeviceAsync(ws.argbook_host.data(), accQueue);
    offset=0;
    for(unsigned i=0; i<narrays; ++i) {
      *dst[i]=ws.argbook.devicePtr()+offset;
      offset+=sizes[i];
    }
  }
  return dev;
}

//use the __PLUMED_USE_OPENACC_TASKSMINE macro to debug the ptm ins a single file
//so that compiling witha a small modification will be faster (the ptm is included nearly everywhere)
#ifndef __PLUMED_USE_OPENACC_TASKSMINE
template <class T, typename precision>
void runAllTasksACC(typename T::input_type actiondata,
                    ParActionsInput<precision> myinput,
                    std::vector<precision>& value_stash,
                    const std::vector<unsigned> & partialTaskList,
                    const unsigned nactive_tasks,
                    const std::size_t nderivatives_per_task,
                    const std::size_t workspace_size,
                    AccDeviceWorkspace<precision>& ws
                   ) {
  // from here on the pointers in myinput are device addresses
  myinput = toDeviceInput(myinput, ws);
  auto actiondata_acc = OpenACC::fromToDataHelper(actiondata);

  // the device buffers are reused between steps: resize() reallocates only if they grow
  ws.value_stash.resize(value_stash.size());
  ws.value_stash.copyToDeviceAsync(value_stash.data(), accQueue);
  auto value_stash_data = ws.value_stash.devicePtr();

  ws.taskList.resize(partialTaskList.size());
  ws.taskList.copyToDeviceAsync(partialTaskList.data(), accQueue);
  auto partialTaskList_data = ws.taskList.devicePtr();

  ws.buffer.resize(workspace_size*nactive_tasks);
  auto buffer = ws.buffer.devicePtr();
  ws.derivatives.resize(nderivatives_per_task*nactive_tasks);
  auto derivatives = ws.derivatives.devicePtr();
  // myinput and the scalars are passed by value as kernel arguments
#pragma acc parallel loop async(accQueue) present(actiondata) \
                     firstprivate(myinput, \
                                  nactive_tasks, \
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
    auto myout = ParActionsOutput<precision>::create (myinput.nscalars,
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
  ws.value_stash.copyFromDeviceAsync(value_stash.data(), accQueue);
  // the only synchronization of the call
  acc_wait(accQueue);
}
#else
template <class T, typename precision>
void runAllTasksACC(typename T::input_type actiondata,
                    ParActionsInput<precision> myinput,
                    std::vector<precision>& value_stash,
                    const std::vector<unsigned> & partialTaskList,
                    const unsigned nactive_tasks,
                    const std::size_t nderivatives_per_task,
                    const std::size_t workspace_size,
                    AccDeviceWorkspace<precision>& ws
                   ) ;
#endif //__PLUMED_USE_OPENACC_TASKSMINE

template <class T>
void AccParallelTaskManager<T>::runAllTasks() {
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
  if (comm.Get_rank()== 0) {// no multigpu shenanigans until this works
    runAllTasksACC<T>(
      actiondata,
      myinput,
      value_stash,
      partialTaskList,
      nactive_tasks,
      nderivatives_per_task,
      workspace_size,
      devws
    );
  }
  comm.Bcast( value_stash.data(), value_stash.size(), 0);

  // MPI Gather everything
  if( !serial ) {
    this->comm.Sum( this->value_stash );
  }

// And transfer the value to the output values
  action->transferStashToValues( partialTaskList,  value_stash );
}

//use the __PLUMED_USE_OPENACC_FORCESMINE macro to debug the ptm ins a single file
//so that compiling witha a small modification will be faster (the ptm is included nearly everywhere)
#ifndef __PLUMED_USE_OPENACC_FORCESMINE
template <class T, typename precision>
void applyForcesWithACC(PLMD::View<precision> forcesForApply,
                        typename T::input_type actiondata,
                        ParActionsInput<precision> myinput,
                        const std::vector<precision>& value_stash,
                        const std::vector<unsigned> & partialTaskList,
                        const unsigned nactive_tasks,
                        const std::size_t nderivatives_per_task,
                        const std::size_t workspace_size,
                        AccDeviceWorkspace<precision>& ws
                       ) {
  // from here on the pointers in myinput are device addresses
  myinput = toDeviceInput(myinput, ws);
  auto actiondata_acc = OpenACC::fromToDataHelper(actiondata);

  // the device buffers are reused between steps: resize() reallocates only if they grow
  ws.value_stash.resize(value_stash.size());
  ws.value_stash.copyToDeviceAsync(value_stash.data(), accQueue);
  auto value_stash_data = ws.value_stash.devicePtr();

  ws.taskList.resize(partialTaskList.size());
  ws.taskList.copyToDeviceAsync(partialTaskList.data(), accQueue);
  auto partialTaskList_data = ws.taskList.devicePtr();

  // forcesForApply has been zeroed on the host: this also clears the device buffer
  ws.forces.resize(forcesForApply.size());
  ws.forces.copyToDeviceAsync(forcesForApply.data(), accQueue);
  auto forcesForApply_data = ws.forces.devicePtr();
  const auto forcesForApply_size = ws.forces.size();
  const std::size_t nind_per_scalar = ForceIndexHolder::indexesPerScalar(myinput);

  // Forces on entries that several tasks can touch (the "thread-unsafe" ones,
  // with index >= threadunsafe_forces_start) are summed as follows:
  // - the last nshared entries of the force array (all the thread-unsafe ones
  //   if they are at most maxShared) are summed by each task in its own column
  //   of a buffer, then a second kernel sums the columns. There are no atomics
  //   and no search: this is the case of multicolvars, where the shared
  //   entries are the virial and possibly an origin atom, and it also covers
  //   the virial of matrices and the atoms of a small last group (e.g. GROUPB
  //   of a contact matrix with a single atom).
  // - the other thread-unsafe entries are added with atomics. Within a task,
  //   the contributions of consecutive values that go to the same entry
  //   (e.g. the central atom of a row of a matrix) are first summed in a
  //   small private accumulator, so that only one atomic per entry and task
  //   is done instead of one per value.
  // With atomics the order of the sums is not fixed, so results are
  // reproducible only up to rounding.
  constexpr std::size_t maxShared = 32;
  const std::size_t n_unsafe = forcesForApply_size - myinput.threadunsafe_forces_start;
  const std::size_t nshared = (n_unsafe < maxShared) ? n_unsafe : maxShared;
  const std::size_t sharedStart = forcesForApply_size - nshared;
  // the forces are added as soon as the indices are known, so each task needs
  // space for the indices of one value at a time
  const std::size_t ind_stride_task = nind_per_scalar;

  ws.derivatives.resize(nderivatives_per_task*nactive_tasks);
  auto derivatives = ws.derivatives.devicePtr();
  ws.indices.resize(ind_stride_task*nactive_tasks);
  auto indices = ws.indices.devicePtr();
  ws.valstmp.resize(myinput.sizeOfFakeVals()*nactive_tasks);
  auto valstmp = ws.valstmp.devicePtr();
  ws.buffer.resize(workspace_size*nactive_tasks);
  auto buffer = ws.buffer.devicePtr();
  ws.sharedForces.resize(nshared*nactive_tasks);
  auto sharedForces = ws.sharedForces.devicePtr();

#define derivativeDrift(taskID,scalarID)  taskID*nderivatives_per_task \
                           + scalarID*myinput.ncomponents*myinput.nderivatives_per_scalar
#define stashDrift(taskID,scalarID) taskID*myinput.nscalars \
                           + scalarID*myinput.ncomponents

  // myinput and the scalars are passed by value as kernel arguments.
  // One task per thread (gang vector), and every inner loop is sequential:
  // without this, nvc++ makes the task loop gang-only and vectorizes the
  // short inner loops, i.e. one block of 128 threads per task, mostly idle
#pragma acc parallel loop gang vector async(accQueue) present(actiondata) \
                     firstprivate(myinput, \
                                  nactive_tasks, \
                                  forcesForApply_size, \
                                  nderivatives_per_task, ind_stride_task, \
                                  nshared, sharedStart, \
                                  workspace_size) \
                        deviceptr(derivatives, \
                                  indices, \
                                  value_stash_data, \
                                  partialTaskList_data, \
                                  forcesForApply_data, \
                                  valstmp, \
                                  sharedForces, \
                                  buffer) \
                          default(none)
  for(unsigned t=0; t<nactive_tasks; ++t) {
    std::size_t task_index = partialTaskList_data[t];
    auto myout = ParActionsOutput<precision>::create( myinput.nscalars,
                 valstmp+myinput.nscalars*t,
                 nderivatives_per_task,
                 derivatives+nderivatives_per_task*t,
                 workspace_size,
                 (workspace_size>0)?buffer+workspace_size*t:nullptr);
    // Calculate the stuff in the loop for this action
    T::performTask( task_index, actiondata, myinput, myout );

    // this task's column of the buffer of the shared entries
#pragma acc loop seq
    for(std::size_t k=0; k<nshared; ++k) {
      sharedForces[k*nactive_tasks+t]=0.0;
    }
    // private accumulators of the other thread-unsafe entries: slot s holds
    // the s-th thread-unsafe contribution of the current value; it is summed
    // to the one of the next value if they go to the same entry, and added
    // to the forces with an atomic otherwise
    constexpr unsigned nslots=32;
    constexpr std::size_t noSlot=~std::size_t(0);
    std::size_t slotIndex[nslots];
    precision slotForce[nslots];
#pragma acc loop seq
    for(unsigned s=0; s<nslots; ++s) {
      slotIndex[s]=noSlot;
      slotForce[s]=0.0;
    }

    // If this is a matrix this returns a number that isn't one as we have to loop over the columns
    const std::size_t nvpt = T::getNumberOfValuesPerTask( task_index, actiondata );
#pragma acc loop seq
    for(unsigned vID=0; vID<nvpt; ++vID) {
      auto force_indices = ForceIndexHolder::create(myinput, indices + t*ind_stride_task);
      // Get the indices for forces
      T::getForceIndices( task_index,
                          vID,
                          forcesForApply_size,
                          actiondata,
                          myinput,
                          force_indices );

      // Create a force input object
      auto finput = ForcesInput<precision>::create ( myinput.nscalars,
                    value_stash_data + stashDrift(task_index,vID),
                    myinput.nderivatives_per_scalar,
                    derivatives + derivativeDrift(t,vID));

      // Gather forces that can be gathered locally
      ParallelTaskManager<T>::gatherThreadSafeForces( myinput,
          force_indices,
          finput,
          View<precision>(forcesForApply_data,
                          forcesForApply_size));
      // Add the forces that other tasks can also touch
      unsigned s=0;
#pragma acc loop seq
      for(unsigned i=0; i<myinput.ncomponents; ++i) {
        const precision ff = finput.force[i];
#pragma acc loop seq
        for(std::size_t d=force_indices.threadsafe_derivatives_end[i];
            d<force_indices.tot_indices[i]; ++d, ++s) {
          const std::size_t fi = force_indices.indices[i][d];
          const precision contribution = ff*finput.deriv[i][d];
          if(fi>=sharedStart) {
            sharedForces[(fi-sharedStart)*nactive_tasks+t] += contribution;
          } else if(s<nslots) {
            if(slotIndex[s]==fi) {
              slotForce[s]+=contribution;
            } else {
              if(slotIndex[s]!=noSlot) {
                const precision old=slotForce[s];
#pragma acc atomic update
                forcesForApply_data[slotIndex[s]] += old;
              }
              slotIndex[s]=fi;
              slotForce[s]=contribution;
            }
          } else {
#pragma acc atomic update
            forcesForApply_data[fi] += contribution;
          }
        }
      }
    }
    // flush the private accumulators
#pragma acc loop seq
    for(unsigned s=0; s<nslots; ++s) {
      if(slotIndex[s]!=noSlot) {
        const precision old=slotForce[s];
#pragma acc atomic update
        forcesForApply_data[slotIndex[s]] += old;
      }
    }
  }

  // sum the columns of the buffer of the shared entries: one gang per entry,
  // the threads of the gang sum over the tasks. It runs after the previous
  // kernel, on the same queue, and only it writes these entries
#pragma acc parallel loop gang async(accQueue) \
                     firstprivate(nactive_tasks, nshared, sharedStart) \
                        deviceptr(sharedForces, forcesForApply_data) \
                          default(none)
  for(std::size_t k=0; k<nshared; ++k) {
    double tmp = 0.0;
#pragma acc loop vector reduction(+:tmp)
    for(unsigned t=0; t<nactive_tasks; ++t) {
      tmp += sharedForces[k*nactive_tasks+t];
    }
    forcesForApply_data[sharedStart+k] += tmp;
  }
#undef derivativeDrift
#undef stashDrift
  ws.forces.copyFromDeviceAsync(forcesForApply.data(), accQueue);
  // the only synchronization of the call
  acc_wait(accQueue);
}
#else
template <class T, typename precision>
void applyForcesWithACC(PLMD::View<precision> forcesForApply,
                        typename T::input_type actiondata,
                        ParActionsInput<precision> myinput,
                        const std::vector<precision>& value_stash,
                        const std::vector<unsigned> & partialTaskList,
                        const unsigned nactive_tasks,
                        const std::size_t nderivatives_per_task,
                        const std::size_t workspace_size,
                        AccDeviceWorkspace<precision>& ws
                       );
#endif //__PLUMED_USE_OPENACC_FORCESMINE

template <class T>
void AccParallelTaskManager<T>::applyForces( std::vector<double>& forcesForApply ) {
  // Get the list of active tasks
  std::vector<unsigned> & partialTaskList= action->getListOfActiveTasks( action );
  unsigned nactive_tasks=partialTaskList.size();
  forceData<precision> forces(forcesForApply);
  // Clear force buffer
  std::fill (forces.ffa.begin(),forces.ffa.end(), precision(0.0));
  // Get all the input data so we can broadcast it to the GPU
  myinput.noderiv = false;
  // Retrieve the forces from the values
  action->transferForcesToStash( partialTaskList, value_stash );

  std::fill (omp_forces[0].begin(),omp_forces[0].end(), 0.0);
  if (comm.Get_rank() == 0) {
    applyForcesWithACC<T>(
      PLMD::View<precision> { forces.ffa.data(), forces.ffa.size() },
      actiondata,
      myinput,
      value_stash,
      partialTaskList,
      nactive_tasks,
      nderivatives_per_task,
      workspace_size,
      devws
    );
  }
  forces.update();
  // MPI Gather everything (this must be extended to the gpu thing, after makning it mpi-aware)
  if( !serial ) {
    this->comm.Sum( forcesForApply );
  }
}

} // namespace PLMD
#endif // __PLUMED_acc_ParallelTaskManager_h


