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
#include "ACCDevice.h"

#include <algorithm>
#include <cstring>
#include <vector>

namespace PLMD {

/// The async queue used by the task manager. All the transfers and the kernels
/// of a call are queued on it and the host waits only once, at the end.
inline constexpr int accQueue=1;

/// A device buffer that remembers what was last copied into it.
/// update() queues the copy only if the data differs (byte by byte) from the
/// previous one, so data that does not change between calls (the argument
/// bookkeeping and the task list, which runAllTasks and applyForces of the
/// same step share and which are usually the same at every step, and the box
/// in NVT) is transferred only once.
/// The copy is made from the host mirror, which is not modified until the
/// next call: the caller does not need to keep the data alive until the
/// queue is waited.
/// The mirror costs a host copy of the data when it changes, so this is
/// meant for small arrays or arrays that rarely change.
template <typename T>
class AccCachedBuffer {
  OpenACC::memoryManager<T> device;
  std::vector<unsigned char> host;
  bool valid=false;
public:
  T* update(const T* data, std::size_t n) {
    const std::size_t nbytes=n*sizeof(T);
    const bool same = valid && host.size()==nbytes
                      && (nbytes==0 || std::memcmp(host.data(), data, nbytes)==0);
    if(!same) {
      host.resize(nbytes);
      if(nbytes>0) {
        std::memcpy(host.data(), data, nbytes);
      }
      device.resize(n);
      device.copyToDeviceAsync(reinterpret_cast<const T*>(host.data()), accQueue);
      valid=true;
    }
    return device.devicePtr();
  }
};

/// Sets n elements of a device array to zero, on accQueue
template <typename T>
void accZeroAsync(T* p, const std::size_t n) {
  if(n==0) {
    return;
  }
#pragma acc parallel loop gang vector async(accQueue) firstprivate(n) deviceptr(p) default(none)
  for(std::size_t i=0; i<n; ++i) {
    p[i]=T(0);
  }
}

/// Device buffers used by the OpenACC task manager.
/// They are kept between calls (one set per action) and grown when needed, so
/// that the device memory is not allocated and freed at every step.
/// Allocation happens on first use, so ranks that do not run on the GPU do not allocate.
template <typename precision>
struct AccDeviceWorkspace {
  /// device copies of the data pointed by ParActionsInput (see toDeviceInput)
  OpenACC::memoryManager<precision> inputdata;
  /// the host array inputdata was last copied from, and its size
  const precision* inputdata_src=nullptr;
  std::size_t inputdata_size=0;
  AccCachedBuffer<Pbc> pbc;
  AccCachedBuffer<std::size_t> argbook;
  /// host staging area for the argument bookkeeping, packed in a single array
  std::vector<std::size_t> argbook_host;
  AccCachedBuffer<unsigned> taskList;
  OpenACC::memoryManager<precision> value_stash;
  OpenACC::memoryManager<precision> buffer;
  OpenACC::memoryManager<precision> derivatives;
  /// address and size of derivatives when it was last cleared
  const precision* derivatives_cleared=nullptr;
  std::size_t derivatives_cleared_size=0;
  OpenACC::memoryManager<precision> forces;
  OpenACC::memoryManager<std::size_t> indices;
  OpenACC::memoryManager<precision> valstmp;
  /// per-task partial sums of the forces on the last entries of the force array
  OpenACC::memoryManager<precision> sharedForces;
  /// Resizes the buffer of the derivatives and returns its device address.
  /// The CPU task manager gives each thread a zeroed buffer at every call, so
  /// the entries that an action never writes (and that it may still list in
  /// getForceIndices) are zero; device memory is not initialized, and garbage
  /// there (even multiplied by a zero force, if it is a NaN) corrupts the
  /// forces. The buffer is therefore cleared whenever it is reallocated or its
  /// size changes: the entries that are written are rewritten at every call,
  /// the others stay zero.
  precision* getDerivatives(const std::size_t n) {
    derivatives.resize(n);
    if(derivatives.devicePtr()!=derivatives_cleared || n!=derivatives_cleared_size) {
      accZeroAsync(derivatives.devicePtr(), n);
      derivatives_cleared=derivatives.devicePtr();
      derivatives_cleared_size=n;
    }
    return derivatives.devicePtr();
  }
};

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
  using ParallelTaskManager<T>::input_buffer;
  using ParallelTaskManager<T>::useacc;
  using ParallelTaskManager<T>::getValueStashSize;
  using ParallelTaskManager<T>::serial;
public:
  typedef typename ParallelTaskManager<T>::ParallelActionsInput ParallelActionsInput;
  typedef typename ParallelTaskManager<T>::ParallelActionsOutput ParallelActionsOutput;
  typedef typename ParallelTaskManager<T>::ForceInput ForceInput;
  typedef typename ParallelTaskManager<T>::ForceOutput ForceOutput;
  typedef typename ParallelTaskManager<T>::precision precision;
private:
/// The GPU of this action (see ACCDevice.h)
  ACC::Device device;
/// Device buffers kept between steps, freed on their GPU
  std::unique_ptr<AccDeviceWorkspace<precision>> devws;
/// The active tasks of this process, when they are shared between the processes
  std::vector<unsigned> my_tasks;
/// MULTIGPU: with MPI the tasks are shared between the processes, each using its GPU
  bool multigpu=false;
/// A copy of actiondata that stays on the device between the calls (see residentActionData)
  std::unique_ptr<typename T::input_type> resident_actiondata;
/// The generation of actiondata (see ParallelTaskManager::getActionInputGeneration) of resident_actiondata
  unsigned long resident_generation=0;
/// Returns a copy of actiondata that is present on the device (of the current device scope).
/// It is copied and uploaded only when actiondata may have changed since the last upload:
/// before, actiondata was copied on the host (it was passed by value) and moved with enter/exit
/// data at every call, twice per step. For large action data (e.g. the 370k segments of the
/// secondary structure variables of 144 peptides: 45 MB) this was most of the time of the action
  const typename T::input_type* residentActionData() {
    const unsigned long generation=this->getActionInputGeneration();
    if( !resident_actiondata || generation!=resident_generation ) {
      if( resident_actiondata ) {
        resident_actiondata->removeFromACCDevice();
      }
      resident_actiondata=std::make_unique<typename T::input_type>( actiondata );
      resident_actiondata->toACCDevice();
      resident_generation=generation;
    }
    return resident_actiondata.get();
  }
/// The GPU of this action (chosen at the first call), -1 without GPUs
  int getDevice() {
    return device.resolve(*action);
  }
/// True if the tasks are shared between the MPI processes, each using its GPU
/// (MULTIGPU, without SERIAL). Otherwise the first process computes all of them
/// and sends the results to the others: by default, because each call has costs
/// on the CPU (copies of the input, of the results and of the forces) that do
/// not decrease with the number of processes, while the sums over the processes
/// add more (see the README of the plugin)
  bool splitTasks() const {
    return multigpu && !serial && comm.Get_size()>1;
  }
/// The active tasks computed by this process, a contiguous block of the list:
/// [first,first+count). Contiguous blocks keep the part of the value stash
/// written by a process contiguous, so that only that part is transferred
  void getMyTasks(const unsigned ntasks, unsigned& first, unsigned& count) const {
    if(!splitTasks()) {
      first=0;
      count=(comm.Get_rank()==0) ? ntasks : 0;
      return;
    }
    const unsigned long long size=comm.Get_size(), rank=comm.Get_rank();
    first=static_cast<unsigned>((ntasks*rank)/size);
    count=static_cast<unsigned>((ntasks*(rank+1))/size)-first;
  }
public:
  static void registerKeywords( Keywords& keys ) {
    ParallelTaskManager<T>::registerKeywords(keys);
    keys.addFlag("MULTIGPU",false,"with MPI, share the tasks between the processes, each using its own GPU, and sum the results "
                 "(by default only the first process computes, on its GPU, and sends the results to the others)");
  }
  AccParallelTaskManager(ActionWithVector* av):
    PLMD::ParallelTaskManager<T>(av),
    devws(std::make_unique<AccDeviceWorkspace<precision>>()) {
    useacc=true;
    action->log.printf("  using GPU to calculate this action\n");
    action->parseFlag("MULTIGPU",multigpu);
    if( multigpu && !serial ) {
      action->log.printf("  with MPI the tasks are shared between the processes, each using its GPU\n");
    } else if( multigpu ) {
      action->log.printf("  MULTIGPU is ignored with SERIAL\n");
    }
  }
  ~AccParallelTaskManager() {
    // members are destroyed after the scope ends: the buffers are freed here
    ACC::DeviceScope scope(device.id());
    if( resident_actiondata ) {
      resident_actiondata->removeFromACCDevice();
    }
    devws.reset();
  }
  AccParallelTaskManager(const AccParallelTaskManager&)=delete;
  AccParallelTaskManager& operator=(const AccParallelTaskManager&)=delete;
/// This runs all the tasks
  void runAllTasks();
/// Apply the forces on the parallel object
  void applyForces( std::vector<double>& forcesForApply, bool reset_forces=true );
};

struct ACCPTM {
  template <typename T>
  using PTM=AccParallelTaskManager<T>;
};

/// Returns a copy of the input whose pointers are device addresses.
/// The data (input data, Pbc and argument bookkeeping) is copied in the
/// device buffers of the workspace, with one transfer per array (the
/// bookkeeping arrays are packed in a single buffer). The Pbc and the
/// bookkeeping are copied only if they changed since the previous call.
/// The input data (positions or arguments, which change at every step) is
/// copied when called from runAllTasks (refresh=true). In applyForces
/// (refresh=false) it is copied only if it was not uploaded from the same
/// host array: the task manager fills input_buffer only in runAllTasks,
/// which always uploads it, so the device copy is then still current.
/// The returned object can be passed to the kernels as firstprivate, i.e. as
/// a kernel argument.
/// The copies are queued on accQueue: in.inputdata must not change until the
/// queue is waited.
/// Before this, the input was put on the device with "enter data" and removed
/// with "exit data" at every call: each pointer member costs a transfer and
/// an attach/detach (another small transfer), i.e. about 40 transfers of a
/// few microseconds each per call, which dominated the time of small actions.
template <typename precision>
ParActionsInput<precision> toDeviceInput(const ParActionsInput<precision>& in,
    AccDeviceWorkspace<precision>& ws, const bool refresh) {
  ParActionsInput<precision> dev=in;

  if(refresh || ws.inputdata_src!=in.inputdata || ws.inputdata_size!=in.dataSize) {
    ws.inputdata.resize(in.dataSize);
    ws.inputdata.copyToDeviceAsync(in.inputdata, accQueue);
    ws.inputdata_src=in.inputdata;
    ws.inputdata_size=in.dataSize;
  }
  dev.inputdata=ws.inputdata.devicePtr();

  // Pbc only holds fixed-size members, so a bitwise copy is valid on the device
  dev.pbc=ws.pbc.update(in.pbc, 1);

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
    const std::size_t* argbook=ws.argbook.update(ws.argbook_host.data(), total);
    offset=0;
    for(unsigned i=0; i<narrays; ++i) {
      *dst[i]=argbook+offset;
      offset+=sizes[i];
    }
  }
  return dev;
}

//use the __PLUMED_USE_OPENACC_TASKSMINE macro to debug the ptm ins a single file
//so that compiling witha a small modification will be faster (the ptm is included nearly everywhere)
#ifndef __PLUMED_USE_OPENACC_TASKSMINE
/// Runs the tasks tasks[0..nactive_tasks) (nactive_tasks>0) on the device.
/// Task t writes the entries [t*nscalars,(t+1)*nscalars) of the value stash:
/// the part of the stash between the first and the last task is cleared on
/// the device (tasks may leave some entries untouched) and copied back, the
/// rest of the host stash is set to zero (it belongs to the tasks of the
/// other processes, or to inactive tasks)
template <class T, typename precision>
void runAllTasksACC(const typename T::input_type* ad,
                    ParActionsInput<precision> myinput,
                    std::vector<precision>& value_stash,
                    const unsigned* tasks,
                    const unsigned nactive_tasks,
                    const std::size_t nderivatives_per_task,
                    const std::size_t workspace_size,
                    AccDeviceWorkspace<precision>& ws
                   ) {
  const auto [tmin,tmax] = std::minmax_element(tasks,tasks+nactive_tasks);
  const std::size_t stash_lo = std::min(std::size_t(*tmin)*myinput.nscalars, value_stash.size());
  const std::size_t stash_hi = std::min((std::size_t(*tmax)+1)*myinput.nscalars, value_stash.size());
  // from here on the pointers in myinput are device addresses
  myinput = toDeviceInput(myinput, ws, true);

  // the device buffers are reused between steps: resize() reallocates only if they grow
  ws.value_stash.resize(value_stash.size());
  auto value_stash_data = ws.value_stash.devicePtr();
  accZeroAsync(value_stash_data+stash_lo, stash_hi-stash_lo);

  auto partialTaskList_data = ws.taskList.update(tasks, nactive_tasks);

  ws.buffer.resize(workspace_size*nactive_tasks);
  auto buffer = ws.buffer.devicePtr();
  auto derivatives = ws.getDerivatives(nderivatives_per_task*nactive_tasks);
  // myinput and the scalars are passed by value as kernel arguments.
  // One task per thread, as in the force kernel
  {
#pragma acc parallel loop gang vector async(accQueue) present(ad[0:1]) \
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
      T::performTask( task_index, *ad, myinput, myout );
    }
  }
  if(stash_hi>stash_lo) {
    acc_memcpy_from_device_async(value_stash.data()+stash_lo, value_stash_data+stash_lo,
                                 (stash_hi-stash_lo)*sizeof(precision), accQueue);
  }
  std::fill(value_stash.begin(), value_stash.begin()+stash_lo, precision(0));
  std::fill(value_stash.begin()+stash_hi, value_stash.end(), precision(0));
  // the only synchronization of the call
  acc_wait(accQueue);
}
#else
template <class T, typename precision>
void runAllTasksACC(const typename T::input_type* ad,
                    ParActionsInput<precision> myinput,
                    std::vector<precision>& value_stash,
                    const unsigned* tasks,
                    const unsigned nactive_tasks,
                    const std::size_t nderivatives_per_task,
                    const std::size_t workspace_size,
                    AccDeviceWorkspace<precision>& ws
                   ) ;
#endif //__PLUMED_USE_OPENACC_TASKSMINE

template <class T>
void AccParallelTaskManager<T>::runAllTasks() {
  // Get the list of active tasks
  std::vector<unsigned> & partialTaskList( action->getListOfActiveTasks() );
  const unsigned nactive_tasks=partialTaskList.size();
  // Reset the values at the start of the task loop
  std::size_t totalvals=getValueStashSize();
  if( value_stash.size()!=totalvals ) {
    value_stash.resize(totalvals);
  }
  // With MULTIGPU each process computes a block of the tasks on its GPU and
  // the stashes are summed; otherwise the first process computes all the
  // tasks, and the stash is sent to the others
  const bool split = splitTasks();
  unsigned first, count;
  getMyTasks( nactive_tasks, first, count );
  if( split || comm.Get_rank()==0 ) {
    // The derivatives are only ever needed during the force calculation
    // Please stop adding a check on whether any of the methods need the derivatives
    // they do not need the derivatives during the forward loop.
    myinput.noderiv = true;
    // Get all the input data so we can broadcast it to the GPU
    action->getInputData( input_buffer, argumentsMap );
    myinput.dataSize = input_buffer.size();
    myinput.inputdata = input_buffer.data();
    // Transfer all the bookeeping information about the arguments
    // (argumentsMap was filled by getInputData)
    myinput.setupArguments( argumentsMap );
    if( count>0 ) {
      ACC::DeviceScope scope(getDevice());
      runAllTasksACC<T>(
        residentActionData(),
        myinput,
        value_stash,
        partialTaskList.data()+first,
        count,
        nderivatives_per_task,
        workspace_size,
        *devws
      );
    } else {
      std::fill( value_stash.begin(), value_stash.end(), precision(0) );
      // the input data was not uploaded: applyForces must upload it if it has tasks
      // (it can have some, when it works on the columns of a matrix)
      devws->inputdata_src=nullptr;
    }
  }
  if( split ) {
    comm.Sum( value_stash );
  } else if( comm.Get_size()>1 ) {
    // the first process holds the complete result: it must be broadcast, not summed
    // (summing would multiply the values by the number of processes)
    comm.Bcast( value_stash.data(), value_stash.size(), 0);
  }

// And transfer the value to the output values
  action->transferStashToValues( partialTaskList,  value_stash );
}

//use the __PLUMED_USE_OPENACC_FORCESMINE macro to debug the ptm ins a single file
//so that compiling witha a small modification will be faster (the ptm is included nearly everywhere)
#ifndef __PLUMED_USE_OPENACC_FORCESMINE
/// Adds to forcesForApply the forces of the tasks tasks[0..nactive_tasks)
/// (nactive_tasks>0). Task t reads the entries [t*nforcescalars,(t+1)*nforcescalars)
/// of the stash of the forces on the values: only the part of the stash between
/// the first and the last task is copied to the device
template <class T, typename precision>
void applyForcesWithACC(PLMD::View<precision> forcesForApply,
                        const typename T::input_type* ad,
                        ParActionsInput<precision> myinput,
                        const std::vector<precision>& value_stash,
                        const unsigned* tasks,
                        const unsigned nactive_tasks,
                        const std::size_t nderivatives_per_task,
                        const std::size_t workspace_size,
                        AccDeviceWorkspace<precision>& ws,
                        const bool reset_forces
                       ) {
  const auto [tmin,tmax] = std::minmax_element(tasks,tasks+nactive_tasks);
  const std::size_t stash_lo = std::min(std::size_t(*tmin)*myinput.nforcescalars, value_stash.size());
  const std::size_t stash_hi = std::min((std::size_t(*tmax)+1)*myinput.nforcescalars, value_stash.size());
  // from here on the pointers in myinput are device addresses.
  // The input data is the one of the last runAllTasks, so nothing is copied
  // here unless something changed in between
  myinput = toDeviceInput(myinput, ws, false);

  // the device buffers are reused between steps: resize() reallocates only if they grow
  ws.value_stash.resize(value_stash.size());
  auto value_stash_data = ws.value_stash.devicePtr();
  if(stash_hi>stash_lo) {
    acc_memcpy_to_device_async(value_stash_data+stash_lo, const_cast<precision*>(value_stash.data())+stash_lo,
                               (stash_hi-stash_lo)*sizeof(precision), accQueue);
  }

  auto partialTaskList_data = ws.taskList.update(tasks, nactive_tasks);

  // the forces are accumulated on the device and copied back as a whole.
  // They start from zero, or from the forces already in forcesForApply when
  // the caller is accumulating over several calls (reset_forces=false, used
  // when the forces on the columns of a matrix are gathered separately)
  ws.forces.resize(forcesForApply.size());
  auto forcesForApply_data = ws.forces.devicePtr();
  const auto forcesForApply_size = ws.forces.size();
  if( reset_forces ) {
    accZeroAsync(forcesForApply_data, forcesForApply_size);
  } else {
    ws.forces.copyToDeviceAsync(forcesForApply.data(), accQueue);
  }
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
  const std::size_t n_unsafe = (forcesForApply_size > myinput.threadunsafe_forces_start) ?
                               forcesForApply_size - myinput.threadunsafe_forces_start : 0;
  const std::size_t nshared = (n_unsafe < maxShared) ? n_unsafe : maxShared;
  const std::size_t sharedStart = forcesForApply_size - nshared;
  // the forces are added as soon as the indices are known, so each task needs
  // space for the indices of one value at a time
  const std::size_t ind_stride_task = nind_per_scalar;

  auto derivatives = ws.getDerivatives(nderivatives_per_task*nactive_tasks);
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
#define stashDrift(taskID,scalarID) taskID*myinput.nforcescalars \
                           + scalarID*myinput.ncomponents

  {
    // myinput and the scalars are passed by value as kernel arguments.
    // One task per thread (gang vector), and every inner loop is sequential:
    // without this, nvc++ makes the task loop gang-only and vectorizes the
    // short inner loops, i.e. one block of 128 threads per task, mostly idle
#pragma acc parallel loop gang vector async(accQueue) present(ad[0:1]) \
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
      T::performTask( task_index, *ad, myinput, myout );

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
      const std::size_t nvpt = T::getNumberOfValuesPerTask( task_index, *ad );
#pragma acc loop seq
      for(unsigned vID=0; vID<nvpt; ++vID) {
        auto force_indices = ForceIndexHolder::create(myinput, indices + t*ind_stride_task);
        // Get the indices for forces
        T::getForceIndices( task_index,
                            vID,
                            forcesForApply_size,
                            *ad,
                            myinput,
                            force_indices );

        // Create a force input object
        auto finput = ForcesInput<precision>::create ( myinput.nforcescalars,
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
    if(nshared>0) {
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
    }
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
                        const typename T::input_type* ad,
                        ParActionsInput<precision> myinput,
                        const std::vector<precision>& value_stash,
                        const unsigned* tasks,
                        const unsigned nactive_tasks,
                        const std::size_t nderivatives_per_task,
                        const std::size_t workspace_size,
                        AccDeviceWorkspace<precision>& ws,
                        const bool reset_forces
                       );
#endif //__PLUMED_USE_OPENACC_FORCESMINE

template <class T>
void AccParallelTaskManager<T>::applyForces( std::vector<double>& forcesForApply, bool reset_forces ) {
  // Get the list of active tasks
  std::vector<unsigned> & partialTaskList= action->getListOfActiveTasks();
  const unsigned nactive_tasks=partialTaskList.size();
  // As in runAllTasks: with MULTIGPU each process computes the forces of a block
  // of the tasks on its GPU and the forces are summed; otherwise the first
  // process computes all of them and sends them to the others
  const bool split = splitTasks();
  unsigned first, count;
  getMyTasks( nactive_tasks, first, count );
  // With reset_forces=false the forces already in forcesForApply are kept: the
  // caller is accumulating over several calls (the forces on the columns of a
  // matrix are gathered separately). They are the same on all the processes,
  // so only the first one keeps them, and the sum counts them once
  const bool from_zero = reset_forces || (split && comm.Get_rank()!=0);
  if( count>0 ) {
    ACC::DeviceScope scope(getDevice());
    // the force buffer is cleared (or not) on the device
    forceData<precision> forces(forcesForApply);
    myinput.noderiv = false;
    // Retrieve the forces from the values (only those of the tasks of this process are used)
    if( split ) {
      my_tasks.assign( partialTaskList.begin()+first, partialTaskList.begin()+first+count );
      action->transferForcesToStash( my_tasks, value_stash );
    } else {
      action->transferForcesToStash( partialTaskList, value_stash );
    }
    applyForcesWithACC<T>(
      PLMD::View<precision> { forces.ffa.data(), forces.ffa.size() },
      residentActionData(),
      myinput,
      value_stash,
      partialTaskList.data()+first,
      count,
      nderivatives_per_task,
      workspace_size,
      *devws,
      from_zero
    );
    forces.update();
  } else if( (split || comm.Get_rank()==0) && from_zero ) {
    std::fill( forcesForApply.begin(), forcesForApply.end(), 0.0 );
  }
  // Every process needs the complete forces, also with SERIAL (where the CPU
  // version computes them on every process)
  if( split ) {
    comm.Sum( forcesForApply );
  } else if( comm.Get_size()>1 ) {
    comm.Bcast( forcesForApply, 0 );
  }
}

} // namespace PLMD
#endif // __PLUMED_acc_ParallelTaskManager_h
