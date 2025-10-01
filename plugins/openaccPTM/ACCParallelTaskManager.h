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

namespace PLMD {

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
                    const std::size_t workspace_size
                   ) {
  auto myinput_acc = OpenACC::fromToDataHelper(myinput);
  auto actiondata_acc = OpenACC::fromToDataHelper(actiondata);

  //template type is deduced
  OpenACC::memoryManager vs{value_stash};
  auto value_stash_data = vs.devicePtr();

  OpenACC::memoryManager ptl{partialTaskList};
  auto partialTaskList_data = ptl.devicePtr();

  OpenACC::memoryManager<precision> buff{workspace_size*nactive_tasks};

  auto buffer = buff.devicePtr();
  OpenACC::memoryManager<precision> dev(nderivatives_per_task*nactive_tasks);
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
  vs.copyFromDevice(value_stash.data());
}
#else
template <class T, typename precision>
void runAllTasksACC(typename T::input_type actiondata,
                    ParActionsInput<precision> myinput,
                    std::vector<precision>& value_stash,
                    const std::vector<unsigned> & partialTaskList,
                    const unsigned nactive_tasks,
                    const std::size_t nderivatives_per_task,
                    const std::size_t workspace_size
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
      workspace_size
    );
  }
  comm.Bcast( value_stash.data(), value_stash.size(), 0);

  // MPI Gather everything
  if( !serial ) {
    this->comm.Sum( this->value_stash );
  }

// And transfer the value to the output values
  action->transferStashToValues( value_stash );
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

  OpenACC::memoryManager<precision> dev{nderivatives_per_task*nactive_tasks};
  auto derivatives = dev.devicePtr();
  OpenACC::memoryManager<std::size_t> ind{nind_per_task*nactive_tasks};
  auto indices = ind.devicePtr();
  OpenACC::memoryManager<precision> vtmp{myinput.sizeOfFakeVals()*nactive_tasks};
  auto valstmp = vtmp.devicePtr();
  OpenACC::memoryManager<precision> buff{workspace_size*nactive_tasks};
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
      auto myout = ParActionsOutput<precision>::create( myinput.nscalars,
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

          auto fdata = ForcesInput<precision>::create( myinput.nscalars,
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
template <class T, typename precision>
void applyForcesWithACC(PLMD::View<precision> forcesForApply,
                        typename T::input_type actiondata,
                        ParActionsInput<precision> myinput,
                        const std::vector<precision>& value_stash,
                        const std::vector<unsigned> & partialTaskList,
                        const unsigned nactive_tasks,
                        const std::size_t nderivatives_per_task,
                        const std::size_t workspace_size
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
  action->transferForcesToStash( value_stash );

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
      workspace_size
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


