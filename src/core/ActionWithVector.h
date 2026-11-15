/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#ifndef __PLUMED_core_ActionWithVector_h
#define __PLUMED_core_ActionWithVector_h

#include "ActionWithValue.h"
#include "ActionAtomistic.h"
#include "ActionWithArguments.h"
#include <limits>
#include <vector>

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
};

class ActionWithVector:
  public ActionAtomistic,
  public ActionWithValue,
  public ActionWithArguments {
  friend class Value;
public:
/// How the rows of the matrices are shared between the MPI processes (see getRowDistribution)
  struct RowDistribution {
/// The number of actions when this was worked out: it is worked out again if actions are added
    std::size_t nactions=0;
/// Incremented every time this is worked out: the buffers that are cleared row by row are then cleared completely
    unsigned generation=0;
/// The tasks are shared by the row of the matrices they work on: row r is done by process r%size
    bool byRow=false;
/// Each row of the output matrices is only used by the process that computes it: the values are not summed over the processes
    bool localOutput=false;
/// The forces on these arguments are only needed by the process that owns each row: they are not summed over the processes
    std::vector<bool> localArgument;
  };
private:
  mutable RowDistribution rowdist;
/// The sizes of the arguments in the forces when they were last cleared, and the generation of rowdist then
  std::vector<std::size_t> clearedForceSizes;
  unsigned clearedForceGeneration=0;
  static constexpr int NoMasksUsed =-1;
/// Check if there is a mask value
  int nmask=NoMasksUsed;
/// The list of active tasks
  std::vector<unsigned> active_tasks;
protected:
/// Turn off the flag that says this action has a masked input
  void ignoreMaskArguments();
/// Accumulate the forces from the Values
  bool checkForForces();
/// Copy all the data to the bookeeping array
  void copyArgumentBookeeping( ArgumentsBookkeeping& argumentsMap ) const ;
/// Implementation of transferStashToValues for both precisions
  template <typename T>
  void stashToValuesImpl( const std::vector<unsigned>& partialTaskList, const std::vector<T>& stash );
/// Implementation of transferForcesToStash for both precisions
  template <typename T>
  void forcesToStashImpl( const std::vector<unsigned>& partialTaskList, std::vector<T>& stash ) const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit ActionWithVector(const ActionOptions&);
  virtual ~ActionWithVector() {}
  void lockRequests() override;
  void unlockRequests() override;
  virtual void prepare() override;
/// Check if a mask has been set
  int getNumberOfMasks() const ;
  void calculateNumericalDerivatives(ActionWithValue* av) override;
/// Get the list of tasks that are active
  virtual std::vector<unsigned>& getListOfActiveTasks();
  const std::vector<unsigned>& getConstListOfActiveTasks() const ;
/// Find out how many tasks we need to perform in this loop
  virtual void getNumberOfTasks( unsigned& ntasks );
/// Determine if a particular task is active based on the values of the input argument
  virtual int checkTaskIsActive( const unsigned& itask ) const ;
/// This is so we can parallelize with GPU
  virtual void getInputData( std::vector<double>& inputdata, ArgumentsBookkeeping& argumentsMap ) const ;
  virtual void getInputData( std::vector<float>& inputdata, ArgumentsBookkeeping& argumentsMap ) const ;
/// Get the atomic input data
  virtual void getAtomicInputData( std::vector<double>& inputdata ) const ;
  virtual void getAtomicInputData( std::vector<float>& inputdata ) const ;
/// Get the argument input data
  virtual void getArgumentInputData( std::vector<double>& inputdata ) const ;
  virtual void getArgumentInputData( std::vector<float>& inputdata ) const ;
/// Get the size of the value in the stash
  virtual unsigned getSizeOfValueInStash( const Value* myval );
/// This is so we an transfer data gathered in the parallel task manager to the underlying values
  virtual void transferStashToValues( const std::vector<unsigned>& partialTaskList, const std::vector<double>& stash );
/// This is so we an transfer data gathered in the parallel task manager to the underlying values
  virtual void transferStashToValues( const std::vector<unsigned>& partialTaskList, const std::vector<float>& stash );
/// This is so we can transfer forces from the values to the parallel task manager
  virtual void transferForcesToStash( const std::vector<unsigned>& partialTaskList, std::vector<double>& stash ) const ;
/// This is so we can transfer forces from the values to the parallel task manager
  virtual void transferForcesToStash( const std::vector<unsigned>& partialTaskList, std::vector<float>& stash ) const ;
/// Get the number of forces to use
  unsigned getNumberOfForceDerivatives() const ;
/// Apply the forces on this data
  void apply() override;
/// Apply the forces on non-zero rank objects
  virtual void applyNonZeroRankForces( std::vector<double>& outforces ) {
    plumed_error();
  }
/// The number of consecutive tasks that work on the same row of the matrices (1 when a task is a row)
  virtual std::size_t getNumberOfTasksPerRow() const {
    return 1;
  }
/// True if the rows of the output matrices can be computed each by one MPI process only
  virtual bool canDistributeRows() const {
    return false;
  }
/// True if task t only reads row t/getNumberOfTasksPerRow() of the matrix argument iarg, and writes forces only on that row
  virtual bool readsOnlyOwnRowOfArgument( unsigned iarg ) const {
    return false;
  }
/// Find out, at the first call, which matrices need not be summed over the MPI processes
/// (it is const because it caches a property of the input, which does not change)
  const RowDistribution& getRowDistribution() const ;
/// Sum the forces over the MPI processes, except those on the arguments that are owned by row
  void sumForcesOverProcesses( std::vector<double>& outforces );
/// Set the forces to zero before the backward loop: on the arguments owned by row only the rows of this process are used
  void clearForcesForApply( std::vector<double>& outforces );
/// With localOutput only the rows of this process are cleared
  void clearInputForces( const bool& force=false ) override;
  void clearDerivatives( const bool& force=false ) override;
};

inline
int ActionWithVector::getNumberOfMasks() const {
  return nmask;
}

inline
void ActionWithVector::ignoreMaskArguments() {
  plumed_assert( nmask<=0 );
  nmask=NoMasksUsed;
}

inline
const std::vector<unsigned>& ActionWithVector::getConstListOfActiveTasks() const {
  return active_tasks;
}

inline
unsigned ActionWithVector::getSizeOfValueInStash( const Value* myval ) {
  if( myval->hasDerivatives() ) {
    return myval->getNumberOfStoredValues()*(1+getNumberOfDerivatives());
  }
  return myval->getNumberOfStoredValues();
}

}

#endif
