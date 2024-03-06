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
#include "tools/MultiValue.h"
#include <vector>

namespace PLMD {

class ActionWithVector:
  public ActionAtomistic,
  public ActionWithValue,
  public ActionWithArguments
{
  friend class Value;
private:
/// Is the calculation to be done in serial
  bool serial;
/// The buffer that we use (we keep a copy here to avoid resizing)
  std::vector<double> buffer;
/// The list of active tasks
  std::vector<unsigned> active_tasks;
  /// Action that must be done before this one
  ActionWithVector* action_to_do_before;
/// Actions that must be done after this one
  ActionWithVector* action_to_do_after;
/// This is the list of actions that control the tasks that we do here
  std::vector<ActionWithVector*> task_control_list;
/// Work backwards through the chain to find an action that has either stored arguments or derivatives
  ActionWithVector* getActionWithDerivatives( ActionWithVector* depaction );
/// Check if there are any grids in the stream
  bool checkForGrids(unsigned& nder) const ;
///  Run the task
  void runTask( const unsigned& taskno, MultiValue& myvals ) const ;
/// Gather the values that we intend to store in the buffer
  void gatherAccumulators( const unsigned& taskCode, const MultiValue& myvals, std::vector<double>& buffer ) const ;
/// Gather the forces on non-scalar quantities
  void gatherForces( const unsigned& i, const MultiValue& myvals, std::vector<double>& forces ) const ;
/// Get the size of the buffer array that holds the data we are gathering over the MPI loop
  void getSizeOfBuffer( const unsigned& nactive_tasks, unsigned& bufsize );
/// Get the number of quantities in the stream
  void getNumberOfStreamedQuantities( const std::string& headstr, unsigned& nquants, unsigned& nmat, unsigned& maxcol, unsigned& nbookeeping );
/// Get the number of stored values in the stream
  bool getNumberOfStoredValues( Value* startat, unsigned& nvals, const unsigned& astart, const std::vector<Value*>& stopat );
/// Add this action to the recursive chain
  bool addActionToChain( const std::vector<std::string>& alabels, ActionWithVector* act );
/// Check the chain for non scalar forces
  bool checkChainForNonScalarForces() const ;
/// Check if a force has been added on one of the components of this action
  bool checkComponentsForForce() const ;
/// Get the tasks that we need for forces
  void getForceTasks( std::vector<unsigned>& force_tasks ) const ;
/// Add the gathered forces to the inputs across the whole chain
  void addForcesToInput( const std::vector<double>& forcesToApply, unsigned& ind );
/// Check if this ation can reduce the number of tasks we perform
  void canReduceTasks( std::vector<ActionWithVector*>& task_reducing_actions );
/// Send information to arguments that tasks are reduced in depedent actions
  void broadcastThatTasksAreReduced( ActionWithVector* aselect );
/// Turn on task reduction flag in dependent actions
  void updateTaskReductionFlag( bool& head_reduce_tasks );
/// Check if a particular task is active at this time
  void taskIsActive( const unsigned& current, int& flag ) const ;
/// This is turned on if there is some action that needs all the tasks
  bool never_reduce_tasks;
/// Are we allowed to reduce the number of tasks being performed
  bool reduce_tasks;
/// Were the atoms retrieved in some earlier action
  bool atomsWereRetrieved;
/// This is used to build the argument store when we cannot use the chain
  unsigned reallyBuildArgumentStore( const unsigned& argstart );
protected:
/// A vector that contains the start point for the argument derivatives
  std::vector<unsigned> arg_deriv_starts;
/// Assert if this action is part of a chain
  bool done_in_chain;
/// This updates whether or not we are using all the task reduction stuff
  void updateTaskListReductionStatus();
/// Run all calculations in serial
  bool runInSerial() const ;
/// Get the list of tasks that are active
  std::vector<unsigned>& getListOfActiveTasks( ActionWithVector* action );
/// Check if the arguments of this action depend on thearg
  bool argumentDependsOn( const std::string& headstr, ActionWithVector* faction, Value* thearg );
/// This sets up the arguments at the start of the calculation
  unsigned buildArgumentStore( const unsigned& argstart );
/// Run all the tasks in the list
  void runAllTasks();
/// Accumulate the forces from the Values
  bool checkForForces();
public:
  static void registerKeywords( Keywords& keys );
  explicit ActionWithVector(const ActionOptions&);
  virtual ~ActionWithVector();
  void lockRequests() override;
  void unlockRequests() override;
  virtual void prepare() override;
  void retrieveAtoms( const bool& force=false ) override;
  void calculateNumericalDerivatives(ActionWithValue* av) override;
/// Are we running this command in a chain
  bool actionInChain() const ;
/// This is overwritten within ActionWithMatrix and is used to build the chain of just matrix actions
  virtual void finishChainBuild( ActionWithVector* act );
/// Check if there are any stored values in arguments
  bool hasStoredArguments() const ;
/// Return a pointer to the first action in the chain
  const ActionWithVector* getFirstActionInChain() const ;
  ActionWithVector* getFirstActionInChain();
/// This is overridden in ActionWithMatrix
  virtual void getAllActionLabelsInMatrixChain( std::vector<std::string>& matchain ) const {}
/// Get the number of derivatives in the stream
  void getNumberOfStreamedDerivatives( unsigned& nderivatives, Value* stopat );
/// Get every the label of every value that is calculated in this chain
  void getAllActionLabelsInChain( std::vector<std::string>& mylabels ) const ;
/// We override clearInputForces here to ensure that forces are deleted from all values
  void clearInputForces( const bool& force=false ) override;
/// We override clearDerivatives here to prevent data in streams from being deleted
  void clearDerivatives( const bool& force=false ) override;
/// Check if we can be after another ActionWithVector
  virtual bool canBeAfterInChain( ActionWithVector* av ) { return true; }
/// Do we always need to do all the tasks for this action
  virtual void areAllTasksRequired( std::vector<ActionWithVector*>& task_reducing_actions ) {}
/// Find out how many tasks we need to perform in this loop
  virtual void getNumberOfTasks( unsigned& ntasks );
/// Check the status of the ith task
  virtual int checkTaskStatus( const unsigned& taskno, int& flag ) const { return flag; }
/// Check if we are in a subchain
  virtual bool isInSubChain( unsigned& nder ) { return false; }
/// Get the additional tasks that are required here
  virtual void getAdditionalTasksRequired( ActionWithVector* action, std::vector<unsigned>& atasks );
/// setup the streamed quantities
  virtual void setupStreamedComponents( const std::string& headstr, unsigned& nquants, unsigned& nmat, unsigned& maxcol, unsigned& nbookeeping );
/// This we override to perform each individual task
  virtual void performTask( const unsigned& current, MultiValue& myvals ) const = 0;
/// This is used to ensure that all indices are updated when you do local average
  virtual void updateAdditionalIndices( const unsigned& ostrn, MultiValue& myvals ) const {}
/// Gather the data from all the OpenMP threads
  virtual void gatherThreads( const unsigned& nt, const unsigned& bufsize, const std::vector<double>& omp_buffer, std::vector<double>& buffer, MultiValue& myvals );
/// Can be used to reduce the number of tasks that are performed when you use an ation from elsewhere
  virtual void switchTaskReduction( const bool& task_reduction, ActionWithVector* aselect ) {}
/// Gather all the data from the MPI processes
  virtual void gatherProcesses( std::vector<double>& buffer );
/// Gather the values that we intend to store in the buffer
  virtual void gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals, const unsigned& bufstart, std::vector<double>& buffer ) const ;
/// Get the force tasks that are active for this action
  virtual void updateForceTasksFromValue( const Value* myval, std::vector<unsigned>& force_tasks ) const ;
/// Check if there is a force that needs to be accumulated on the ith task
  virtual bool checkForTaskForce( const unsigned& itask, const Value* myval ) const ;
/// Gather the forces on a particular value
  virtual void gatherForcesOnStoredValue( const Value* myval, const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const ;
/// This is to transfer data from the buffer to the final value
  void finishComputations( const std::vector<double>& buf );
/// Apply the forces on this data
  virtual void apply();
};

inline
bool ActionWithVector::actionInChain() const {
  return (action_to_do_before!=NULL);
}

inline
bool ActionWithVector::runInSerial() const {
  return serial;
}

}

#endif
