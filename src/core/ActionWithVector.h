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
/// Check if there is a mask value
  int nmask;
/// Is the calculation to be done in serial
  bool serial;
/// Are we in the forward pass through the calculation
  bool forwardPass;
/// The buffer that we use (we keep a copy here to avoid resizing)
  std::vector<double> buffer;
/// A tempory vector of MultiValue so we can avoid doing lots of resizes
  std::vector<MultiValue> myvals;
/// A tempory set of vectors for holding forces over threads
  std::vector<std::vector<double> > omp_forces;
/// The list of active tasks
  std::vector<unsigned> active_tasks;
/// Clear all the bookeeping arrays
  void clearMatrixBookeeping();
///  Run the task
  void runTask( const unsigned& taskno, MultiValue& myvals ) const ;
/// Gather the values that we intend to store in the buffer
  void gatherAccumulators( const unsigned& taskCode, const MultiValue& myvals, std::vector<double>& buffer ) const ;
/// Get the number of stored values in the stream
  bool getNumberOfStoredValues( Value* startat, unsigned& nvals, const unsigned& astart, const std::vector<Value*>& stopat );
/// Check the chain for non scalar forces
  bool checkChainForNonScalarForces() const ;
protected:
/// Turn off the flag that says this action has a masked input
  void ignoreMaskArguments();
/// Run all calculations in serial
  bool runInSerial() const ;
/// Run all the tasks in the list
  void runAllTasks();
/// Check if a force has been added on one of the components of this action
  bool checkComponentsForForce() const ;
/// Accumulate the forces from the Values
  bool checkForForces();
/// Gather the forces on a particular value
  void gatherForcesOnStoredValue( const unsigned& ival, const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const ;
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
/// This is for resizing the task list
  virtual unsigned getNumberOfAtomsPerTask() const { return 0; }
/// Turn off the calculation of the derivatives during the forward pass through a calculation
  bool doNotCalculateDerivatives() const override ;
/// Get the list of tasks that are active
  virtual std::vector<unsigned>& getListOfActiveTasks( ActionWithVector* action );
/// Do we always need to do all the tasks for this action
  virtual void areAllTasksRequired( std::vector<ActionWithVector*>& task_reducing_actions ) {}
/// Find out how many tasks we need to perform in this loop
  virtual void getNumberOfTasks( unsigned& ntasks );
/// Check the status of the ith task
  virtual int checkTaskStatus( const unsigned& taskno, int& flag ) const { return flag; }
/// Determine if a particular task is active based on the values of the input argument
  virtual int checkTaskIsActive( const unsigned& itask ) const ;
/// This we override to perform each individual task
  virtual void performTask( const unsigned& current, MultiValue& myvals ) const = 0;
/// This is used to ensure that all indices are updated when you do local average
  virtual void updateAdditionalIndices( const unsigned& ostrn, MultiValue& myvals ) const {}
/// Can be used to reduce the number of tasks that are performed when you use an ation from elsewhere
  virtual void switchTaskReduction( const bool& task_reduction, ActionWithVector* aselect ) {}
/// Gather the values that we intend to store in the buffer
  virtual void gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals, const unsigned& bufstart, std::vector<double>& buffer ) const {}
/// Check if there is a force that needs to be accumulated on the ith task
  virtual bool checkForTaskForce( const unsigned& itask, const Value* myval ) const ;
/// Gather the forces on non-scalar quantities
  virtual void gatherForces( const unsigned& i, const MultiValue& myvals, std::vector<double>& forces ) const ;
/// This is to transfer data from the buffer to the final value
  void finishComputations( const std::vector<double>& buf );
/// Get the number of forces to use
  virtual void getNumberOfForceDerivatives( unsigned& nforces, unsigned& nderiv ) const ;
/// Apply the forces on this data
  virtual void apply();
};

inline
bool ActionWithVector::runInSerial() const {
  return serial;
}

inline
int ActionWithVector::getNumberOfMasks() const {
  return nmask;
}

inline
void ActionWithVector::ignoreMaskArguments() {
  plumed_assert( nmask<=0 ); nmask=-1;
}

}

#endif
