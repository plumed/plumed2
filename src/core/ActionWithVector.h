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
/// Work backwards through the chain to find an action that has either stored arguments or derivatives
  const ActionWithVector* getActionWithDerivatives() const ;
/// Check if there are any grids in the stream 
  bool checkForGrids(unsigned& nder) const ;
/// Find out how many tasks we need to perform in this loop
  void getNumberOfTasks( unsigned& ntasks );
///  Run the task
  void runTask( const unsigned& taskno, MultiValue& myvals ) const ;
/// Gather the values that we intend to store in the buffer
  void gatherAccumulators( const unsigned& taskCode, const MultiValue& myvals, std::vector<double>& buffer ) const ;
/// Get the size of the buffer array that holds the data we are gathering over the MPI loop
  void getSizeOfBuffer( const unsigned& nactive_tasks, unsigned& bufsize );
/// Get the number of quantities in the stream
  void getNumberOfStreamedQuantities( unsigned& nquants, unsigned& nmat, unsigned& maxcol, unsigned& nbookeeping );
/// Get the number of derivatives in the stream
  void getNumberOfStreamedDerivatives( unsigned& nderivatives, const std::string& stopat );
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
protected:
/// A vector that contains the start point for the argument derivatives
  std::vector<unsigned> arg_deriv_starts;
/// Assert if this action is part of a chain
  bool done_in_chain;
/// Run all calculations in serial
  bool runInSerial() const ;
/// Get the list of tasks that are active
  std::vector<unsigned>& getListOfActiveTasks();
/// This sets up the arguments at the start of the calculation
  unsigned buildArgumentStore( const unsigned& argstart );
/// Get the position of the argument in the streamm and set it if we need to
  unsigned getArgumentPositionInStream( const unsigned& jder, MultiValue& myvals ) const ;
/// Run all the tasks in the list
  void runAllTasks();  
/// Accumulate the forces from the Values
  bool checkForForces();
public:
  static void registerKeywords( Keywords& keys );
  explicit ActionWithVector(const ActionOptions&);
  virtual ~ActionWithVector() {}
  void lockRequests() override;
  void unlockRequests() override;
  void retrieveAtoms() override;
  void calculateNumericalDerivatives(ActionWithValue* av) override; 
/// Are we running this command in a chain
  bool actionInChain() const ;
/// This is overwritten within ActionWithMatrix and is used to build the chain of just matrix actions
  virtual void finishChainBuild( ActionWithVector* act ) {}
/// Check if there are any stored values in arguments
  bool hasStoredArguments() const ;
/// Return a pointer to the first action in the chain
  const ActionWithVector* getFirstActionInChain() const ;
  ActionWithVector* getFirstActionInChain();
/// Get every the label of every value that is calculated in this chain
  void getAllActionLabelsInChain( std::vector<std::string>& mylabels ) const ;
/// We override clearDerivatives here to prevent data in streams from being deleted
  void clearDerivatives( const bool& force=false ) override;
/// Check if we can be after another ActionWithVector
  virtual bool canBeAfterInChain( ActionWithVector* av ) { return true; }
/// setup the streamed quantities
  virtual void setupStreamedComponents( unsigned& nquants, unsigned& nmat, unsigned& maxcol, unsigned& nbookeeping );
/// This we override to perform each individual task
  virtual void performTask( const unsigned& current, MultiValue& myvals ) const = 0;
/// Gather the data from all the OpenMP threads
  virtual void gatherThreads( const unsigned& nt, const unsigned& bufsize, const std::vector<double>& omp_buffer, std::vector<double>& buffer, MultiValue& myvals );
/// Gather all the data from the MPI processes
  virtual void gatherProcesses( std::vector<double>& buffer );
/// Gather the values that we intend to store in the buffer
  virtual void gatherStoredValue( const unsigned& valindex, const unsigned& code, const MultiValue& myvals, const unsigned& bufstart, std::vector<double>& buffer ) const ;
/// Check if there is a force that needs to be accumulated on the ith task
  virtual bool checkForTaskForce( const unsigned& itask, const Value* myval ) const ;
/// Gather the forces on non-scalar quantities
  virtual void gatherForces( const unsigned& i, const MultiValue& myvals, std::vector<double>& forces ) const ;
/// This is to transfer data from the buffer to the final value
  void finishComputations( const std::vector<double>& buf );
/// Apply the forces on this data
  void apply() override;
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
