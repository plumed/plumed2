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

class ActionWithVector:
  public ActionAtomistic,
  public ActionWithValue,
  public ActionWithArguments {
  friend class Value;
private:
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
  virtual std::vector<unsigned>& getListOfActiveTasks( ActionWithVector* action );
/// Find out how many tasks we need to perform in this loop
  virtual void getNumberOfTasks( unsigned& ntasks );
/// Determine if a particular task is active based on the values of the input argument
  virtual int checkTaskIsActive( const unsigned& itask ) const ;
/// This is so we can parallelize with GPU
  virtual void getInputData( std::vector<double>& inputdata ) const ;
/// This is so we an transfer data gathered in the parallel task manager to the underlying values
  virtual void transferStashToValues( const std::vector<unsigned>& partialTaskList, const std::vector<double>& stash );
/// This is so we can transfer forces from the values to the parallel task manager
  virtual void transferForcesToStash( const std::vector<unsigned>& partialTaskList, std::vector<double>& stash ) const ;
/// Get the number of forces to use
  unsigned getNumberOfForceDerivatives() const ;
/// Apply the forces on this data
  void apply() override;
/// Apply the forces on non-zero rank objects
  virtual void applyNonZeroRankForces( std::vector<double>& outforces ) {
    plumed_error();
  }
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

}

#endif
