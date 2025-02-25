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
#ifndef __PLUMED_core_ActionWithMatrix_h
#define __PLUMED_core_ActionWithMatrix_h

#include "ActionWithVector.h"

namespace PLMD {

class ActionWithMatrix : public ActionWithVector {
private:
/// This is used to clear up the matrix elements
  void clearMatrixElements( MultiValue& myvals ) const ;
/// This does the calculation of a particular matrix element
  void runTask( const std::string& controller, const unsigned& current, const unsigned colno, MultiValue& myvals ) const ;
public:
  static void registerKeywords( Keywords& keys );
  explicit ActionWithMatrix(const ActionOptions&);
///
  virtual bool isAdjacencyMatrix() const {
    return false;
  }
/// This should return the number of columns to help with sparse storage of matrices
  virtual unsigned getNumberOfColumns() const = 0;
//// This does some setup before we run over the row of the matrix
  virtual void setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const = 0;
/// Run over one row of the matrix
  virtual void performTask( const unsigned& task_index, MultiValue& myvals ) const ;
/// This is the virtual that will do the calculation of the task for a particular matrix element
  virtual void performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const = 0;
/// This is the jobs that need to be done when we have run all the jobs in a row of the matrix
  virtual void runEndOfRowJobs( const unsigned& ival, const std::vector<unsigned> & indices, MultiValue& myvals ) const = 0;
/// This is overwritten in Adjacency matrix where we have a neighbour list
  virtual void updateNeighbourList() {}
/// Run the calculation
  virtual void calculate() override;
/// Check if there are forces we need to account for on this task
  bool checkForTaskForce( const unsigned& itask, const Value* myval ) const override ;
/// This gathers the force on a particular value
  void gatherForces( const unsigned& itask, const MultiValue& myvals, std::vector<double>& forces ) const override ;
};

}
#endif
