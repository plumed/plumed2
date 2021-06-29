/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2017 The plumed team
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
#ifndef __PLUMED_adjmat_MatrixProductBase_h
#define __PLUMED_adjmat_MatrixProductBase_h

#include "core/ActionAtomistic.h"
#include "core/ActionWithArguments.h"
#include "core/ActionWithValue.h"

namespace PLMD {
namespace adjmat {

class MatrixProductBase :
  public ActionAtomistic,
  public ActionWithArguments,
  public ActionWithValue
{
private:
  bool skip_ieqj;
  void updateCentralMatrixIndex( const unsigned& ind, MultiValue& myvals ) const ;
protected:
  std::vector<double> forcesToApply;
public:
  static void registerKeywords( Keywords& keys );
  explicit MatrixProductBase(const ActionOptions&);
  unsigned getNumberOfDerivatives() const override;
  virtual unsigned getNumberOfColumns() const ;
  bool mustBeTreatedAsDistinctArguments() const override ;
  void getTasksForParent( const std::string& parent, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags );
  void lockRequests() override;
  void unlockRequests() override;
  void calculateNumericalDerivatives( ActionWithValue* a=NULL ) override;
  void calculate() override;
  void update() override;
  void runFinalJobs() override;
  unsigned getNumberOfFinalTasks() override;
  std::vector<unsigned> getMatrixShapeForFinalTasks() override;
  virtual void performTask( const unsigned& task_index, MultiValue& myvals ) const ;
  virtual bool performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const ;
  virtual double computeVectorProduct( const unsigned& index1, const unsigned& index2,
                                       const std::vector<double>& vec1, const std::vector<double>& vec2,
                                       std::vector<double>& dvec1, std::vector<double>& dvec2, MultiValue& myvals ) const = 0;
  void apply() override;
};

}
}
#endif

