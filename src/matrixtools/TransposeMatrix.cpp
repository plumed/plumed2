/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2017 The plumed team
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
#include "ActionWithInputMatrices.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace matrixtools {

class TransposeMatrix : public ActionWithInputMatrices {
private:
  std::vector<std::string> aLabelsInChain;
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit TransposeMatrix(const ActionOptions&);
///
  unsigned getNumberOfColumns() const override;
///
  void getTasksForParent( const std::string& parent, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags );
/// Do the calculation
  void completeMatrixOperations() override;
///
  void apply() override;
///
  double getForceOnMatrixElement( const unsigned& imat, const unsigned& jrow, const unsigned& kcol ) const override;
};

PLUMED_REGISTER_ACTION(TransposeMatrix,"TRANSPOSE")

void TransposeMatrix::registerKeywords( Keywords& keys ) {
  ActionWithInputMatrices::registerKeywords( keys );
}

TransposeMatrix::TransposeMatrix(const ActionOptions& ao):
  Action(ao),
  ActionWithInputMatrices(ao)
{
  if( getNumberOfArguments()!=1 ) error("should only be one argument for this action");
  if( getPntrToArgument(0)->isSymmetric() ) error("input matrix is symmetric.  Transposing will achieve nothing!");
  std::vector<unsigned> shape(2); 
  shape[0]=getPntrToArgument(0)->getShape()[1]; 
  shape[1]=getPntrToArgument(0)->getShape()[0];
  addValue( shape ); 
}

void TransposeMatrix::getTasksForParent( const std::string& parent, std::vector<std::string>& actionsThatSelectTasks, std::vector<unsigned>& tflags ) { 
  if( aLabelsInChain.size()==0 ) getAllActionLabelsInChain( aLabelsInChain );
  bool ignore = checkUsedOutsideOfChain( aLabelsInChain, parent, actionsThatSelectTasks, tflags );
}

unsigned TransposeMatrix::getNumberOfColumns() const {
  return getPntrToArgument(0)->getShape()[0];
}

void TransposeMatrix::completeMatrixOperations() {
  // Retrieve the non-zero pairs
  unsigned nedge=0; retrieveEdgeList( 0, nedge ); std::vector<unsigned> shape( getPntrToOutput(0)->getShape() );
  for(unsigned i=0; i<nedge;++i) getPntrToOutput(0)->set( pairs[i].second*shape[1] + pairs[i].first, vals[i] ); 
}

void TransposeMatrix::apply() {
  // Apply force on the matrix
  if( getPntrToOutput(0)->forcesWereAdded() ) applyForceOnMatrix(0);
}

double TransposeMatrix::getForceOnMatrixElement( const unsigned& imat, const unsigned& jrow, const unsigned& kcol ) const {
  return getPntrToOutput(0)->getForce(kcol*getPntrToOutput(0)->getShape()[1]+jrow);
}


}
}
