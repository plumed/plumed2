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
#include "ActionWithMatrix.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace adjmat {

class MatrixTimesMatrix : public ActionWithMatrix {
private:
  unsigned nderivatives;
public:
  static void registerKeywords( Keywords& keys );
  explicit MatrixTimesMatrix(const ActionOptions&);
  unsigned getNumberOfDerivatives();
  unsigned getNumberOfColumns() const override { return getConstPntrToComponent(0)->getShape()[1]; }
  void setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const ;
  void performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const override;
  void runEndOfRowJobs( const unsigned& ival, const std::vector<unsigned> & indices, MultiValue& myvals ) const override {}
};

PLUMED_REGISTER_ACTION(MatrixTimesMatrix,"MATRIX_PRODUCT")

void MatrixTimesMatrix::registerKeywords( Keywords& keys ) {
  ActionWithMatrix::registerKeywords(keys); keys.use("ARG");
}

MatrixTimesMatrix::MatrixTimesMatrix(const ActionOptions&ao):
Action(ao),
ActionWithMatrix(ao)
{
  if( getNumberOfArguments()!=2 ) error("should be two arguments to this action, a matrix and a vector");
  if( getPntrToArgument(0)->getRank()!=2 || getPntrToArgument(0)->hasDerivatives() ) error("first argument to this action should be a matrix");
  if( getPntrToArgument(1)->getRank()!=2 || getPntrToArgument(1)->hasDerivatives() ) error("first argument to this action should be a matrix");
  if( getPntrToArgument(0)->getShape()[1]!=getPntrToArgument(1)->getShape()[0] ) error("number of columns in first matrix does not equal number of rows in second matrix");
  std::vector<unsigned> shape(2); shape[0]=getPntrToArgument(0)->getShape()[0]; shape[1]=getPntrToArgument(1)->getShape()[1];
  addValue( shape ); setNotPeriodic(); nderivatives = buildArgumentStore(0);
}

unsigned MatrixTimesMatrix::getNumberOfDerivatives() {
  return nderivatives;
}

void MatrixTimesMatrix::setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const {
  unsigned start_n = getPntrToArgument(0)->getShape()[0], size_v = getPntrToArgument(1)->getShape()[1]; 
  if( indices.size()!=size_v+1 ) indices.resize( size_v+1 );
  for(unsigned i=0; i<size_v; ++i) indices[i+1] = start_n + i;
  myvals.setSplitIndex( size_v + 1 );
}

void MatrixTimesMatrix::performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const {
  unsigned ostrn = getConstPntrToComponent(0)->getPositionInStream(), ind2=index2;
  if( index2>=getPntrToArgument(0)->getShape()[0] ) ind2 = index2 - getPntrToArgument(0)->getShape()[0];

  unsigned nmult=getPntrToArgument(0)->getShape()[1]; double matval=0;
  for(unsigned i=0; i<nmult; ++i) {
      double val1 = getElementOfMatrixArgument( 0, index1, i, myvals );
      double val2 = getElementOfMatrixArgument( 1, i, ind2, myvals ); 
      matval+= val1*val2;

      if( doNotCalculateDerivatives() ) continue;

      addDerivativeOnMatrixArgument( 0, 0, index1, i, val2, myvals );
      addDerivativeOnMatrixArgument( 0, 1, i, ind2, val1, myvals ); 
  }
  // And add this part of the product
  myvals.addValue( ostrn, matval );
}

}
}
