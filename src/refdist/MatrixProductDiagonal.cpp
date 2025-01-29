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
#include "core/ActionWithVector.h"
#include "core/ActionRegister.h"

//+PLUMEDOC FUNCTION MATRIX_PRODUCT_DIAGONAL
/*
Calculate the product of two matrices and return a vector that contains the diagonal elements of the ouptut vector

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace refdist {

class MatrixProductDiagonal : public ActionWithVector {
private:
public:
  static void registerKeywords( Keywords& keys );
  explicit MatrixProductDiagonal(const ActionOptions&);
  unsigned getNumberOfDerivatives() override ;
  void calculate() override ;
  void performTask( const unsigned& task_index, MultiValue& myvals ) const override ;
};

PLUMED_REGISTER_ACTION(MatrixProductDiagonal,"MATRIX_PRODUCT_DIAGONAL")

void MatrixProductDiagonal::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords(keys);
  keys.addInputKeyword("compulsory","ARG","vector/matrix","the two vectors/matrices whose product are to be taken");
  keys.setValueDescription("scalar/vector","a vector containing the diagonal elements of the matrix that obtaned by multiplying the two input matrices together");
}

MatrixProductDiagonal::MatrixProductDiagonal(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao) {
  if( getNumberOfArguments()!=2 ) {
    error("should be two arguments to this action, a matrix and a vector");
  }

  unsigned ncols;
  if( getPntrToArgument(0)->getRank()==1 ) {
    if( getPntrToArgument(0)->hasDerivatives() ) {
      error("first argument to this action should be a vector or matrix");
    }
    ncols = 1;
  } else if( getPntrToArgument(0)->getRank()==2 ) {
    if( getPntrToArgument(0)->hasDerivatives() ) {
      error("first argument to this action should be a matrix");
    }
    ncols = getPntrToArgument(0)->getShape()[1];
  }

  if( getPntrToArgument(1)->getRank()==1 ) {
    if( getPntrToArgument(1)->hasDerivatives() ) {
      error("second argument to this action should be a vector or matrix");
    }
    if( ncols!=getPntrToArgument(1)->getShape()[0] ) {
      error("number of columns in first matrix does not equal number of elements in vector");
    }
    if( getPntrToArgument(0)->getShape()[0]!=1 ) {
      error("matrix output by this action must be square");
    }
    addValueWithDerivatives();
    setNotPeriodic();
  } else {
    if( getPntrToArgument(1)->getRank()!=2 || getPntrToArgument(1)->hasDerivatives() ) {
      error("second argument to this action should be a vector or a matrix");
    }
    if( ncols!=getPntrToArgument(1)->getShape()[0] ) {
      error("number of columns in first matrix does not equal number of rows in second matrix");
    }
    if( getPntrToArgument(0)->getShape()[0]!=getPntrToArgument(1)->getShape()[1] ) {
      error("matrix output by this action must be square");
    }
    std::vector<unsigned> shape(1);
    shape[0]=getPntrToArgument(0)->getShape()[0];
    addValue( shape );
    setNotPeriodic();
  }
  getPntrToArgument(0)->buildDataStore();
  getPntrToArgument(1)->buildDataStore();
}

unsigned MatrixProductDiagonal::getNumberOfDerivatives() {
  if( doNotCalculateDerivatives() ) {
    return 0;
  }
  return getPntrToArgument(0)->getNumberOfValues() + getPntrToArgument(1)->getNumberOfValues();;
}

void MatrixProductDiagonal::performTask( const unsigned& task_index, MultiValue& myvals ) const {
  unsigned ostrn = getConstPntrToComponent(0)->getPositionInStream();
  Value* arg1 = getPntrToArgument(0);
  Value* arg2 = getPntrToArgument(1);
  if( arg1->getRank()==1 ) {
    double val1 = arg1->get( task_index );
    double val2 = arg2->get( task_index );
    myvals.addValue( ostrn, val1*val2 );

    if( doNotCalculateDerivatives() ) {
      return;
    }

    myvals.addDerivative( ostrn, task_index, val2 );
    myvals.updateIndex( ostrn, task_index );
    unsigned nvals = getPntrToArgument(0)->getNumberOfValues();
    myvals.addDerivative( ostrn, nvals + task_index, val1 );
    myvals.updateIndex( ostrn, nvals + task_index );
  } else {
    unsigned nmult = arg1->getRowLength(task_index);
    unsigned nrowsA = getPntrToArgument(0)->getShape()[1];
    unsigned nrowsB = 1;
    if( getPntrToArgument(1)->getRank()>1 ) {
      nrowsB = getPntrToArgument(1)->getShape()[1];
    }
    unsigned nvals1 = getPntrToArgument(0)->getNumberOfValues();

    double matval = 0;
    for(unsigned i=0; i<nmult; ++i) {
      unsigned kind = arg1->getRowIndex( task_index, i );
      double val1 = arg1->get( task_index*nrowsA + kind );
      double val2 = arg2->get( kind*nrowsB + task_index );
      matval += val1*val2;

      if( doNotCalculateDerivatives() ) {
        continue;
      }

      myvals.addDerivative( ostrn, task_index*nrowsA + kind, val2 );
      myvals.updateIndex( ostrn, task_index*nrowsA + kind );
      myvals.addDerivative( ostrn, nvals1 + kind*nrowsB + task_index, val1 );
      myvals.updateIndex( ostrn, nvals1 + kind*nrowsB + task_index );
    }
    // And add this part of the product
    myvals.addValue( ostrn, matval );
  }
}

void MatrixProductDiagonal::calculate() {
  if( getPntrToArgument(1)->getRank()==1 ) {
    unsigned nder = getNumberOfDerivatives();
    MultiValue myvals( 1, nder, 0, 0, 0 );
    performTask( 0, myvals );

    Value* myval=getPntrToComponent(0);
    myval->set( myvals.get(0) );
    for(unsigned i=0; i<nder; ++i) {
      myval->setDerivative( i, myvals.getDerivative(0,i) );
    }
  } else {
    runAllTasks();
  }
}

}
}
