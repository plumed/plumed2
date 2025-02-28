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
#include "core/ParallelTaskManager.h"
#include "core/ActionRegister.h"
#include "core/MatrixView.h"

//+PLUMEDOC FUNCTION MATRIX_PRODUCT_DIAGONAL
/*
Calculate the product of two matrices and return a vector that contains the diagonal elements of the ouptut vector

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace refdist {

class MatrixProductDiagonalInput {
public:
  MatrixView A;
  MatrixView B;
  MatrixProductDiagonalInput& operator=( const MatrixProductDiagonalInput& m ) { 
    // Don't need to copy A and B here as they will be set in calculate
    return *this;
  }
};

class MatrixProductDiagonal : public ActionWithVector {
private:
  ParallelTaskManager<MatrixProductDiagonal,MatrixProductDiagonalInput> taskmanager;
public:
  static void registerKeywords( Keywords& keys );
  explicit MatrixProductDiagonal(const ActionOptions&);
  unsigned getNumberOfDerivatives() override ;
  void calculate() override ;
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
  void performTask( const unsigned& task_index, MultiValue& myvals ) const override { plumed_error(); }
  static void performTask( std::size_t task_index, const MatrixProductDiagonalInput& actiondata, ParallelActionsInput& input, ParallelActionsOutput& output );
  static void gatherForces( std::size_t task_index, const MatrixProductDiagonalInput& actiondata, const ParallelActionsInput& input, const ForceInput& fdata, ForceOutput& forces );
};

PLUMED_REGISTER_ACTION(MatrixProductDiagonal,"MATRIX_PRODUCT_DIAGONAL")

void MatrixProductDiagonal::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords(keys);
  keys.addInputKeyword("compulsory","ARG","vector/matrix","the two vectors/matrices whose product are to be taken");
  keys.setValueDescription("scalar/vector","a vector containing the diagonal elements of the matrix that obtaned by multiplying the two input matrices together");
  ParallelTaskManager<MatrixProductDiagonal,MatrixProductDiagonalInput>::registerKeywords( keys );
}

MatrixProductDiagonal::MatrixProductDiagonal(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  taskmanager(this) {
  if( getNumberOfArguments()!=2 ) {
    error("should be two vectors or matrices in argument to this action");
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
    std::vector<std::size_t> shape(1);
    shape[0]=getPntrToArgument(0)->getShape()[0];
    addValue( shape );
    setNotPeriodic();
    taskmanager.setupParallelTaskManager( 1, ncols + getPntrToArgument(1)->getShape()[0], 0 );
    taskmanager.setActionInput( MatrixProductDiagonalInput() );
  }
}

unsigned MatrixProductDiagonal::getNumberOfDerivatives() {
  if( doNotCalculateDerivatives() ) {
    return 0;
  }
  return getPntrToArgument(0)->getNumberOfValues() + getPntrToArgument(1)->getNumberOfValues();;
}

void MatrixProductDiagonal::calculate() {
  if( getPntrToArgument(0)->getRank()==1 && getPntrToArgument(1)->getRank()==1 ) {
    double val1 = getPntrToArgument(0)->get(0);
    double val2 = getPntrToArgument(1)->get(0);
    Value* myval = getPntrToComponent(0);
    myval->set( val1*val2 );
  
    if( doNotCalculateDerivatives() ) {
      return;
    }
  
    myval->setDerivative( 0, val2 );
    myval->setDerivative( 1, val1 ); 
  } else if( getPntrToArgument(1)->getRank()==1 ) {
    Value* arg1 = getPntrToArgument(0);
    Value* arg2 = getPntrToArgument(1);
    unsigned nmult = arg1->getRowLength(0);
    unsigned nvals1 = arg1->getNumberOfValues();

    double matval=0;
    Value* myval = getPntrToComponent(0);
    for(unsigned i=0; i<nmult; ++i) {
        unsigned kind = arg1->getRowIndex( 0, i );
        double val1 = arg1->get( i, false );
        double val2 = arg2->get( kind );
        matval += val1*val2;

        if( doNotCalculateDerivatives() ) {
          continue;
        }

        myval->addDerivative( kind, val2 );
        myval->addDerivative( nvals1 + kind, val1 );
    }  
    myval->set( matval );
  } else {
    taskmanager.getActionInput().A.setup( 0, getPntrToArgument(0) );
    taskmanager.getActionInput().B.setup( getPntrToArgument(0)->getNumberOfStoredValues(), getPntrToArgument(1) );
    taskmanager.runAllTasks();
  }
}

void MatrixProductDiagonal::applyNonZeroRankForces( std::vector<double>& outforces ) {
   taskmanager.applyForces( outforces );
}

void MatrixProductDiagonal::performTask( std::size_t task_index, const MatrixProductDiagonalInput& actiondata, ParallelActionsInput& input, ParallelActionsOutput& output ){
   output.values[0] = 0;
   std::size_t ibase = task_index*(actiondata.A.ncols+1);
   std::size_t n = actiondata.A.bookeeping[ibase];
   for(unsigned i=0; i<n; ++i) {
       double val1 = input.inputdata[task_index*actiondata.A.ncols+i];
       double val2 = MatrixView::getElement( actiondata.A.bookeeping[ibase+1+i], task_index, actiondata.B, input.inputdata );
       output.values[0] += val1*val2;
       output.derivatives[i] = val2;
       output.derivatives[n +i] = val1;
   } 
}

void MatrixProductDiagonal::gatherForces( std::size_t task_index, const MatrixProductDiagonalInput& actiondata, const ParallelActionsInput& input, const ForceInput& fdata, ForceOutput& forces ){
   double ff = fdata.force[0];
   std::size_t ibase = task_index*(actiondata.A.ncols+1);
   std::size_t n = actiondata.A.bookeeping[ibase];
   for(unsigned i=0; i<n; ++i) {
       std::size_t ival = actiondata.A.bookeeping[ibase+1+i];
       forces.thread_unsafe[ task_index*actiondata.A.shape[1] + ival ] = ff*fdata.deriv[0][i];
       forces.thread_unsafe[ actiondata.B.start + ival*actiondata.B.shape[1] + task_index ] = ff*fdata.deriv[0][n+i];
   }
}

}
}
