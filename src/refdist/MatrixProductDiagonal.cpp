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
#include <limits>

//+PLUMEDOC FUNCTION MATRIX_PRODUCT_DIAGONAL
/*
Calculate the product of two matrices and return a vector that contains the diagonal elements of the ouptut vector

This action allows you to [multiply](https://en.wikipedia.org/wiki/Matrix_multiplication) two matrices and recover the diagonal of the
product of the matrices. The following input shows an example where two contact matrices are multiplied together.

```plumed
c1: CONTACT_MATRIX GROUP=1-100 SWITCH={RATIONAL R_0=0.1}
c2: CONTACT_MATRIX GROUP=1-100 SWITCH={RATIONAL R_0=0.2}
m: MATRIX_PRODUCT_DIAGONAL ARG=c1,c2
PRINT ARG=m FILE=colvar
```

As you can see from the documentation for the shortcut [EUCLIDEAN_DISTANCE](EUCLIDEAN_DISTANCE.md) this action is useful for computing distances between
vectors.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace refdist {

class MatrixProductDiagonalInput {
public:
  MatrixProductDiagonalInput& operator=( const MatrixProductDiagonalInput& m ) {
    // Don't need to copy A and B here as they will be set in calculate
    return *this;
  }
};

class MatrixProductDiagonal : public ActionWithVector {
public:
  using input_type = MatrixProductDiagonalInput;
  using PTM = ParallelTaskManager<MatrixProductDiagonal>;
private:
  PTM taskmanager;
public:
  static void registerKeywords( Keywords& keys );
  explicit MatrixProductDiagonal(const ActionOptions&);
  unsigned getNumberOfDerivatives() override ;
  void calculate() override ;
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
  static void performTask( std::size_t task_index, const MatrixProductDiagonalInput& actiondata, ParallelActionsInput& input, ParallelActionsOutput& output );
  static int getNumberOfValuesPerTask( std::size_t task_index, const MatrixProductDiagonalInput& actiondata );
  static void getForceIndices( std::size_t task_index, std::size_t colno, std::size_t ntotal_force, const MatrixProductDiagonalInput& actiondata, const ParallelActionsInput& input, ForceIndexHolder force_indices );
};

PLUMED_REGISTER_ACTION(MatrixProductDiagonal,"MATRIX_PRODUCT_DIAGONAL")

void MatrixProductDiagonal::registerKeywords( Keywords& keys ) {
  ActionWithVector::registerKeywords(keys);
  keys.addInputKeyword("compulsory","ARG","vector/matrix","the two vectors/matrices whose product are to be taken");
  keys.setValueDescription("scalar/vector","a vector containing the diagonal elements of the matrix that obtaned by multiplying the two input matrices together");
  PTM::registerKeywords( keys );
}

MatrixProductDiagonal::MatrixProductDiagonal(const ActionOptions&ao):
  Action(ao),
  ActionWithVector(ao),
  taskmanager(this) {
  if( getNumberOfArguments()!=2 ) {
    error("should be two vectors or matrices in argument to this action");
  }

  const unsigned ncols = [&]()->unsigned{
    if( getPntrToArgument(0)->getRank()==1 ) {
      if( getPntrToArgument(0)->hasDerivatives() ) {
        error("first argument to this action should be a vector or matrix");
      }
      return 1;
    } else if( getPntrToArgument(0)->getRank()==2 ) {
      if( getPntrToArgument(0)->hasDerivatives() ) {
        error("first argument to this action should be a matrix");
      }
      return getPntrToArgument(0)->getShape()[1];
    }
    //you should not be here!
    error("first argument to this action should be a vector or matrix");
    //this is a dummy returns(in fact the code won't get here.
    //this lambda is there due to a -Wmaybe-uninitialized
    return std::numeric_limits<unsigned>::max();
  }();

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
    taskmanager.setupParallelTaskManager( ncols + getPntrToArgument(1)->getShape()[0], 0 );
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
    taskmanager.runAllTasks();
  }
}

void MatrixProductDiagonal::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

void MatrixProductDiagonal::performTask( std::size_t task_index,
    const MatrixProductDiagonalInput& actiondata,
    ParallelActionsInput& input,
    ParallelActionsOutput& output ) {
  auto arg0=ArgumentBookeepingHolder::create( 0, input );
  auto arg1=ArgumentBookeepingHolder::create( 1, input );
  std::size_t fpos = task_index*(1+arg0.ncols);
  std::size_t nmult = arg0.bookeeping[fpos];
  std::size_t vstart = task_index*arg0.ncols;
  output.values[0] = 0;
  if( arg1.ncols<arg1.shape[1] ) {
    plumed_merror("multiplying by a sparse matrix is not implemented as I don't think it is needed");
  } else {
    std::size_t base = arg1.start + task_index;
    for(unsigned i=0; i<nmult; ++i) {
      double val1 = input.inputdata[vstart+i];
      double val2 = input.inputdata[ base + arg1.ncols*arg0.bookeeping[fpos+1+i] ];
      output.values[0] += val1*val2;
      output.derivatives[i] = val2;
      output.derivatives[nmult + i] = val1;
    }
  }
}

int MatrixProductDiagonal::getNumberOfValuesPerTask( std::size_t task_index, const MatrixProductDiagonalInput& actiondata ) {
  return 1;
}

void MatrixProductDiagonal::getForceIndices( std::size_t task_index,
    std::size_t colno,
    std::size_t ntotal_force,
    const MatrixProductDiagonalInput& actiondata,
    const ParallelActionsInput& input,
    ForceIndexHolder force_indices ) {
  auto arg0=ArgumentBookeepingHolder::create( 0, input );
  auto arg1=ArgumentBookeepingHolder::create( 1, input );
  std::size_t fpos = task_index*(1+arg0.ncols);
  std::size_t nmult = arg0.bookeeping[fpos];
  std::size_t vstart = task_index*arg0.ncols;
  if( arg1.ncols<arg1.shape[1] ) {
    plumed_merror("multiplying by a sparse matrix is not implemented as I don't think it is needed");
  } else {
    std::size_t base = arg1.start + task_index;
    for(unsigned i=0; i<nmult; ++i) {
      force_indices.indices[0][i] = vstart + i;
      force_indices.indices[0][nmult+i] = base + arg1.ncols*arg0.bookeeping[fpos+1+i];
    }
    force_indices.threadsafe_derivatives_end[0] = 2*nmult;
    force_indices.tot_indices[0] = 2*nmult;
  }
}

}
}
