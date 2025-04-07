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
#ifndef __PLUMED_matrixtools_MatrixTimesMatrix_h
#define __PLUMED_matrixtools_MatrixTimesMatrix_h

#include "core/ActionWithMatrix.h"
#include "core/ParallelTaskManager.h"

namespace PLMD {
namespace matrixtools {

template <class T>
class MatrixTimesMatrixInput {
public:
  T funcinput;
  RequiredMatrixElements outmat;
};

class InputVectors {
public:
  std::size_t nelem;
  std::vector<double> arg1, arg2;
  InputVectors( std::size_t n ) : nelem(n), arg1(n), arg2(n) {}
};

template <class T>
class MatrixTimesMatrix : public ActionWithMatrix {
public:
  using input_type = MatrixTimesMatrixInput<T>;
  using PTM = ParallelTaskManager<MatrixTimesMatrix<T>>;
private:
  PTM taskmanager;
public:
  static void registerKeywords( Keywords& keys );
  explicit MatrixTimesMatrix(const ActionOptions&);
  void prepare() override ;
  unsigned getNumberOfDerivatives() override;
  void calculate() override ;
  void applyNonZeroRankForces( std::vector<double>& outforces ) override ;
  static void performTask( std::size_t task_index, const MatrixTimesMatrixInput<T>& actiondata, ParallelActionsInput& input, ParallelActionsOutput& output );
  static void gatherForces( std::size_t task_index, const MatrixTimesMatrixInput<T>& actiondata, const ParallelActionsInput& input, const ForceInput& fdata, ForceOutput forces );
};

template <class T>
void MatrixTimesMatrix<T>::registerKeywords( Keywords& keys ) {
  ActionWithMatrix::registerKeywords(keys);
  keys.addInputKeyword("optional","MASK","matrix","a matrix that is used to used to determine which elements of the output matrix to compute");
  keys.addInputKeyword("compulsory","ARG","matrix","the label of the two matrices from which the product is calculated");
  keys.addFlag("ELEMENTS_ON_DIAGONAL_ARE_ZERO",false,"set all diagonal elements to zero");
  T::registerKeywords( keys );
  PTM::registerKeywords( keys );
}

template <class T>
MatrixTimesMatrix<T>::MatrixTimesMatrix(const ActionOptions&ao):
  Action(ao),
  ActionWithMatrix(ao),
  taskmanager(this) {
  int nm=getNumberOfMasks();
  if( nm<0 ) {
    nm = 0;
  }
  if( getNumberOfArguments()-nm!=2 ) {
    error("should be two arguments to this action, a matrix and a vector");
  }
  if( getPntrToArgument(0)->getRank()!=2 || getPntrToArgument(0)->hasDerivatives() ) {
    error("first argument to this action should be a matrix");
  }
  if( getPntrToArgument(1)->getRank()!=2 || getPntrToArgument(1)->hasDerivatives() ) {
    error("second argument to this action should be a matrix");
  }
  if( getPntrToArgument(0)->getShape()[1]!=getPntrToArgument(1)->getShape()[0] ) {
    error("number of columns in first matrix does not equal number of rows in second matrix");
  }
  std::vector<std::size_t> shape(2);
  shape[0]=getPntrToArgument(0)->getShape()[0];
  shape[1]=getPntrToArgument(1)->getShape()[1];
  addValue( shape );
  setNotPeriodic();
  getPntrToComponent(0)->reshapeMatrixStore( shape[1] );
  if( getName()!="DISSIMILARITIES" && getPntrToArgument(0)->isDerivativeZeroWhenValueIsZero() && getPntrToArgument(1)->isDerivativeZeroWhenValueIsZero() ) {
    getPntrToComponent(0)->setDerivativeIsZeroWhenValueIsZero();
  }

  if( nm>0 ) {
    unsigned iarg = getNumberOfArguments()-1;
    if( getPntrToArgument(iarg)->getRank()!=2 || getPntrToArgument(0)->hasDerivatives() ) {
      error("argument passed to MASK keyword should be a matrix");
    }
    if( getPntrToArgument(iarg)->getShape()[0]!=shape[0] || getPntrToArgument(iarg)->getShape()[1]!=shape[1] ) {
      error("argument passed to MASK keyword has the wrong shape");
    }
  }
  MatrixTimesMatrixInput<T> actdata;
  actdata.funcinput.setup( this, getPntrToArgument(0) );
  taskmanager.setupParallelTaskManager( 1, 2, getPntrToArgument(0)->getNumberOfStoredValues() );
  taskmanager.setActionInput( actdata );
}

template <class T>
unsigned MatrixTimesMatrix<T>::getNumberOfDerivatives() {
  return getPntrToArgument(0)->getNumberOfStoredValues() + getPntrToArgument(1)->getNumberOfStoredValues();
}

template <class T>
void MatrixTimesMatrix<T>::prepare() {
  ActionWithVector::prepare();
  Value* myval = getPntrToComponent(0);
  if( myval->getShape()[0]==getPntrToArgument(0)->getShape()[0] && myval->getShape()[1]==getPntrToArgument(1)->getShape()[1] ) {
    return;
  }
  std::vector<std::size_t> shape(2);
  shape[0]=getPntrToArgument(0)->getShape()[0];
  shape[1]=getPntrToArgument(1)->getShape()[1];
  myval->setShape(shape);
  myval->reshapeMatrixStore( shape[1] );
}

template <class T>
void MatrixTimesMatrix<T>::calculate() {
  if( !getPntrToComponent(0)->isDerivativeZeroWhenValueIsZero() ) {
    if( getPntrToArgument(0)->getNumberOfColumns()<getPntrToArgument(0)->getShape()[1] ) {
      if( !doNotCalculateDerivatives() ) {
        error("cannot calculate derivatives for this action with sparse matrices");
      } else if( getName()=="DISSIMILARITIES" ) {
        error("cannot calculate dissimilarities for sparse matrices");
      }
    }
    if( getPntrToArgument(1)->getNumberOfColumns()<getPntrToArgument(1)->getShape()[1] ) {
      if( !doNotCalculateDerivatives() ) {
        error("cannot calculate derivatives for this action with sparse matrices");
      } else if( getName()=="DISSIMILARITIES" ) {
        error("cannot calculate dissimilarities for sparse matrices");
      }
    }
  }
  updateBookeepingArrays( taskmanager.getActionInput().outmat );
  unsigned nvals = getPntrToComponent(0)->getNumberOfColumns();
  taskmanager.setupParallelTaskManager( 1, 2*nvals*getPntrToArgument(0)->getNumberOfColumns(), getPntrToArgument(1)->getNumberOfStoredValues(), nvals );
  taskmanager.runAllTasks();
}

template <class T>
void MatrixTimesMatrix<T>::performTask( std::size_t task_index,
                                        const MatrixTimesMatrixInput<T>& actiondata,
                                        ParallelActionsInput& input,
                                        ParallelActionsOutput& output ) {
  ArgumentBookeepingHolder arg0( 0, input ), arg1( 1, input );
  std::size_t fpos = task_index*(1+arg0.ncols);
  std::size_t nmult = arg0.bookeeping[fpos];
  std::size_t vstart = task_index*arg0.ncols;
  InputVectors vectors( nmult );
  if( arg1.ncols<arg1.shape[1] ) {
    std::size_t fstart = task_index*(1+actiondata.outmat.ncols);
    std::size_t nelements = actiondata.outmat.bookeeping[fstart];
    for(unsigned i=0; i<nelements; ++i) {
      std::size_t nm = 0;
      for(unsigned j=0; j<nmult; ++j) {
        std::size_t kind = arg0.bookeeping[fpos+1+j];
        std::size_t bstart = kind*(arg1.ncols + 1);
        std::size_t nr = arg1.bookeeping[bstart];
        for(unsigned k=0; k<nr; ++k) {
          if( arg1.bookeeping[bstart+1+k]==actiondata.outmat.bookeeping[fstart+1+i] ) {
            nm++;
            break;
          }
        }
      }
      vectors.nelem = nm;
      nm = 0;
      for(unsigned j=0; j<nmult; ++j) {
        std::size_t kind = arg0.bookeeping[fpos+1+j];
        std::size_t bstart = kind*(arg1.ncols + 1);
        std::size_t nr = arg1.bookeeping[bstart];
        for(unsigned k=0; k<nr; ++k) {
          if( arg1.bookeeping[bstart+1+k]==actiondata.outmat.bookeeping[fstart+1+i] ) {
            vectors.arg1[nm] = input.inputdata[ vstart + j ];
            vectors.arg2[nm] = input.inputdata[ arg1.start + kind*arg1.ncols + k ];
            nm++;
            break;
          }
        }
      }
      MatrixElementOutput elem( 1, 2*nmult, output.values.data() + i, output.derivatives.data() + 2*nmult*i );
      T::calculate( input.noderiv, actiondata.funcinput, vectors, elem );
    }
  } else {
    // Retrieve the row of the first matrix
    for(unsigned i=0; i<nmult; ++i) {
      vectors.arg1[i] = input.inputdata[ vstart + i ];
    }

    // Now do our multiplications
    std::size_t fstart = task_index*(1+actiondata.outmat.ncols);
    std::size_t nelements = actiondata.outmat.bookeeping[fstart];
    for(unsigned i=0; i<nelements; ++i) {
      std::size_t base = arg1.start + actiondata.outmat.bookeeping[fstart+1+i];
      for(unsigned j=0; j<nmult; ++j) {
        vectors.arg2[j] = input.inputdata[ base + arg1.ncols*arg0.bookeeping[fpos+1+j] ];
      }
      MatrixElementOutput elem( 1, 2*nmult, output.values.data() + i, output.derivatives.data() + 2*nmult*i );
      T::calculate( input.noderiv, actiondata.funcinput, vectors, elem );
    }
  }
}

template <class T>
void MatrixTimesMatrix<T>::applyNonZeroRankForces( std::vector<double>& outforces ) {
  taskmanager.applyForces( outforces );
}

template <class T>
void MatrixTimesMatrix<T>::gatherForces( std::size_t task_index,
    const MatrixTimesMatrixInput<T>& actiondata,
    const ParallelActionsInput& input,
    const ForceInput& fdata,
    ForceOutput forces ) {
  ArgumentBookeepingHolder arg0( 0, input ), arg1( 1, input );
  std::size_t fpos = task_index*(1+arg0.ncols);
  std::size_t nmult = arg0.bookeeping[fpos];
  if( arg1.ncols<arg1.shape[1] ) {
    std::size_t fstart = task_index*(1+actiondata.outmat.ncols);
    std::size_t nelements = actiondata.outmat.bookeeping[fstart];
    for(unsigned i=0; i<nelements; ++i) {
      std::size_t nm = 0, ind = actiondata.outmat.bookeeping[fstart+1+i];
      double force = fdata.force[i];
      View<const double,helpers::dynamic_extent> deriv( fdata.deriv.data() + 2*nmult*i, 2*nmult );
      for(unsigned j=0; j<nmult; ++j) {
        std::size_t kind = arg0.bookeeping[fpos+1+j];
        std::size_t bstart = kind*(arg1.ncols + 1);
        std::size_t nr = arg1.bookeeping[bstart];
        for(unsigned k=0; k<nr; ++k) {
          if( arg1.bookeeping[bstart+1+k]==actiondata.outmat.bookeeping[fstart+1+i] ) {
            forces.thread_unsafe[task_index*arg0.ncols + j] += force*deriv[nm];
            forces.thread_safe[kind*arg1.ncols + k] += force*deriv[nmult+nm];
            nm++;
            break;
          }
        }
      }
    }
  } else {
    std::size_t fstart = task_index*(1+actiondata.outmat.ncols);
    std::size_t nelements = actiondata.outmat.bookeeping[fstart];
    for(unsigned i=0; i<nelements; ++i) {
      std::size_t ind = actiondata.outmat.bookeeping[fstart+1+i];
      double force = fdata.force[i];
      View<const double,helpers::dynamic_extent> deriv( fdata.deriv.data() + 2*nmult*i, 2*nmult );
      for(unsigned j=0; j<nmult; ++j) {
        forces.thread_unsafe[task_index*arg0.ncols + j] += force*deriv[j];
        forces.thread_safe[arg0.bookeeping[fpos+1+j]*arg1.ncols + ind] += force*deriv[nmult+j];
      }
    }
  }
}

}
}
#endif
