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
#include "core/ActionWithMatrix.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MCOLVAR MATRIX_PRODUCT
/*
Calculate the product of two matrices

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC ANALYSIS DISSIMILARITIES
/*
Calculate the matrix of dissimilarities between a trajectory of atomic configurations.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace matrixtools {

class MatrixTimesMatrix : public ActionWithMatrix {
private:
  bool squared;
  unsigned nderivatives;
  bool stored_matrix1, stored_matrix2;
public:
  static void registerKeywords( Keywords& keys );
  explicit MatrixTimesMatrix(const ActionOptions&);
  void prepare() override ;
  unsigned getNumberOfDerivatives();
  unsigned getNumberOfColumns() const override {
    return getConstPntrToComponent(0)->getShape()[1];
  }
  void getAdditionalTasksRequired( ActionWithVector* action, std::vector<unsigned>& atasks ) override ;
  void setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const ;
  void performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const override;
  void runEndOfRowJobs( const unsigned& ival, const std::vector<unsigned> & indices, MultiValue& myvals ) const override ;
};

PLUMED_REGISTER_ACTION(MatrixTimesMatrix,"MATRIX_PRODUCT")
PLUMED_REGISTER_ACTION(MatrixTimesMatrix,"DISSIMILARITIES")

void MatrixTimesMatrix::registerKeywords( Keywords& keys ) {
  ActionWithMatrix::registerKeywords(keys);
  keys.addInputKeyword("compulsory","ARG","matrix","the label of the two matrices from which the product is calculated");
  keys.addFlag("SQUARED",false,"calculate the squares of the dissimilarities (this option cannot be used with MATRIX_PRODUCT)");
  keys.setValueDescription("matrix","the product of the two input matrices");
}

MatrixTimesMatrix::MatrixTimesMatrix(const ActionOptions&ao):
  Action(ao),
  ActionWithMatrix(ao) {
  if( getNumberOfArguments()!=2 ) {
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
  std::vector<unsigned> shape(2);
  shape[0]=getPntrToArgument(0)->getShape()[0];
  shape[1]=getPntrToArgument(1)->getShape()[1];
  addValue( shape );
  setNotPeriodic();
  nderivatives = buildArgumentStore(0);
  std::string headstr=getFirstActionInChain()->getLabel();
  stored_matrix1 = getPntrToArgument(0)->ignoreStoredValue( headstr );
  stored_matrix2 = getPntrToArgument(1)->ignoreStoredValue( headstr );
  if( getName()=="DISSIMILARITIES" ) {
    parseFlag("SQUARED",squared);
    if( squared ) {
      log.printf("  calculating the squares of the dissimilarities \n");
    }
  } else {
    squared=true;
  }
}

unsigned MatrixTimesMatrix::getNumberOfDerivatives() {
  return nderivatives;
}

void MatrixTimesMatrix::prepare() {
  Value* myval = getPntrToComponent(0);
  if( myval->getShape()[0]==getPntrToArgument(0)->getShape()[0] && myval->getShape()[1]==getPntrToArgument(1)->getShape()[1] ) {
    return;
  }
  std::vector<unsigned> shape(2);
  shape[0]=getPntrToArgument(0)->getShape()[0];
  shape[1]=getPntrToArgument(1)->getShape()[1];
  myval->setShape(shape);
  if( myval->valueIsStored() ) {
    myval->reshapeMatrixStore( shape[1] );
  }
}

void MatrixTimesMatrix::getAdditionalTasksRequired( ActionWithVector* action, std::vector<unsigned>& atasks ) {

  ActionWithMatrix* adj=dynamic_cast<ActionWithMatrix*>( getPntrToArgument(0)->getPntrToAction() );
  if( !adj->isAdjacencyMatrix() ) {
    return;
  }
  adj->retrieveAtoms();
  adj->getAdditionalTasksRequired( action, atasks );
}

void MatrixTimesMatrix::setupForTask( const unsigned& task_index, std::vector<unsigned>& indices, MultiValue& myvals ) const {
  unsigned start_n = getPntrToArgument(0)->getShape()[0], size_v = getPntrToArgument(1)->getShape()[1];
  if( indices.size()!=size_v+1 ) {
    indices.resize( size_v+1 );
  }
  for(unsigned i=0; i<size_v; ++i) {
    indices[i+1] = start_n + i;
  }
  myvals.setSplitIndex( size_v + 1 );
}

void MatrixTimesMatrix::performTask( const std::string& controller, const unsigned& index1, const unsigned& index2, MultiValue& myvals ) const {
  unsigned ostrn = getConstPntrToComponent(0)->getPositionInStream(), ind2=index2;
  if( index2>=getPntrToArgument(0)->getShape()[0] ) {
    ind2 = index2 - getPntrToArgument(0)->getShape()[0];
  }

  Value* myarg = getPntrToArgument(0);
  unsigned nmult=myarg->getRowLength(index1);
  double matval=0;
  std::vector<double>  dvec1(nmult), dvec2(nmult);
  for(unsigned i=0; i<nmult; ++i) {
    unsigned kind = myarg->getRowIndex( index1, i );
    double val1 = getElementOfMatrixArgument( 0, index1, kind, myvals );
    double val2 = getElementOfMatrixArgument( 1, kind, ind2, myvals );
    if( getName()=="DISSIMILARITIES" ) {
      double tmp = getPntrToArgument(0)->difference(val2, val1);
      matval += tmp*tmp;
      if( !squared ) {
        dvec1[i] = 2*tmp;
        dvec2[i] = -2*tmp;
        continue;
      } else {
        val2 = -2*tmp;
        val1 = 2*tmp;
      }
    } else {
      matval+= val1*val2;
    }

    if( doNotCalculateDerivatives() ) {
      continue;
    }

    addDerivativeOnMatrixArgument( stored_matrix1, 0, 0, index1, kind, val2, myvals );
    addDerivativeOnMatrixArgument( stored_matrix2, 0, 1, kind, ind2, val1, myvals );
  }
  // And add this part of the product
  if( !squared ) {
    matval = sqrt(matval);
  }
  myvals.addValue( ostrn, matval );
  if( squared || doNotCalculateDerivatives() ) {
    return;
  }

  for(unsigned i=0; i<nmult; ++i) {
    unsigned kind = myarg->getRowIndex( index1, i );
    addDerivativeOnMatrixArgument( stored_matrix1, 0, 0, index1, kind, dvec1[i]/(2*matval), myvals );
    addDerivativeOnMatrixArgument( stored_matrix2, 0, 1, kind, ind2, dvec2[i]/(2*matval), myvals );
  }
}

void MatrixTimesMatrix::runEndOfRowJobs( const unsigned& ival, const std::vector<unsigned> & indices, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() || !matrixChainContinues() ) {
    return ;
  }

  unsigned mat1s = ival*getPntrToArgument(0)->getShape()[1];
  unsigned nmult = getPntrToArgument(0)->getShape()[1], ss = getPntrToArgument(1)->getShape()[1];
  unsigned nmat = getConstPntrToComponent(0)->getPositionInMatrixStash(), nmat_ind = myvals.getNumberOfMatrixRowDerivatives( nmat );
  std::vector<unsigned>& matrix_indices( myvals.getMatrixRowDerivativeIndices( nmat ) );
  unsigned ntwo_atoms = myvals.getSplitIndex();
  for(unsigned j=0; j<nmult; ++j) {
    matrix_indices[nmat_ind] = mat1s + j;
    nmat_ind++;
    for(unsigned i=1; i<ntwo_atoms; ++i) {
      unsigned ind2 = indices[i];
      if( ind2>=getPntrToArgument(0)->getShape()[0] ) {
        ind2 = indices[i] - getPntrToArgument(0)->getShape()[0];
      }
      matrix_indices[nmat_ind] = arg_deriv_starts[1] + j*ss + ind2;
      nmat_ind++;
    }
  }
  myvals.setNumberOfMatrixRowDerivatives( nmat, nmat_ind );
}

}
}
