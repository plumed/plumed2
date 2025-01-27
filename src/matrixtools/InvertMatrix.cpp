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
#include "MatrixOperationBase.h"
#include "core/ActionSetup.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MCOLVAR INVERT_MATRIX
/*
Calculate the inverse of the input matrix

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace matrixtools {

class InvertMatrix : public MatrixOperationBase {
private:
  bool input_is_constant;
  Matrix<double> mymatrix;
  Matrix<double> inverse;
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit InvertMatrix(const ActionOptions&);
///
  unsigned getNumberOfDerivatives() override {
    return 0;
  }
/// Do the calculation
  void calculate() override;
///
  void apply() override;
  double getForceOnMatrixElement( const unsigned& jrow, const unsigned& krow ) const override {
    plumed_error();
  }
};

PLUMED_REGISTER_ACTION(InvertMatrix,"INVERT_MATRIX")

void InvertMatrix::registerKeywords( Keywords& keys ) {
  MatrixOperationBase::registerKeywords( keys );
  keys.setValueDescription("matrix","the inverse of the input matrix");
}

InvertMatrix::InvertMatrix(const ActionOptions& ao):
  Action(ao),
  MatrixOperationBase(ao),
  input_is_constant(false) {
  if( getPntrToArgument(0)->getShape()[0]!=getPntrToArgument(0)->getShape()[1] ) {
    error("input matrix should be square");
  }

  ActionSetup* as = dynamic_cast<ActionSetup*>( getPntrToArgument(0)->getPntrToAction() );
  if(as) {
    input_is_constant=true;
  }

  std::vector<unsigned> shape(2);
  shape[0]=shape[1]=getPntrToArgument(0)->getShape()[0];
  addValue( shape );
  setNotPeriodic();
  getPntrToComponent(0)->buildDataStore();
  getPntrToComponent(0)->reshapeMatrixStore( shape[1] );
  mymatrix.resize( shape[0], shape[1] );
  inverse.resize( shape[0], shape[1] );
}

void InvertMatrix::calculate() {
  // Retrieve the matrix from input
  retrieveFullMatrix( mymatrix );
  // Now invert the matrix
  Invert( mymatrix, inverse );
  // And set the inverse
  unsigned k = 0;
  Value* myval=getPntrToComponent(0);
  for(unsigned i=0; i<mymatrix.nrows(); ++i) {
    for(unsigned j=0; j<mymatrix.ncols(); ++j) {
      myval->set( k, inverse(i,j) );
      k++;
    }
  }

  if( !doNotCalculateDerivatives() && !input_is_constant ) {
    error("derivatives of inverse matrix have not been implemented");
  }
}

void InvertMatrix::apply() {
  if( doNotCalculateDerivatives() || input_is_constant ) {
    return;
  }
  error("derivatives of inverse matrix have not been implemented");
}

}
}
