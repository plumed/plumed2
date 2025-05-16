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

This action allows you to calculate the inverse of a real symmetric matrix.
If the matrix is not symmetric then the [dgetrf](https://www.netlib.org/lapack/explore-html-3.6.1/dd/d9a/group__double_g_ecomputational_ga0019443faea08275ca60a734d0593e60.html)
and [dgetri](https://www.netlib.org/lapack/explore-html-3.6.1/dd/d9a/group__double_g_ecomputational_ga56d9c860ce4ce42ded7f914fdb0683ff.html)
functions from the [LAPACK](https://www.netlib.org/lapack/explore-html/) library are used.
If the matrix is symmetric then we use [dsyevr](https://www.netlib.org/lapack/explore-html/d1/d56/group__heevr_gaa334ac0c11113576db0fc37b7565e8b5.html#gaa334ac0c11113576db0fc37b7565e8b5)
to find the eigenvalues. The inverse matrix of the input matrix $M$ with [eigendecomposition](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix).:

$$
M = Q \Lambda Q^{-1}
$$

is then found as:

$$
M^{-1} = Q \Lambda^{-1} Q^{-1}
$$

where $\Lambda^{-1}$ is the inverse of the (diagonal) matrix of eigenvalues $\Lambda$ and $Q$ is the matrix of eigenvectors.

The following example shows how this action is used in practise:

```plumed
c: DISTANCE_MATRIX ATOMS=1-4
ci: INVERT_MATRIX ARG=c
PRINT ARG=ci FILE=colvar
```

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

  std::vector<std::size_t> shape(2);
  shape[0]=shape[1]=getPntrToArgument(0)->getShape()[0];
  addValue( shape );
  setNotPeriodic();
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
