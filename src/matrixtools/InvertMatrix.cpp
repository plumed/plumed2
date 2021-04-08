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

class InvertMatrix : public ActionWithInputMatrices {
private:
  Matrix<double> mymatrix;
  Matrix<double> inverse;
  std::vector<double> forcesToApply;
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit InvertMatrix(const ActionOptions&);
/// Do the calculation
  void completeMatrixOperations() override;
///
  void apply();
};

PLUMED_REGISTER_ACTION(InvertMatrix,"INVERT_MATRIX")

void InvertMatrix::registerKeywords( Keywords& keys ) {
  ActionWithInputMatrices::registerKeywords( keys );
}

InvertMatrix::InvertMatrix(const ActionOptions& ao):
  Action(ao),
  ActionWithInputMatrices(ao)
{
  if( getNumberOfArguments()!=1 ) error("should only be one argument for this action");
  if( getPntrToArgument(0)->getShape()[0]!=getPntrToArgument(0)->getShape()[1] ) error("input matrix should be square");
 
  std::vector<unsigned> shape(2); shape[0]=shape[1]=getPntrToArgument(0)->getShape()[0];
  addValue( shape ); setNotPeriodic(); getPntrToComponent(0)->alwaysStoreValues();

  mymatrix.resize( shape[0], shape[1] ); inverse.resize( shape[0], shape[1] );
  // Now request the arguments to make sure we store things we need
  forcesToApply.resize( shape[0]*shape[1] );
}

void InvertMatrix::completeMatrixOperations() {
  // Retrieve the matrix from input
  retrieveFullMatrix( 0, mymatrix );
  // Now invert the matrix
  Invert( mymatrix, inverse );
  // And set the inverse
  unsigned k = 0;
  for(unsigned i=0; i<mymatrix.nrows(); ++i) {
    for(unsigned j=0; j<mymatrix.ncols(); ++j) {
        getPntrToOutput(0)->set( k, inverse(i,j) ); k++;
    }
  } 

  if( !doNotCalculateDerivatives() ) error("derivatives of inverse matrix have not been implemented");
}

void InvertMatrix::apply() {
  if( doNotCalculateDerivatives() ) return;
  error("derivatives of inverse matrix have not been implemented");
}

}
}
