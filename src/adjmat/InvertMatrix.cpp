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
#include "core/ActionWithArguments.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"
#include "tools/Matrix.h"

namespace PLMD {
namespace adjmat {

class InvertMatrix :
  public ActionWithArguments,
  public ActionWithValue
{
private:
  Matrix<double> mymatrix;
  Matrix<double> inverse;
  std::vector<double> forcesToApply;
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit InvertMatrix(const ActionOptions&);
/// Get the numebr of derivatives
  unsigned getNumberOfDerivatives() const { return getPntrToArgument(0)->getNumberOfValues(getLabel()); }
/// Do the calculation
  void calculate();
///
  void apply();
};

PLUMED_REGISTER_ACTION(InvertMatrix,"INVERT_MATRIX")

void InvertMatrix::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys ); keys.use("ARG");
}

InvertMatrix::InvertMatrix(const ActionOptions& ao):
  Action(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao)
{
  if( getNumberOfArguments()!=1 ) error("should only be one argument for this action");
  if( getPntrToArgument(0)->getRank()!=2 ) error("input argument for this action should be a matrix");
  if( getPntrToArgument(0)->getShape()[0]!=getPntrToArgument(0)->getShape()[1] ) error("input matrix should be square");
 
  std::vector<unsigned> shape(2); shape[0]=shape[1]=getPntrToArgument(0)->getShape()[0];
  addValue( shape ); setNotPeriodic();

  mymatrix.resize( shape[0], shape[1] ); inverse.resize( shape[0], shape[1] );
  // Now request the arguments to make sure we store things we need
  std::vector<Value*> args( getArguments() ); requestArguments(args, false );
  forcesToApply.resize( shape[0]*shape[1] );
}

void InvertMatrix::calculate() {
  // Retrieve the matrix from input
  unsigned k = 0;
  for(unsigned i=0; i<mymatrix.nrows(); ++i) {
    for(unsigned j=0; j<mymatrix.ncols(); ++j) {
      mymatrix(i,j) = getPntrToArgument(0)->get( k ); k++;
    }
  }
  // Now invert the matrix
  Invert( mymatrix, inverse );
  // And set the inverse
  k = 0;
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
