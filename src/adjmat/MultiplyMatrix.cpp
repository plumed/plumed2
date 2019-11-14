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

class MultiplyMatrix :
  public ActionWithArguments,
  public ActionWithValue
{
private:
  Matrix<double> mymatrix;
  Matrix<double> inverse;
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit MultiplyMatrix(const ActionOptions&);
/// Get the numebr of derivatives
  unsigned getNumberOfDerivatives() const { return getPntrToArgument(0)->getNumberOfValues(getLabel()) + getPntrToArgument(1)->getNumberOfValues(getLabel()); }
/// Do the calculation
  void calculate();
///
  void apply();
};

PLUMED_REGISTER_ACTION(MultiplyMatrix,"MULTIPLY_MATRICES")

void MultiplyMatrix::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys ); keys.use("ARG");
}

MultiplyMatrix::MultiplyMatrix(const ActionOptions& ao):
  Action(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao)
{
  if( getNumberOfArguments()!=2 ) error("should be two arguments for this action");
  if( getPntrToArgument(0)->getRank()!=2 ) error("first input argument for this action should be a matrix");
  if( getPntrToArgument(1)->getRank()!=2 ) error("second input argument for this action should be a matrix");
  if( getPntrToArgument(0)->getShape()[1]!=getPntrToArgument(1)->getShape()[0] ) error("number of columns in first matrix should equal number of rows in second");
 
  std::vector<unsigned> shape(2); shape[0]=getPntrToArgument(0)->getShape()[0]; 
  shape[1]=getPntrToArgument(1)->getShape()[1]; addValue( shape ); setNotPeriodic();
  // Now request the arguments to make sure we store things we need
  std::vector<Value*> args( getArguments() ); requestArguments(args, false );
}

void MultiplyMatrix::calculate() {
  // Retrieve the matrix from input
  unsigned k = 0; Value* myout=getPntrToOutput(0); 
  unsigned nc = myout->getShape()[0], nr = myout->getShape()[1];
  Value* mat1=getPntrToArgument(0); Value* mat2=getPntrToArgument(1); unsigned nn=mat1->getShape()[1];
  for(unsigned i=0; i<nc; ++i) {
    for(unsigned j=0; j<nr; ++j) {
      double val=0; 
      for(unsigned k=0; k<nn; ++k) val += mat1->get( i*nc+k )*mat2->get( k*nn + j );
      myout->set( i*nc+j, val );
    }
  }

  if( !doNotCalculateDerivatives() ) error("derivatives of inverse matrix have not been implemented");
}

void MultiplyMatrix::apply() {
  if( doNotCalculateDerivatives() ) return;
  error("derivatives of inverse matrix have not been implemented");
}

}
}
