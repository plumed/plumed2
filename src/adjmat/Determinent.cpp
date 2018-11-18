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

class Determinant :
  public ActionWithArguments,
  public ActionWithValue
{
private:
  Matrix<double> mymatrix;
  std::vector<double> eigvals;
  Matrix<double> eigvecs;
  std::vector<double> forcesToApply;
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit Determinant(const ActionOptions&);
/// Get the numebr of derivatives
  unsigned getNumberOfDerivatives() const { return getPntrToArgument(0)->getNumberOfValues(getLabel()); }
/// Do the calculation
  void calculate();
///
  void apply();
};

PLUMED_REGISTER_ACTION(Determinant,"DETERMINANT")

void Determinant::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys ); keys.use("ARG");
}

Determinant::Determinant(const ActionOptions& ao):
  Action(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao)
{
  if( getNumberOfArguments()!=1 ) error("should only be one argument for this action");
  if( getPntrToArgument(0)->getRank()!=2 ) error("input argument for this action should be a matrix");
  if( getPntrToArgument(0)->getShape()[0]!=getPntrToArgument(0)->getShape()[1] ) error("input matrix should be square");

  std::vector<unsigned> shape; addValue( shape ); setNotPeriodic();
  std::vector<unsigned> eigvecs_shape(2); eigvecs_shape[0]=eigvecs_shape[1]=getPntrToArgument(0)->getShape()[0];
  mymatrix.resize( eigvecs_shape[0], eigvecs_shape[1] ); eigvals.resize( eigvecs_shape[0] );
  eigvecs.resize( eigvecs_shape[0], eigvecs_shape[1] );
  // Now request the arguments to make sure we store things we need
  std::vector<Value*> args( getArguments() ); requestArguments(args, false );
  forcesToApply.resize( eigvecs_shape[0]*eigvecs_shape[0] );
}

void Determinant::calculate() {
  // Retrieve the matrix from input
  unsigned k = 0;
  for(unsigned i=0; i<mymatrix.nrows(); ++i) {
    for(unsigned j=0; j<mymatrix.ncols(); ++j) {
      mymatrix(i,j) = getPntrToArgument(0)->get( k ); k++;
    }
  }
  // Now diagonalize the matrix
  diagMat( mymatrix, eigvals, eigvecs );
  // The determinent is the product of the eigenvalues
  double det = 1; for(unsigned i=0;i<eigvals.size();++i) det *= eigvals[i];
  getPntrToOutput(0)->set(det);

  if( !doNotCalculateDerivatives() ) {
    Value* valout = getPntrToOutput(0);
    for(unsigned i=0; i<mymatrix.nrows(); ++i) {
      for(unsigned j=0; j<mymatrix.ncols(); ++j) {
        unsigned nplace = i*mymatrix.nrows()+j;
        for(unsigned k=0; k<eigvals.size(); ++k) valout->addDerivative( nplace, det*eigvecs(k,i)*eigvecs(k,j)/eigvals[k] );
      }
    }
  }
}

void Determinant::apply() {
  if( doNotCalculateDerivatives() ) return;
  // Forces on determinent 
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned ss=0;
  if( getForcesFromValues( forcesToApply ) ) setForcesOnArguments( 0, forcesToApply, ss );
}

}
}
