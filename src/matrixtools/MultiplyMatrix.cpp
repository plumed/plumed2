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

class MultiplyMatrix : public ActionWithInputMatrices {
private:
  Matrix<double> mat1, mat2, matout;
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit MultiplyMatrix(const ActionOptions&);
/// Do the calculation
  void completeMatrixOperations() override;
///
  void apply() override;
  double getForceOnMatrixElement( const unsigned& imat, const unsigned& jrow, const unsigned& krow ) const override { plumed_error(); }
};

PLUMED_REGISTER_ACTION(MultiplyMatrix,"MULTIPLY_MATRICES")

void MultiplyMatrix::registerKeywords( Keywords& keys ) {
  ActionWithInputMatrices::registerKeywords( keys );
}

MultiplyMatrix::MultiplyMatrix(const ActionOptions& ao):
  Action(ao),
  ActionWithInputMatrices(ao)
{
  if( getNumberOfArguments()!=2 ) error("should be two arguments for this action");
  if( getPntrToArgument(0)->getShape()[1]!=getPntrToArgument(1)->getShape()[0] ) error("number of columns in first matrix should equal number of rows in second");
 
  std::vector<unsigned> shape(2); shape[0]=getPntrToArgument(0)->getShape()[0]; 
  shape[1]=getPntrToArgument(1)->getShape()[1]; addValue( shape ); 
  mat1.resize( shape[0], getPntrToArgument(0)->getShape()[1] ); 
  mat2.resize( getPntrToArgument(1)->getShape()[0], shape[1] ); 
  matout.resize( shape[0], shape[1] );
}

void MultiplyMatrix::completeMatrixOperations() {
  // Retrieve the matrices from input and do the multiplication
  retrieveFullMatrix( 0, mat1 ); retrieveFullMatrix( 1, mat2 ); mult( mat1, mat2, matout );
  Value* myout=getPntrToOutput(0); unsigned nc = myout->getShape()[0], nr = myout->getShape()[1];
  for(unsigned i=0; i<nc; ++i) {
     for(unsigned j=0; j<nr; ++j) myout->set( i*nc+j, matout(i,j) );
  }
}

void MultiplyMatrix::apply() {
  // Do nothing if no forces were added
  if( !getPntrToOutput(0)->forcesWereAdded() ) return ;
  // Apply the forces on the matrix elements 
  Value* outval=getPntrToOutput(0); std::vector<unsigned> shape( outval->getShape() );
  Value* arg1=getPntrToArgument(0); Value* arg2=getPntrToArgument(1); unsigned ktot=arg1->getShape()[1];
  for(unsigned i=0; i<shape[0]; ++i) {
      for(unsigned j=0; j<shape[1]; ++j) {
          double force = outval->getForce( i*shape[1] + j );
          for(unsigned k=0; k<ktot; ++k ) {
              arg1->addForce( i*ktot + k, force*arg2->get(k*shape[1] + j) );
              arg2->addForce( k*shape[1]+j, force*arg1->get(i*ktot + k)  );
          }
      }
  }
}

}
}
