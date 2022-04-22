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

class DiagonalizeMatrix : public ActionWithInputMatrices {
private:
  std::vector<unsigned> desired_vectors;
  Matrix<double> mymatrix;
  std::vector<double> eigvals;
  Matrix<double> eigvecs;
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit DiagonalizeMatrix(const ActionOptions&);
/// This is required to set the number of derivatives for the eigenvalues
  unsigned getNumberOfDerivatives() const override { return getPntrToArgument(0)->getNumberOfValues(); }
/// 
  unsigned getNumberOfFinalTasks() override { return getPntrToArgument(0)->getShape()[0]; }
/// Do the calculation
  void completeMatrixOperations() override;
///
  void apply();
///
  double getForceOnMatrixElement( const unsigned& imat, const unsigned& jrow, const unsigned& krow ) const override;
};

PLUMED_REGISTER_ACTION(DiagonalizeMatrix,"DIAGONALIZE")

void DiagonalizeMatrix::registerKeywords( Keywords& keys ) {
  ActionWithInputMatrices::registerKeywords( keys );
  keys.add("compulsory","VECTORS","all","the eigenvalues and vectors that you would like to calculate.  1=largest, 2=second largest and so on");
  keys.addOutputComponent("vals","default","the eigevalues of the input matrix");
  keys.addOutputComponent("vecs","default","the eigenvectors of the input matrix");
}

DiagonalizeMatrix::DiagonalizeMatrix(const ActionOptions& ao):
  Action(ao),
  ActionWithInputMatrices(ao)
{
  if( getNumberOfArguments()!=1 ) error("should only be one argument for this action");
  if( getPntrToArgument(0)->getShape()[0]!=getPntrToArgument(0)->getShape()[1] ) error("input matrix should be square");

  std::vector<std::string> eigv; parseVector("VECTORS",eigv);
  if( eigv.size()>1 ) {
    Tools::interpretRanges(eigv); desired_vectors.resize( eigv.size() );
    for(unsigned i=0; i<eigv.size(); ++i) Tools::convert( eigv[i], desired_vectors[i] );
  } else  {
    if( eigv.size()==0 ) error("missing input to VECTORS keyword");
    unsigned ivec;
    if( Tools::convert( eigv[0], ivec ) ) {
      desired_vectors.resize(1); desired_vectors[0]=ivec;
    } else if( eigv[0]=="all") {
      desired_vectors.resize( getPntrToArgument(0)->getShape()[0] );
      for(unsigned i=0; i<desired_vectors.size(); ++i) desired_vectors[i] = i + 1;
    } else error("input to VECTOR keyword should be list of numbers or all");
  }

  std::string num; std::vector<unsigned> eval_shape(0);
  std::vector<unsigned> evec_shape(1); evec_shape[0] = getPntrToArgument(0)->getShape()[0];
  for(unsigned i=0; i<desired_vectors.size(); ++i) {
    Tools::convert( desired_vectors[i], num );
    addComponent( "vals-" + num, eval_shape ); componentIsNotPeriodic( "vals-" + num );
    addComponent( "vecs-" + num, evec_shape ); componentIsNotPeriodic( "vecs-" + num );
    // Make sure eigenvalues are always stored
    getPntrToComponent( 2*i+1 )->alwaysStoreValues();
  }

  std::vector<unsigned> eigvecs_shape(2); eigvecs_shape[0]=eigvecs_shape[1]=getPntrToArgument(0)->getShape()[0];
  mymatrix.resize( eigvecs_shape[0], eigvecs_shape[1] ); eigvals.resize( eigvecs_shape[0] ); eigvecs.resize( eigvecs_shape[0], eigvecs_shape[1] ); 
}

void DiagonalizeMatrix::completeMatrixOperations() {
  if( getPntrToArgument(0)->getShape()[0]==0 ) return ;
  // Resize stuff that might need resizing
  unsigned nvals=getPntrToArgument(0)->getShape()[0]; 
  if( eigvals.size()!=nvals ) { 
      mymatrix.resize( nvals, nvals ); eigvals.resize( nvals ); 
      eigvecs.resize( nvals, nvals ); std::vector<unsigned> shape(1); shape[0]=nvals;
      for(unsigned i=0; i<desired_vectors.size(); ++i) {
         if( getPntrToComponent( 2*i+1 )->getShape()[0]!=nvals ) getPntrToComponent( 2*i+1 )->setShape( shape );
      }
  }

  // Retrieve the matrix from input
  retrieveFullMatrix( 0, mymatrix );
  // Now diagonalize the matrix
  diagMat( mymatrix, eigvals, eigvecs );
  // And set the eigenvalues and eigenvectors
  for(unsigned i=0; i<desired_vectors.size(); ++i) {
    getPntrToOutput(2*i)->set( eigvals[ mymatrix.ncols()-desired_vectors[i]] );
    Value* evec_out = getPntrToOutput(2*i+1); unsigned vreq = mymatrix.ncols()-desired_vectors[i];
    for(unsigned j=0; j<mymatrix.ncols(); ++j) evec_out->set( j, eigvecs( vreq, j ) );
  }
}

void DiagonalizeMatrix::apply() {
  if( doNotCalculateDerivatives() ) return;

  // Check for forces on eigenvectors or eigenvalues
  bool forces=false;
  for(unsigned i=0; i<desired_vectors.size(); ++i) {
    if( getPntrToOutput(2*i)->forcesWereAdded() || getPntrToOutput(2*i+1)->forcesWereAdded() ) { forces=true; break; }
  }
  if( !forces ) return;
  // Get forces from matrix
  applyForceOnMatrix( 0 );
}

double DiagonalizeMatrix::getForceOnMatrixElement( const unsigned& imat, const unsigned& jrow, const unsigned& kcol ) const {
  double ff = 0;
  for(unsigned i=0; i<desired_vectors.size(); ++i) {
    // Deal with forces on eigenvalues
    if( getPntrToOutput(2*i)->forcesWereAdded() ) {
        unsigned ncol = mymatrix.ncols()-desired_vectors[i]; 
        ff += getPntrToOutput(2*i)->getForce(0)*eigvecs(ncol,jrow)*eigvecs(ncol,kcol);
    }
    // And forces on eigenvectors
    if( !getPntrToOutput(2*i+1)->forcesWereAdded() ) continue;

    unsigned ncol = mymatrix.ncols()-desired_vectors[i];
    for(unsigned n=0; n<mymatrix.nrows(); ++n) {
      double tmp2 = 0;
      for(unsigned m=0; m<mymatrix.nrows(); ++m) {
        if( m==ncol ) continue;
        tmp2 += eigvecs(m,n)*eigvecs(m,jrow)*eigvecs(ncol,kcol) / (eigvals[ncol]-eigvals[m]);
      }
      ff += getPntrToOutput(2*i+1)->getForce(n) * tmp2;
    }
  }
  return ff;
}

}
}
