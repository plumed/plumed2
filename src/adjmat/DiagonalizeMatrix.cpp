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

class DiagonalizeMatrix :
  public ActionWithArguments,
  public ActionWithValue
{
private:
  std::vector<unsigned> desired_vectors;
  Matrix<double> mymatrix;
  std::vector<double> eigvals;
  Matrix<double> eigvecs;
  std::vector<double> forcesToApply;
  void diagonalizeMatrix();
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit DiagonalizeMatrix(const ActionOptions&);
/// Get the numebr of derivatives
  unsigned getNumberOfDerivatives() const { return getPntrToArgument(0)->getNumberOfValues(getLabel()); }
/// Do the calculation
  void calculate();
  void update();
  void runFinalJobs();
///
  void apply();
};

PLUMED_REGISTER_ACTION(DiagonalizeMatrix,"DIAGONALIZE")

void DiagonalizeMatrix::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys ); keys.use("ARG");
  keys.add("compulsory","VECTORS","all","the eigenvalues and vectors that you would like to calculate.  1=largest, 2=second largest and so on");
  keys.addOutputComponent("vals","default","the eigevalues of the input matrix");
  keys.addOutputComponent("vecs","default","the eigenvectors of the input matrix");
}

DiagonalizeMatrix::DiagonalizeMatrix(const ActionOptions& ao):
  Action(ao),
  ActionWithArguments(ao),
  ActionWithValue(ao)
{
  if( getNumberOfArguments()!=1 ) error("should only be one argument for this action");
  if( getPntrToArgument(0)->getRank()!=2 ) error("input argument for this action should be a matrix");
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
    addComponentWithDerivatives( "vals-" + num, eval_shape ); componentIsNotPeriodic( "vals-" + num );
    addComponent( "vecs-" + num, evec_shape ); componentIsNotPeriodic( "vecs-" + num );
    // Make sure eigenvalues are always stored
    getPntrToComponent( 2*i+1 )->alwaysStoreValues();
  }

  std::vector<unsigned> eigvecs_shape(2); eigvecs_shape[0]=eigvecs_shape[1]=getPntrToArgument(0)->getShape()[0];
  mymatrix.resize( eigvecs_shape[0], eigvecs_shape[1] ); eigvals.resize( eigvecs_shape[0] );
  eigvecs.resize( eigvecs_shape[0], eigvecs_shape[1] );
  // Now request the arguments to make sure we store things we need
  std::vector<Value*> args( getArguments() ); requestArguments(args, false );
  forcesToApply.resize( evec_shape[0]*evec_shape[0] );
}

void DiagonalizeMatrix::diagonalizeMatrix() {
  // Retrieve the matrix from input
  unsigned k = 0;
  for(unsigned i=0; i<mymatrix.nrows(); ++i) {
    for(unsigned j=0; j<mymatrix.ncols(); ++j) {
      mymatrix(i,j) = getPntrToArgument(0)->get( k ); k++;
    }
  }
  // Now diagonalize the matrix
  diagMat( mymatrix, eigvals, eigvecs );
  // And set the eigenvalues and eigenvectors
  for(unsigned i=0; i<desired_vectors.size(); ++i) {
    getPntrToOutput(2*i)->set( eigvals[ mymatrix.ncols()-desired_vectors[i]] );
    Value* evec_out = getPntrToOutput(2*i+1); unsigned vreq = mymatrix.ncols()-desired_vectors[i];
    for(unsigned j=0; j<mymatrix.ncols(); ++j) evec_out->set( j, eigvecs( vreq, j ) );
  }

  if( !doNotCalculateDerivatives() ) {
    for(unsigned i=0; i<mymatrix.nrows(); ++i) {
      for(unsigned j=0; j<mymatrix.ncols(); ++j) {
        unsigned nplace = i*mymatrix.nrows()+j;
        for(unsigned k=0; k<desired_vectors.size(); ++k) {
          unsigned ncol = mymatrix.ncols()-desired_vectors[k];
          getPntrToOutput(2*k)->addDerivative( nplace, eigvecs(ncol,i)*eigvecs(ncol,j) );
        }
      }
    }
  }
}

void DiagonalizeMatrix::calculate() {
  diagonalizeMatrix();
}

void DiagonalizeMatrix::update() {
  if( skipUpdate() ) return;
  diagonalizeMatrix();
}

void DiagonalizeMatrix::runFinalJobs() {
  if( skipUpdate() ) return;
  diagonalizeMatrix();
}

void DiagonalizeMatrix::apply() {
  if( doNotCalculateDerivatives() ) return;

  // Forces on eigenvalues
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned ss=0;
  if( getForcesFromValues( forcesToApply ) ) setForcesOnArguments( 0, forcesToApply, ss );

  // Check for forces on eigenvectors
  bool eigvec_forces=false;
  for(unsigned i=0; i<desired_vectors.size(); ++i) {
    if( getPntrToOutput(2*i+1)->forcesWereAdded() ) { eigvec_forces=true; break; }
  }
  if( !eigvec_forces ) return;

  // Forces on eigenvectors
  unsigned nn=0;
  for(unsigned j=0; j<mymatrix.nrows(); ++j) {
    for(unsigned k=0; k<mymatrix.ncols(); ++k) {
      double tmp1=0;
      for(unsigned i=0; i<desired_vectors.size(); ++i) {
        if( !getPntrToOutput(2*i+1)->forcesWereAdded() ) continue;

        unsigned ncol = mymatrix.ncols()-desired_vectors[i];
        for(unsigned n=0; n<mymatrix.nrows(); ++n) {
          double tmp2 = 0;
          for(unsigned m=0; m<mymatrix.nrows(); ++m) {
            if( m==ncol ) continue;
            tmp2 += eigvecs(m,n)*eigvecs(m,j)*eigvecs(ncol,k) / (eigvals[ncol]-eigvals[m]);
          }
          tmp1 += getPntrToOutput(2*i+1)->getForce(n) * tmp2;
        }
      }
      getPntrToArgument(0)->addForce( nn, tmp1 ); nn++;
    }
  }

}

}
}
