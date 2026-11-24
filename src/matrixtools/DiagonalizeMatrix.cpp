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
#include "core/ActionRegister.h"

#include <algorithm>

//+PLUMEDOC ANALYSIS DIAGONALIZE
/*
Calculate the eigenvalues and eigenvectors of a square matrix

This action allows you to use the [dsyevr](https://www.netlib.org/lapack/explore-html/d1/d56/group__heevr_gaa334ac0c11113576db0fc37b7565e8b5.html#gaa334ac0c11113576db0fc37b7565e8b5)
function from the [LAPACK](https://www.netlib.org/lapack/explore-html/) library to calculate the [eigenvalues and eigenvectors](https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix)
 of a real symmetric matrix. For example, the following input can be used to calculate and print all four eigenvalues of the input [DISTANCE_MATRIX](DISTANCE_MATRIX.md).

```plumed
d: DISTANCE_MATRIX GROUP=1-4
diag: DIAGONALIZE ARG=d
PRINT ARG=diag.vals-1,diag.vals-2,diag.vals-3,diag.vals-4 FILE=colvar
```

If you wish to calculate only a subset of the eigenvalues and eigenvectors you would use an input like the one shown below.  This input calculates and outputs the largest eigenvalue
and its corresponding eigenvector.

```plumed
d: DISTANCE_MATRIX GROUP=1-4
diag: DIAGONALIZE ARG=d VECTORS=1
PRINT ARG=diag.vals-1,diag.vecs-1 FILE=colvar
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace matrixtools {

class DiagonalizeMatrix : public MatrixOperationBase {
private:
  std::vector<unsigned> desired_vectors;
  Matrix<double> mymatrix;
  std::vector<double> eigvals;
  Matrix<double> eigvecs;
/// For each requested eigenpair that has forces: the index of the eigenpair in eigvecs,
/// the force on the eigenvalue and a vector used to compute the forces from the eigenvector
  struct ForceTerm {
    unsigned ncol;
    double evalForce;
    bool hasEvecForce;
    std::vector<double> w;
  };
  std::vector<ForceTerm> forceTerms;
public:
  static void registerKeywords( Keywords& keys );
/// Constructor
  explicit DiagonalizeMatrix(const ActionOptions&);
/// This is required to set the number of derivatives for the eigenvalues
  unsigned getNumberOfDerivatives() override {
    return getPntrToArgument(0)->getNumberOfValues();
  }
///
  void prepare() override ;
///
  void calculate() override ;
/// Prepare the quantities needed to compute the forces on the matrix elements
  void apply() override ;
///
  double getForceOnMatrixElement( const unsigned& jrow, const unsigned& krow ) const override;
};

PLUMED_REGISTER_ACTION(DiagonalizeMatrix,"DIAGONALIZE")

void DiagonalizeMatrix::registerKeywords( Keywords& keys ) {
  MatrixOperationBase::registerKeywords( keys );
  keys.add("compulsory","VECTORS","all","the eigenvalues and vectors that you would like to calculate.  1=largest, 2=second largest and so on");
  keys.addOutputComponent("vals","default","scalar","the eigevalues of the input matrix");
  keys.addOutputComponent("vecs","default","vector","the eigenvectors of the input matrix");
}

DiagonalizeMatrix::DiagonalizeMatrix(const ActionOptions& ao):
  Action(ao),
  MatrixOperationBase(ao) {
  if( getPntrToArgument(0)->getShape()[0]!=getPntrToArgument(0)->getShape()[1] ) {
    error("input matrix should be square");
  }

  std::vector<std::string> eigv;
  parseVector("VECTORS",eigv);
  if( eigv.size()>1 ) {
    Tools::interpretRanges(eigv);
    desired_vectors.resize( eigv.size() );
    for(unsigned i=0; i<eigv.size(); ++i) {
      Tools::convert( eigv[i], desired_vectors[i] );
    }
  } else  {
    if( eigv.size()==0 ) {
      error("missing input to VECTORS keyword");
    }
    unsigned ivec;
    if( eigv[0]=="all" ) {
      desired_vectors.resize( getPntrToArgument(0)->getShape()[0] );
      for(unsigned i=0; i<desired_vectors.size(); ++i) {
        desired_vectors[i] = i + 1;
      }
    } else {
      Tools::convert( eigv[0], ivec );
      desired_vectors.resize(1);
      desired_vectors[0]=ivec;
    }
  }

  std::string num;
  std::vector<std::size_t> eval_shape(0);
  std::vector<std::size_t> evec_shape(1);
  evec_shape[0] = getPntrToArgument(0)->getShape()[0];
  for(unsigned i=0; i<desired_vectors.size(); ++i) {
    Tools::convert( desired_vectors[i], num );
    addComponent( "vals-" + num, eval_shape );
    componentIsNotPeriodic( "vals-" + num );
    addComponent( "vecs-" + num, evec_shape );
    componentIsNotPeriodic( "vecs-" + num );
    // Make sure eigenvalues are always stored
  }

  std::vector<unsigned> eigvecs_shape(2);
  eigvecs_shape[0]=eigvecs_shape[1]=getPntrToArgument(0)->getShape()[0];
  mymatrix.resize( eigvecs_shape[0], eigvecs_shape[1] );
  eigvals.resize( eigvecs_shape[0] );
  eigvecs.resize( eigvecs_shape[0], eigvecs_shape[1] );
}

void DiagonalizeMatrix::prepare() {
  std::vector<std::size_t> shape(1);
  shape[0]=getPntrToArgument(0)->getShape()[0];
  for(unsigned i=0; i<desired_vectors.size(); ++i) {
    if( getPntrToComponent( 2*i+1 )->getShape()[0]!=shape[0] ) {
      getPntrToComponent( 2*i+1 )->setShape( shape );
    }
  }

}

void DiagonalizeMatrix::calculate() {
  if( getPntrToArgument(0)->getShape()[0]==0 ) {
    return ;
  }
  // Resize stuff that might need resizing
  unsigned nvals=getPntrToArgument(0)->getShape()[0];
  if( eigvals.size()!=nvals ) {
    mymatrix.resize( nvals, nvals );
    eigvals.resize( nvals );
    eigvecs.resize( nvals, nvals );
  }

  // Retrieve the matrix from input
  retrieveFullMatrix( mymatrix );
  // Now diagonalize the matrix. The forces on an eigenvector depend on all the eigenvectors,
  // so they are all computed when derivatives are needed. Otherwise only the eigenpairs from
  // the smallest requested one to the largest eigenvalue are computed, which is faster
  // when only a few of the largest eigenvalues are requested (e.g. VECTORS=1)
  if( doNotCalculateDerivatives() ) {
    const unsigned maxvec=*std::max_element( desired_vectors.begin(), desired_vectors.end() );
    diagMatRange( mymatrix, nvals-maxvec, nvals-1, eigvals, eigvecs );
  } else {
    diagMat( mymatrix, eigvals, eigvecs );
  }
  // And set the eigenvalues and eigenvectors
  for(unsigned i=0; i<desired_vectors.size(); ++i) {
    getPntrToComponent(2*i)->set( eigvals[ mymatrix.ncols()-desired_vectors[i]] );
    Value* evec_out = getPntrToComponent(2*i+1);
    unsigned vreq = mymatrix.ncols()-desired_vectors[i];
    for(unsigned j=0; j<mymatrix.ncols(); ++j) {
      evec_out->set( j, eigvecs( vreq, j ) );
    }
  }
}

void DiagonalizeMatrix::apply() {
  if( doNotCalculateDerivatives() ) {
    return;
  }
  // The derivative of eigenvector c with respect to element (j,k) of the matrix is
  //   d v_c[n] / d A(j,k) = sum_{m!=c} v_m[n] v_m[j] v_c[k] / (lambda_c-lambda_m)
  // so the force on (j,k) from the forces f on v_c is v_c[k]*w[j], with
  //   w[j] = sum_{m!=c} v_m[j] (v_m.f) / (lambda_c-lambda_m)
  // computing w once here costs O(N^2), instead of O(N^2) for each matrix element
  const unsigned nvals=eigvals.size();
  forceTerms.clear();
  for(unsigned i=0; i<desired_vectors.size(); ++i) {
    const bool evalForce=getConstPntrToComponent(2*i)->forcesWereAdded();
    const bool evecForce=getConstPntrToComponent(2*i+1)->forcesWereAdded();
    if( !evalForce && !evecForce ) {
      continue;
    }
    ForceTerm ft;
    ft.ncol=nvals-desired_vectors[i];
    ft.evalForce=( evalForce ? getConstPntrToComponent(2*i)->getForce(0) : 0.0 );
    ft.hasEvecForce=evecForce;
    if( evecForce ) {
      const Value* evec=getConstPntrToComponent(2*i+1);
      ft.w.assign(nvals,0.0);
      for(unsigned m=0; m<nvals; ++m) {
        if( m==ft.ncol ) {
          continue;
        }
        double vf=0;
        for(unsigned n=0; n<nvals; ++n) {
          vf += eigvecs(m,n)*evec->getForce(n);
        }
        vf /= (eigvals[ft.ncol]-eigvals[m]);
        for(unsigned j=0; j<nvals; ++j) {
          ft.w[j] += vf*eigvecs(m,j);
        }
      }
    }
    forceTerms.push_back(std::move(ft));
  }
  MatrixOperationBase::apply();
}

double DiagonalizeMatrix::getForceOnMatrixElement( const unsigned& jrow, const unsigned& kcol ) const {
  double ff = 0;
  for(const auto & ft : forceTerms) {
    // Forces on eigenvalues
    ff += ft.evalForce*eigvecs(ft.ncol,jrow)*eigvecs(ft.ncol,kcol);
    // Forces on eigenvectors
    if( ft.hasEvecForce ) {
      ff += ft.w[jrow]*eigvecs(ft.ncol,kcol);
    }
  }
  return ff;
}

}
}
