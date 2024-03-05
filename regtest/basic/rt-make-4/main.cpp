#include <iostream>
#include "plumed/tools/Matrix.h"
#include "plumed/tools/Tensor.h"

int main () {

  // Define symmetric matrix
  PLMD::Matrix<double> mat1(3,3); PLMD::OFile out; out.open("output");
  mat1(0,0)=1.0; mat1(0,1)=0.2; mat1(0,2)=0.3;
  mat1(1,0)=0.2; mat1(1,1)=0.2; mat1(1,2)=0.6;
  mat1(2,0)=0.3; mat1(2,1)=0.6; mat1(2,2)=0.4;

  // Test diagonalize
  std::vector<double> eigval(3); PLMD::Matrix<double> eigvec(3,3); 
  diagMat( mat1, eigval, eigvec ); 
  out<<"Eigenvalues "<<eigval[0]<<" "<<eigval[1]<<" "<<eigval[2]<<"\n";
  out<<"Eigenvectors : \n";
  for(unsigned i=0;i<3;++i){
      out<<eigvec(i,0)<<" "<<eigvec(i,1)<<" "<<eigvec(i,2)<<"\n";
  }

  // Test inverse
  out<<"Inverse : \n";
  PLMD::Matrix<double> inverse(3,3); Invert( mat1, inverse );
  for(unsigned i=0;i<3;++i){ 
      for(unsigned j=0;j<3;++j) out<<inverse(i,j)<<" ";
      out<<"\n";
  }

  // Test pseudoinverse 
  out<<"Pseudoinverse : \n";
  PLMD::Matrix<double> mat(3,2);
  mat(0,0)=0.1; mat(0,1)=0.2; 
  mat(1,0)=0.3; mat(1,1)=0.5;
  mat(2,0)=0.4; mat(2,1)=0.6;
  PLMD::Matrix<double> pseu(2,3);
  pseudoInvert( mat, pseu );
  for(unsigned i=0;i<pseu.nrows();++i){
     for(unsigned j=0;j<pseu.ncols();++j) out<<" "<<pseu(i,j);
     out<<"\n";
  }


/// The following is to test the fact that
/// dsyevr actually touches more elements of eval
/// than needed if there are identical eigenvectors.
/// This would trigger an error before fix
/// f1da0a9b3a13f904bd97d6f92a2fb5e1b6479ac0
  {
    PLMD::TensorGeneric<4,4> mat;
    PLMD::TensorGeneric<1,4> evec;
    PLMD::VectorGeneric<4> eval_underlying;

    auto eval = new(&eval_underlying[0]) PLMD::VectorGeneric<1>;

    mat[1][0]=mat[0][1]=3.0;
    mat[1][1]=5.0;
    mat[3][2]=mat[2][3]=3.0;
    mat[3][3]=5.0;
    PLMD::diagMatSym(mat,*eval,evec);

    out<<"Test diagmat (identical eigenvalues)\n";
/// Also eval_underlying[1] will be modified
    out<<eval_underlying[0]<<" "<<eval_underlying[1];
    out<<"\n";
  }

  {
    PLMD::TensorGeneric<4,4> mat;
    PLMD::TensorGeneric<1,4> evec;
    PLMD::VectorGeneric<4> eval_underlying;

    auto eval = new(&eval_underlying[0]) PLMD::VectorGeneric<1>;

    mat[1][0]=mat[0][1]=3.0;
    mat[1][1]=5.0;
    mat[3][2]=mat[2][3]=3.0;
    mat[3][3]=6.0;
    PLMD::diagMatSym(mat,*eval,evec);

    out<<"Test diagmat (different eigenvalues)\n";
    out<<eval_underlying[0]<<" "<<eval_underlying[1];
    out<<"\n";
  }

/// The following is to test the fact that
/// dsyevr actually modifies the diagonalized matrix
/// This would trigger an error before fix
/// f1da0a9b3a13f904bd97d6f92a2fb5e1b6479ac0
  {
    PLMD::TensorGeneric<4,4> mat;
    PLMD::TensorGeneric<4,4> evec;
    PLMD::VectorGeneric<4> eval;

    mat[1][0]=mat[0][1]=3.0;
    mat[1][1]=5.0;
    mat[0][3]=mat[3][0]=-1.0;
    mat[3][2]=mat[2][3]=3.0;
    mat[3][3]=6.0;

    out<<"Test diagmat (is matrix modified)\n";
    out<<"Before:\n";
    for(unsigned i=0;i<4;i++) {
      for(unsigned j=0;j<4;j++) out<<" "<<mat[i][j];
      out<<"\n";
    }

    PLMD::diagMatSym(mat,eval,evec);

    out<<"After:\n";
    for(unsigned i=0;i<4;i++) {
      for(unsigned j=0;j<4;j++) out<<" "<<mat[i][j];
      out<<"\n";
    }
  }
  out.close();
  return 0;
}
