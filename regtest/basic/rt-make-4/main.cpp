#include <iostream>
#include "plumed/tools/Matrix.h"
#include "plumed/tools/Tensor.h"

template <typename precision>
void matrixTest (PLMD::OFile& out, const std::string& fmt) {

  // Define symmetric matrix
  PLMD::Matrix<precision> mat1(3,3);
  mat1(0,0)=1.0;
  mat1(0,1)=0.2;
  mat1(0,2)=0.3;
  mat1(1,0)=0.2;
  mat1(1,1)=0.2;
  mat1(1,2)=0.6;
  mat1(2,0)=0.3;
  mat1(2,1)=0.6;
  mat1(2,2)=0.4;

  // Test diagonalize
  std::vector<precision> eigval(3);
  PLMD::Matrix<precision> eigvec(3,3);
  diagMat( mat1, eigval, eigvec );
  out.printf(("Eigenvalues "+fmt+" "+fmt+" "+fmt+"\n").c_str(),
             eigval[0],eigval[1],eigval[2]);
  out<<"Eigenvectors : \n";
  for(unsigned i=0; i<3; ++i) {
    out.printf((fmt+" "+fmt+" "+fmt+"\n").c_str(),
               eigvec(i,0), eigvec(i,1), eigvec(i,2));
  }

  // Test inverse
  out<<"Inverse : \n";
  PLMD::Matrix<precision> inverse(3,3);
  Invert( mat1, inverse );
  for(unsigned i=0; i<3; ++i) {
    for(unsigned j=0; j<3; ++j) {
      out.printf((" "+fmt).c_str(),inverse(i,j));
    }
    out<<"\n";
  }

  // Test pseudoinverse
  out<<"Pseudoinverse : \n";
  PLMD::Matrix<precision> mat(3,2);
  mat(0,0)=0.1;
  mat(0,1)=0.2;
  mat(1,0)=0.3;
  mat(1,1)=0.5;
  mat(2,0)=0.4;
  mat(2,1)=0.6;
  PLMD::Matrix<precision> pseu(2,3);
  pseudoInvert( mat, pseu );
  for(unsigned i=0; i<pseu.nrows(); ++i) {
    for(unsigned j=0; j<pseu.ncols(); ++j) {
      out.printf((" "+fmt).c_str(),pseu(i,j));
    }
    out<<"\n";
  }


/// The following is to test the fact that
/// dsyevr actually touches more elements of eval
/// than needed if there are identical eigenvectors.
/// This would trigger an error before fix
/// f1da0a9b3a13f904bd97d6f92a2fb5e1b6479ac0
  {
    PLMD::TensorTyped<precision,4,4> mat;
    PLMD::TensorTyped<precision,1,4> evec;
    PLMD::VectorTyped<precision,4> eval_underlying;

    //The both of the two following lines are scary, but on my machine are equivalent
    //PLMD::Vector1d* eval=reinterpret_cast<PLMD::Vector1d*>(eval_underlying.data());
    auto eval = new(&eval_underlying[0]) PLMD::VectorTyped<precision,1>;

    mat[1][0]=mat[0][1]=3.0;
    mat[1][1]=5.0;
    mat[3][2]=mat[2][3]=3.0;
    mat[3][3]=5.0;
    PLMD::diagMatSym(mat,*eval,evec);

    out<<"Test diagmat (identical eigenvalues)\n";
/// Also eval_underlying[1] will be modified
    out.printf((fmt+" "+fmt+"\n").c_str(),
               eval_underlying[0], eval_underlying[1]);
  }

  {
    PLMD::TensorTyped<precision,4,4> mat;
    PLMD::TensorTyped<precision,1,4> evec;
    PLMD::VectorTyped<precision,4> eval_underlying;

    auto eval = new(&eval_underlying[0]) PLMD::VectorTyped<precision,1>;

    mat[1][0]=mat[0][1]=3.0;
    mat[1][1]=5.0;
    mat[3][2]=mat[2][3]=3.0;
    mat[3][3]=6.0;
    PLMD::diagMatSym(mat,*eval,evec);

    out<<"Test diagmat (different eigenvalues)\n";
    out.printf((fmt+" "+fmt+"\n").c_str(),
               eval_underlying[0], eval_underlying[1]);
  }

/// The following is to test the fact that
/// dsyevr actually modifies the diagonalized matrix
/// This would trigger an error before fix
/// f1da0a9b3a13f904bd97d6f92a2fb5e1b6479ac0
  {
    PLMD::TensorTyped<precision,4,4> mat;
    PLMD::TensorTyped<precision,4,4> evec;
    PLMD::VectorTyped<precision,4> eval;

    mat[1][0]=mat[0][1]=3.0;
    mat[1][1]=5.0;
    mat[0][3]=mat[3][0]=-1.0;
    mat[3][2]=mat[2][3]=3.0;
    mat[3][3]=6.0;

    out<<"Test diagmat (is matrix modified)\n";
    out<<"Before:\n";
    for(unsigned i=0; i<4; i++) {
      for(unsigned j=0; j<4; j++) {
        out.printf((" "+fmt).c_str(),mat[i][j]);
      }
      out<<"\n";
    }

    PLMD::diagMatSym(mat,eval,evec);

    out<<"After:\n";
    for(unsigned i=0; i<4; i++) {
      for(unsigned j=0; j<4; j++) {
        out.printf((" "+fmt).c_str(),mat[i][j]);
      }
      out<<"\n";
    }
  }
}
int main () {
  PLMD::OFile out;
  out.open("output");
  matrixTest<double>(out,"%8.6f");
  out.close();
  PLMD::OFile out_float;
  out_float.open("output_float");
  matrixTest<float>(out_float,"%8.4f");
  out_float.close();
  return 0;
}
