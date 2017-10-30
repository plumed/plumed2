#include <iostream>
#include "plumed/tools/Matrix.h"

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
  out.close();

  return 0;
}
