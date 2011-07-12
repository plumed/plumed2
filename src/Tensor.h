#ifndef __PLUMED_Tensor_h
#define __PLUMED_Tensor_h

#include "Vector.h"

namespace PLMD{

/// 3x3 matrix of double.
/// Useful to simplify syntax. All the methods are inlined for better optimization.
class Tensor{
  double d[3][3];
public:
/// scale the tensor by a factor s
  friend Tensor operator*(double,const Tensor&);
/// return t1+t2
  friend Tensor operator+(const Tensor&,const Tensor&);
/// initialize the tensor to zero
  Tensor();
/// initialize a tensor as an external product of two Vector
  Tensor(const Vector&v1,const Vector&v2);
/// access element
  double & operator() (int i,int j);
/// access element
  const double & operator() (int i,int j)const;
/// set it to zero
  void clear();
/// returns the determinant
  double determinant()const;
/// return an identity tensor
  static Tensor identity();
/// return the matrix inverse
  Tensor inverse()const;
  Tensor transpose()const;
};

/// matrix multiplication
Tensor matmul(const Tensor&,const Tensor&);

Vector matmul(const Tensor&,const Vector&);
Vector matmul(const Vector&,const Tensor&);


inline
Tensor::Tensor(){
  for(int i=0;i<3;i++)for(int j=0;j<3;j++)d[i][j]=0.0;
}

inline
Tensor::Tensor(const Vector&v1,const Vector&v2){
  for(int i=0;i<3;i++)for(int j=0;j<3;j++)d[i][j]=v1[i]*v2[j];
}

inline
double & Tensor::operator() (int i,int j){
  return d[i][j];
}

inline
const double & Tensor::operator() (int i,int j)const{
  return d[i][j];
}

inline
void Tensor::clear(){
  d[0][0]=0.0;
  d[0][1]=0.0;
  d[0][2]=0.0;
  d[1][0]=0.0;
  d[1][1]=0.0;
  d[1][2]=0.0;
  d[2][0]=0.0;
  d[2][1]=0.0;
  d[2][2]=0.0;
}

inline
double Tensor::determinant()const{
  return
     d[0][0]*d[1][1]*d[2][2]
   + d[0][1]*d[1][2]*d[2][0]
   + d[0][2]*d[1][0]*d[2][1]
   - d[0][0]*d[1][2]*d[2][1]
   - d[0][1]*d[1][0]*d[2][2]
   - d[0][2]*d[1][1]*d[2][0];
}

inline
Tensor Tensor::identity(){
  Tensor t;
  t(0,0)=1.0;
  t(1,1)=1.0;
  t(2,2)=1.0;
  return t;
}

inline
Tensor Tensor::transpose()const{
  Tensor t;
  for(int i=0;i<3;i++)for(int j=0;j<3;j++) t(i,j)=d[j][i];
  return t;
}


inline
Tensor operator*(double s,const Tensor&v){
  Tensor t;
  for(int i=0;i<3;i++)for(int j=0;j<3;j++) t(i,j)=s*v(i,j);
  return t;
}

inline
Tensor Tensor::inverse()const{
  Tensor t;
  double invdet=1.0/determinant();
  for(int i=0;i<3;i++) for(int j=0;j<3;j++)
    t(j,i)=invdet*(   d[(i+1)%3][(j+1)%3]*d[(i+2)%3][(j+2)%3]
                     -d[(i+1)%3][(j+2)%3]*d[(i+2)%3][(j+1)%3]);
  return t;
}

inline
Tensor matmul(const Tensor&a,const Tensor&b){
  Tensor t;
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) {
    t(i,j)+=a(i,k)*b(k,j);
  }
  return t;
}

inline
Vector matmul(const Tensor&a,const Vector&b){
  Vector t;
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) t(i)+=a(i,j)*b(j);
  return t;
}

inline
Vector matmul(const Vector&a,const Tensor&b){
  Vector t;
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) t(i)+=a(i)*b(i,j);
  return t;
}

inline
Tensor operator+(const Tensor&t1,const Tensor&t2){
  Tensor t(t1);
  for(int i=0;i<3;i++) for(int j=0;j<3;j++) t(i,j)+=t2(i,j);
  return t;
}





}

#endif

