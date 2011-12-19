#ifndef __PLUMED_Tensor_h
#define __PLUMED_Tensor_h

#include "Vector.h"

namespace PLMD{

/// 3x3 matrix of double.
/// Useful to simplify syntax. All the methods are inlined for better optimization.
class Tensor{
  double d[3][3];
public:
/// initialize the tensor to zero
  Tensor();
/// initialize a tensor as an external product of two Vector
  Tensor(const Vector&v1,const Vector&v2);
/// initialize a tensor with 9 values, in standard C order
  Tensor(double,double,double,double,double,double,double,double,double);
/// set it to zero
  void clear();
/// access element
  double & operator() (int i,int j);
/// access element
  const double & operator() (int i,int j)const;
  Tensor& operator +=(const Tensor& b);
  Tensor& operator -=(const Tensor& b);
  Tensor& operator *=(double);
  Tensor& operator /=(double);
/// return +t
  Tensor operator +()const;
/// return -t
  Tensor operator -()const;
/// return t1+t2
  friend Tensor operator+(const Tensor&,const Tensor&);
/// return t1+t2
  friend Tensor operator-(const Tensor&,const Tensor&);
/// scale the tensor by a factor s
  friend Tensor operator*(double,const Tensor&);
/// scale the tensor by a factor s
  friend Tensor operator*(const Tensor&,double s);
/// scale the tensor by a factor 1/s
  friend Tensor operator/(const Tensor&,double s);
/// returns the determinant
  double determinant()const;
/// return an identity tensor
  static Tensor identity();
/// return the matrix inverse
  Tensor inverse()const;
  Tensor transpose()const;
/// matrix multiplication
  friend Tensor matmul(const Tensor&,const Tensor&);
  friend Vector matmul(const Tensor&,const Vector&);
  friend Vector matmul(const Vector&,const Tensor&);
};

inline
Tensor::Tensor(){
  for(int i=0;i<3;i++)for(int j=0;j<3;j++)d[i][j]=0.0;
}

inline
Tensor::Tensor(const Vector&v1,const Vector&v2){
  d[0][0]=v1(0)*v2(0);
  d[0][1]=v1(0)*v2(1);
  d[0][2]=v1(0)*v2(2);
  d[1][0]=v1(1)*v2(0);
  d[1][1]=v1(1)*v2(1);
  d[1][2]=v1(1)*v2(2);
  d[2][0]=v1(2)*v2(0);
  d[2][1]=v1(2)*v2(1);
  d[2][2]=v1(2)*v2(2);
}

inline
Tensor::Tensor(double d00,double d01,double d02,double d10,double d11,double d12,double d20,double d21,double d22){
  d[0][0]=d00;
  d[0][1]=d01;
  d[0][2]=d02;
  d[1][0]=d10;
  d[1][1]=d11;
  d[1][2]=d12;
  d[2][0]=d20;
  d[2][1]=d21;
  d[2][2]=d22;
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
double & Tensor::operator() (int i,int j){
  return d[i][j];
}

inline
const double & Tensor::operator() (int i,int j)const{
  return d[i][j];
}

inline
Tensor& Tensor::operator +=(const Tensor& b){
  d[0][0]+=b(0,0);
  d[0][1]+=b(0,1);
  d[0][2]+=b(0,2);
  d[1][0]+=b(1,0);
  d[1][1]+=b(1,1);
  d[1][2]+=b(1,2);
  d[2][0]+=b(2,0);
  d[2][1]+=b(2,1);
  d[2][2]+=b(2,2);
  return *this;
}

inline
Tensor& Tensor::operator -=(const Tensor& b){
  d[0][0]-=b(0,0);
  d[0][1]-=b(0,1);
  d[0][2]-=b(0,2);
  d[1][0]-=b(1,0);
  d[1][1]-=b(1,1);
  d[1][2]-=b(1,2);
  d[2][0]-=b(2,0);
  d[2][1]-=b(2,1);
  d[2][2]-=b(2,2);
  return *this;
}

inline
Tensor& Tensor::operator *=(double s){
  d[0][0]*=s;
  d[0][1]*=s;
  d[0][2]*=s;
  d[1][0]*=s;
  d[1][1]*=s;
  d[1][2]*=s;
  d[2][0]*=s;
  d[2][1]*=s;
  d[2][2]*=s;
  return *this;
}

inline
Tensor& Tensor::operator /=(double s){
  return (*this)*=1.0/s;
}

inline
Tensor Tensor::operator+()const{
  return *this;
}

inline
Tensor Tensor::operator-()const{
  return -1.0*(*this);
}

inline
Tensor operator+(const Tensor&t1,const Tensor&t2){
  Tensor t(t1);
  t+=t2;
  return t;
}

inline
Tensor operator-(const Tensor&t1,const Tensor&t2){
  Tensor t(t1);
  t-=t2;
  return t;
}

inline
Tensor operator*(const Tensor&t1,double s){
  Tensor t(t1);
  t*=s;
  return t;
}

inline
Tensor operator*(double s,const Tensor&t1){
  return t1*s;
}

inline
Tensor operator/(const Tensor&t1,double s){
  return t1*(1.0/s);
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





}

#endif

