#ifndef __PLUMED_Tensor_h
#define __PLUMED_Tensor_h

#include "Vector.h"

namespace PLMD{

/// 3x3 matrix of double.
/// Useful to simplify syntax. All the methods are inlined for better optimization.
/// Several functions are declared as friends even if not necessary so as to
/// properly appear in Doxygen documentation..
class Tensor{
  double d[3][3];
public:
/// Small utility class which just contains a pointer to the tensor so as
/// to be able to use the [][] syntax (const version)
/// See C++ FAQ 13.12
  class Const_row{
    friend class Tensor; // this so as to allow only Tensor to instantiate Const_row
                         // the user should not manipulate it directly
    const Tensor& t;
    const unsigned i;
    Const_row(const Tensor&t,unsigned i); // constructor is private and cannot be manipulated by the user
  public:
  /// access element
    const double & operator[] (unsigned j)const;
  };
/// Small utility class which just contains a pointer to the tensor so as
/// to be able to use the [][] syntax
/// See C++ FAQ 13.12
  class Row{
    friend class Tensor; // this so as to allow only Tensor to instantiate Const_row
                         // the user should not manipulate it directly
    Tensor& t;
    const unsigned i;
    Row(Tensor&t,unsigned i); // constructor is private and cannot be manipulated by the user
  public:
  /// access element
    double & operator[] (unsigned j);
  };
/// initialize the tensor to zero
  Tensor();
/// initialize a tensor as an external product of two Vector
  Tensor(const Vector&v1,const Vector&v2);
/// initialize a tensor with 9 values, in standard C order
  Tensor(double,double,double,double,double,double,double,double,double);
/// set it to zero
  void clear();
/// access element
  double & operator() (unsigned i,unsigned j);
/// access element
  const double & operator() (unsigned i,unsigned j)const;
/// access element (with [][] syntax)
  Row operator[] (unsigned i);
/// access element (with [][] syntax)
  Const_row operator[] (unsigned i)const;
/// increment
  Tensor& operator +=(const Tensor& b);
/// decrement
  Tensor& operator -=(const Tensor& b);
/// multiply
  Tensor& operator *=(double);
/// divide
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
/// return the transpose matrix
  Tensor transpose()const;
/// matrix-matrix multiplication
  friend Tensor matmul(const Tensor&,const Tensor&);
/// matrix-vector multiplication
  friend Vector matmul(const Tensor&,const Vector&);
/// vector-matrix multiplication
  friend Vector matmul(const Vector&,const Tensor&);
/// returns the determinant of a tensor
  friend double determinant(const Tensor&);
/// returns the inverse of a tensor (same as inverse())
  friend Tensor inverse(const Tensor&);
/// returns the transpose of a tensor (same as transpose())
  friend Tensor transpose(const Tensor&);
/// returns the transpose of a tensor (same as Tensor(const Vector&,const Vector&))
  friend Tensor extProduct(const Vector&,const Vector&);
};

inline
Tensor::Const_row::Const_row(const Tensor&t,unsigned i):
  t(t),i(i){}

inline
Tensor::Row::Row(Tensor&t,unsigned i):
  t(t),i(i){}

inline
Tensor::Tensor(){
  for(unsigned i=0;i<3;i++)for(unsigned j=0;j<3;j++)d[i][j]=0.0;
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
double & Tensor::operator() (unsigned i,unsigned j){
  return d[i][j];
}

inline
const double & Tensor::operator() (unsigned i,unsigned j)const{
  return d[i][j];
}

inline
const double & Tensor::Const_row::operator[] (unsigned j)const{
  return t(i,j);
}

inline
double & Tensor::Row::operator[] (unsigned j){
  return t(i,j);
}

inline
Tensor::Row Tensor::operator[] (unsigned i){
  return Row(*this,i);
}

inline
Tensor::Const_row Tensor::operator[] (unsigned i)const{
  return Const_row(*this,i);
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
  return Tensor(
      -d[0][0],
      -d[0][1],
      -d[0][2],
      -d[1][0],
      -d[1][1],
      -d[1][2],
      -d[2][0],
      -d[2][1],
      -d[2][2]);
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
  for(unsigned i=0;i<3;i++)for(unsigned j=0;j<3;j++) t(i,j)=d[j][i];
  return t;
}

inline
Tensor Tensor::inverse()const{
  Tensor t;
  double invdet=1.0/determinant();
  for(unsigned i=0;i<3;i++) for(unsigned j=0;j<3;j++)
    t(j,i)=invdet*(   d[(i+1)%3][(j+1)%3]*d[(i+2)%3][(j+2)%3]
                     -d[(i+1)%3][(j+2)%3]*d[(i+2)%3][(j+1)%3]);
  return t;
}

inline
Tensor matmul(const Tensor&a,const Tensor&b){
  Tensor t;
  for(unsigned i=0;i<3;i++) for(unsigned j=0;j<3;j++) for(unsigned k=0;k<3;k++) {
    t(i,j)+=a(i,k)*b(k,j);
  }
  return t;
}

inline
Vector matmul(const Tensor&a,const Vector&b){
  Vector t;
  for(unsigned i=0;i<3;i++) for(unsigned j=0;j<3;j++) t(i)+=a(i,j)*b(j);
  return t;
}

inline
Vector matmul(const Vector&a,const Tensor&b){
  Vector t;
  for(unsigned i=0;i<3;i++) for(unsigned j=0;j<3;j++) t(i)+=a(j)*b(j,i);
  return t;
}

inline
double determinant(const Tensor&t){
  return t.determinant();
}

inline
Tensor inverse(const Tensor&t){
  return t.inverse();
}

inline
Tensor transpose(const Tensor&t){
  return t.transpose();
}

inline
Tensor extProduct(const Vector&v1,const Vector&v2){
  return Tensor(v1,v2);
}

}

#endif

