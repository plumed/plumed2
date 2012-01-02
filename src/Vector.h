#ifndef __PLUMED_Vector_h
#define __PLUMED_Vector_h

namespace PLMD{

/**
3d vector of double

Class implementing a standard three-component vector of double.
Vector elements are initialized to zero by default,
Useful to simplify syntax.
All the methods are inlined for better optimization.
Accepts both [] and () syntax for access.
Several functions are declared as friends even if not necessary so as to
properly appear in Doxygen documentation..

Example of usage
\verbatim
#include "Vector.h"

using namespace PLMD;

int main(){
  Vector v1;
  v1[0]=3.0;
// use equivalently () and [] syntax:
  v1(1)=5.0;
// initialize with components
  Vector v2=Vector(1.0,2.0,3.0);
  Vector v3=crossProduct(v1,v2);
  double d=dotProduct(v1,v2);
  v3+=v1;
  v2=v1+2.0*v3;
}
\endverbatim

*/
class Vector{
  double d[3];
public:
/// create it with components x0, x1 and x2
  Vector(double x0,double x1,double x2);
/// create it null
  Vector();
/// set it to zero
  void clear();
/// array-like access [i]
  double & operator[](unsigned i);
/// array-like access [i]
  const double & operator[](unsigned i)const;
/// parenthesis access (i)
  double & operator()(unsigned i);
/// parenthesis access (i)
  const double & operator()(unsigned i)const;
/// increment
  Vector& operator +=(const Vector& b);
/// decrement
  Vector& operator -=(const Vector& b);
/// multiply
  Vector& operator *=(double s);
/// divide
  Vector& operator /=(double s);
/// sign +
  Vector operator +()const;
/// sign -
  Vector operator -()const;
/// return v1+v2
  friend Vector operator+(const Vector&,const Vector&);
/// return v1-v2
  friend Vector operator-(const Vector&,const Vector&);
/// return s*v
  friend Vector operator*(double,const Vector&);
/// return v*s
  friend Vector operator*(const Vector&,double);
/// return v/s
  friend Vector operator/(const Vector&,double);
/// return v2-v1
  friend Vector delta(const Vector&v1,const Vector&v2);
/// return v1 .scalar. v2
  friend double dotProduct(const Vector&,const Vector&);
/// return v1 .vector. v2
  friend Vector crossProduct(const Vector&,const Vector&);
/// compute the squared modulo
  double modulo2()const;
/// compute the modulo
  double modulo()const;
};

inline
Vector::Vector(){
  d[0]=0.0;
  d[1]=0.0;
  d[2]=0.0;
}

inline
Vector::Vector(double x0,double x1,double x2){
  d[0]=x0;
  d[1]=x1;
  d[2]=x2;
}

inline
void Vector::clear(){
  d[0]=0.0;
  d[1]=0.0;
  d[2]=0.0;
}

inline
double & Vector::operator[](unsigned i){
  return d[i];
}

inline
const double & Vector::operator[](unsigned i)const{
  return d[i];
}

inline
double & Vector::operator()(unsigned i){
  return d[i];
}

inline
const double & Vector::operator()(unsigned i)const{
  return d[i];
}

inline
Vector& Vector::operator +=(const Vector& b){
  d[0]+=b(0);
  d[1]+=b(1);
  d[2]+=b(2);
  return *this;
}

inline
Vector& Vector::operator -=(const Vector& b){
  d[0]-=b(0);
  d[1]-=b(1);
  d[2]-=b(2);
  return *this;
}

inline
Vector& Vector::operator *=(double s){
  d[0]*=s;
  d[1]*=s;
  d[2]*=s;
  return *this;
}

inline
Vector& Vector::operator /=(double s){
  double x=1.0/s;
  d[0]*=x;
  d[1]*=x;
  d[2]*=x;
  return *this;
}

inline
Vector  Vector::operator +()const{
  return *this;
}

inline
Vector  Vector::operator -()const{
  return Vector(-d[0],-d[1],-d[2]);
}

inline
Vector operator+(const Vector&v1,const Vector&v2){
  return Vector(v1[0]+v2[0],v1[1]+v2[1],v1[2]+v2[2]);
}

inline
Vector operator-(const Vector&v1,const Vector&v2){
  return Vector(v1[0]-v2[0],v1[1]-v2[1],v1[2]-v2[2]);
}

inline
Vector operator*(double s,const Vector&v){
  return Vector(s*v[0],s*v[1],s*v[2]);
}

inline
Vector operator*(const Vector&v,double s){
  return s*v;
}

inline
Vector operator/(const Vector&v,double s){
  return v*(1.0/s);
}

inline
Vector delta(const Vector&v1,const Vector&v2){
  return Vector(v2[0]-v1[0],v2[1]-v1[1],v2[2]-v1[2]);
}

inline
double Vector::modulo2()const{
  return d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
}

inline
double dotProduct(const Vector& v1,const Vector& v2){
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

inline
Vector crossProduct(const Vector& v1,const Vector& v2){
  return Vector(
    v1[1]*v2[2]-v1[2]*v2[1],
    v1[2]*v2[0]-v1[0]*v2[2],
    v1[0]*v2[1]-v1[1]*v2[0]);
}

}


#endif

