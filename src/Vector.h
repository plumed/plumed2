#ifndef __PLUMED_Vector_h
#define __PLUMED_Vector_h

#include <cmath>
#include <cassert>

namespace PLMD{

/**
Class implementing fixed size vectors of doubles

This class implements a vector of doubles with size fixed at
compile time. It is useful for small fixed size objects (e.g.
3d vectors) as it does not waste space to store the vector size.
Moreover, as the compiler knows the size, it can be completely
opimized inline.
All the methods are inlined for better optimization.
Vector elements are initialized to zero by default. Notice that
this means that constructor is a bit slow. This point might change
in future if we find performance issues.
Accepts both [] and () syntax for access.
Several functions are declared as friends even if not necessary so as to
properly appear in Doxygen documentation.

Aliases are defined to simplify common declarations (Vector, Vector2d, Vector3d, Vector4d).
Also notice that some operations are only available for 3 dimensional vectors.

Example of usage
\verbatim
#include "Vector.h"

using namespace PLMD;

int main(){
  VectorGeneric<3> v1;
  v1[0]=3.0;
// use equivalently () and [] syntax:
  v1(1)=5.0;
// initialize with components
  VectorGeneric<3> v2=VectorGeneric<3>(1.0,2.0,3.0);
  VectorGeneric<3> v3=crossProduct(v1,v2);
  double d=dotProduct(v1,v2);
  v3+=v1;
  v2=v1+2.0*v3;
}
\endverbatim

*/

template <unsigned n>
class VectorGeneric{
  double d[n];
public:
/// Create it with preassigned components.
/// Only available for sizes 2, 3 and 4
  VectorGeneric(double,double);
  VectorGeneric(double,double,double);
  VectorGeneric(double,double,double,double);
/// create it null
  VectorGeneric();
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
  VectorGeneric& operator +=(const VectorGeneric& b);
/// decrement
  VectorGeneric& operator -=(const VectorGeneric& b);
/// multiply
  VectorGeneric& operator *=(double s);
/// divide
  VectorGeneric& operator /=(double s);
/// sign +
  VectorGeneric operator +()const;
/// sign -
  VectorGeneric operator -()const;
/// return v1+v2
  template<unsigned m>
  friend VectorGeneric<m> operator+(const VectorGeneric<m>&,const VectorGeneric<m>&);
/// return v1-v2
  template<unsigned m>
  friend VectorGeneric<m> operator-(const VectorGeneric<m>&,const VectorGeneric<m>&);
/// return s*v
  template<unsigned m>
  friend VectorGeneric<m> operator*(double,const VectorGeneric<m>&);
/// return v*s
  template<unsigned m>
  friend VectorGeneric<m> operator*(const VectorGeneric<m>&,double);
/// return v/s
  template<unsigned m>
  friend VectorGeneric<m> operator/(const VectorGeneric<m>&,double);
/// return v2-v1
  template<unsigned m>
  friend VectorGeneric<m> delta(const VectorGeneric<m>&v1,const VectorGeneric<m>&v2);
/// return v1 .scalar. v2
  template<unsigned m>
  friend double dotProduct(const VectorGeneric<m>&,const VectorGeneric<m>&);
/// return v1 .vector. v2
/// Only available for size 3
  friend VectorGeneric<3> crossProduct(const VectorGeneric<3>&,const VectorGeneric<3>&);
/// compute the squared modulo
  double modulo2()const;
/// Compute the modulo.
/// Shortcut for sqrt(v.modulo2())
  double modulo()const;
/// friend version of modulo2 (to simplify some syntax)
  template<unsigned m>
  friend double modulo2(const VectorGeneric<m>&);
/// friend version of modulo (to simplify some syntax)
  template<unsigned m>
  friend double modulo(const VectorGeneric<m>&);
};

template<>
inline
VectorGeneric<2>:: VectorGeneric(double x0,double x1){
  d[0]=x0;
  d[1]=x1;
}

template<>
inline
VectorGeneric<3>:: VectorGeneric(double x0,double x1,double x2){
  d[0]=x0;
  d[1]=x1;
  d[2]=x2;
}

template<>
inline
VectorGeneric<4>:: VectorGeneric(double x0,double x1,double x2,double x3){
  d[0]=x0;
  d[1]=x1;
  d[2]=x2;
  d[3]=x3;
}

template <unsigned n>
VectorGeneric<n>::VectorGeneric(){
  for(unsigned i=0;i<n;i++) d[i]=0.0;
}

template <unsigned n>
void VectorGeneric<n>::clear(){
  for(unsigned i=0;i<n;i++) d[i]=0.0;
}

template <unsigned n>
double & VectorGeneric<n>::operator[](unsigned i){
  return d[i];
}

template <unsigned n>
const double & VectorGeneric<n>::operator[](unsigned i)const{
  return d[i];
}

template <unsigned n>
double & VectorGeneric<n>::operator()(unsigned i){
  return d[i];
}

template <unsigned n>
const double & VectorGeneric<n>::operator()(unsigned i)const{
  return d[i];
}

template <unsigned n>
VectorGeneric<n>& VectorGeneric<n>::operator +=(const VectorGeneric<n>& b){
  for(unsigned i=0;i<n;i++) d[i]+=b.d[i];
  return *this;
}

template <unsigned n>
VectorGeneric<n>& VectorGeneric<n>::operator -=(const VectorGeneric<n>& b){
  for(unsigned i=0;i<n;i++) d[i]-=b.d[i];
  return *this;
}

template <unsigned n>
VectorGeneric<n>& VectorGeneric<n>::operator *=(double s){
  for(unsigned i=0;i<n;i++) d[i]*=s;
  return *this;
}

template <unsigned n>
VectorGeneric<n>& VectorGeneric<n>::operator /=(double s){
  double x=1.0/s;
  return (*this)*=x;
}

template <unsigned n>
VectorGeneric<n>  VectorGeneric<n>::operator +()const{
  return *this;
}

template <unsigned n>
VectorGeneric<n> VectorGeneric<n>::operator -()const{
  VectorGeneric<n> r;
  for(unsigned i=0;i<n;i++) r[i]=-d[i];
  return r;
}

template <unsigned n>
VectorGeneric<n> operator+(const VectorGeneric<n>&v1,const VectorGeneric<n>&v2){
  VectorGeneric<n> v(v1);
  return v+=v2;
}

template <unsigned n>
VectorGeneric<n> operator-(const VectorGeneric<n>&v1,const VectorGeneric<n>&v2){
  VectorGeneric<n> v(v1);
  return v-=v2;
}

template <unsigned n>
VectorGeneric<n> operator*(double s,const VectorGeneric<n>&v){
  VectorGeneric<n> vv(v);
  return vv*=s;
}

template <unsigned n>
VectorGeneric<n> operator*(const VectorGeneric<n>&v,double s){
  return s*v;
}

template <unsigned n>
VectorGeneric<n> operator/(const VectorGeneric<n>&v,double s){
  return v*(1.0/s);
}

template <unsigned n>
VectorGeneric<n> delta(const VectorGeneric<n>&v1,const VectorGeneric<n>&v2){
  return v2-v1;
}

template <unsigned n>
double VectorGeneric<n>::modulo2()const{
  double r=0.0;
  for(unsigned i=0;i<n;i++) r+=d[i]*d[i];
  return r;
}

template <unsigned n>
double dotProduct(const VectorGeneric<n>& v1,const VectorGeneric<n>& v2){
  double r=0.0;
  for(unsigned i=0;i<n;i++) r+=v1[i]*v2[i];
  return r;
}

// faster, specialized version for 2d
inline
double dotProduct(const VectorGeneric<2>& v1,const VectorGeneric<2>& v2){
  return v1[0]*v2[0]+v1[1]*v2[1];
}

template<>
inline
double VectorGeneric<2>::modulo2()const{
  return d[0]*d[0]+d[1]*d[1];
}

// faster, specialized version for 3d
inline
double dotProduct(const VectorGeneric<3>& v1,const VectorGeneric<3>& v2){
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
}

template<>
inline
double VectorGeneric<3>::modulo2()const{
  return d[0]*d[0]+d[1]*d[1]+d[2]*d[2];
}

// faster, specialized version for 4d
inline
double dotProduct(const VectorGeneric<4>& v1,const VectorGeneric<4>& v2){
  return v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]+v1[3]*v2[3];
}

template<>
inline
double VectorGeneric<4>::modulo2()const{
  return d[0]*d[0]+d[1]*d[1]+d[2]*d[2]+d[3]*d[3];
}

inline
VectorGeneric<3> crossProduct(const VectorGeneric<3>& v1,const VectorGeneric<3>& v2){
 return VectorGeneric<3>(
   v1[1]*v2[2]-v1[2]*v2[1],
   v1[2]*v2[0]-v1[0]*v2[2],
   v1[0]*v2[1]-v1[1]*v2[0]);
}

template<unsigned n>
double VectorGeneric<n>::modulo()const{
  return sqrt(modulo2());
}

template<unsigned n>
double modulo2(const VectorGeneric<n>&v){
  return v.modulo2();
}

template<unsigned n>
double modulo(const VectorGeneric<n>&v){
  return v.modulo();
}


/// Alias for two dimensional vectors;
typedef VectorGeneric<2> Vector2d;
/// Alias for three dimensional vectors;
typedef VectorGeneric<3> Vector3d;
/// Alias for four dimensional vectors;
typedef VectorGeneric<4> Vector4d;
/// Alias for three dimensional vectors;
typedef Vector3d Vector;

}

#endif

