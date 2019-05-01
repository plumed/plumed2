/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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
#ifndef __PLUMED_tools_Vector_h
#define __PLUMED_tools_Vector_h

#include <cmath>
#include <iosfwd>
#include <array>
#include "LoopUnroller.h"

namespace PLMD {

/**
\ingroup TOOLBOX
Class implementing fixed size vectors of doubles

\tparam n The number of elements of the vector.

This class implements a vector of doubles with size fixed at
compile time. It is useful for small fixed size objects (e.g.
3d vectors) as it does not waste space to store the vector size.
Moreover, as the compiler knows the size, it can be completely
opimized inline.
All the methods are inlined for better optimization and
all the loops are explicitly unrolled using PLMD::LoopUnroller class.
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
class VectorGeneric {
  std::array<double,n> d;
/// Auxiliary private function for constructor
  void auxiliaryConstructor();
/// Auxiliary private function for constructor
  template<typename... Args>
  void auxiliaryConstructor(double first,Args... arg);
public:
/// Constructor accepting n double parameters.
/// Can be used as Vector<3>(1.0,2.0,3.0) or Vector<2>(2.0,3.0).
/// In case a wrong number of parameters is given, a static assertion will fail.
  template<typename... Args>
  VectorGeneric(double first,Args... arg);
/// create it null
  VectorGeneric();
/// set it to zero
  void zero();
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
/// << operator.
/// Allows printing vector `v` with `std::cout<<v;`
  template<unsigned m>
  friend std::ostream & operator<<(std::ostream &os, const VectorGeneric<m>&);
};

template <unsigned n>
void VectorGeneric<n>::auxiliaryConstructor()
{}

template <unsigned n>
template<typename... Args>
void VectorGeneric<n>::auxiliaryConstructor(double first,Args... arg)
{
  d[n-(sizeof...(Args))-1]=first;
  auxiliaryConstructor(arg...);
}

template <unsigned n>
template<typename... Args>
VectorGeneric<n>::VectorGeneric(double first,Args... arg)
{
  static_assert((sizeof...(Args))+1==n,"you are trying to initialize a Vector with the wrong number of arguments");
  auxiliaryConstructor(first,arg...);
}

template <unsigned n>
void VectorGeneric<n>::zero() {
  LoopUnroller<n>::_zero(d.data());
}

template <unsigned n>
VectorGeneric<n>::VectorGeneric() {
  LoopUnroller<n>::_zero(d.data());
}

template <unsigned n>
double & VectorGeneric<n>::operator[](unsigned i) {
  return d[i];
}

template <unsigned n>
const double & VectorGeneric<n>::operator[](unsigned i)const {
  return d[i];
}

template <unsigned n>
double & VectorGeneric<n>::operator()(unsigned i) {
  return d[i];
}

template <unsigned n>
const double & VectorGeneric<n>::operator()(unsigned i)const {
  return d[i];
}

template <unsigned n>
VectorGeneric<n>& VectorGeneric<n>::operator +=(const VectorGeneric<n>& b) {
  LoopUnroller<n>::_add(d.data(),b.d.data());
  return *this;
}

template <unsigned n>
VectorGeneric<n>& VectorGeneric<n>::operator -=(const VectorGeneric<n>& b) {
  LoopUnroller<n>::_sub(d.data(),b.d.data());
  return *this;
}

template <unsigned n>
VectorGeneric<n>& VectorGeneric<n>::operator *=(double s) {
  LoopUnroller<n>::_mul(d.data(),s);
  return *this;
}

template <unsigned n>
VectorGeneric<n>& VectorGeneric<n>::operator /=(double s) {
  LoopUnroller<n>::_mul(d.data(),1.0/s);
  return *this;
}

template <unsigned n>
VectorGeneric<n>  VectorGeneric<n>::operator +()const {
  return *this;
}

template <unsigned n>
VectorGeneric<n> VectorGeneric<n>::operator -()const {
  VectorGeneric<n> r;
  LoopUnroller<n>::_neg(r.d.data(),d.data());
  return r;
}

template <unsigned n>
VectorGeneric<n> operator+(const VectorGeneric<n>&v1,const VectorGeneric<n>&v2) {
  VectorGeneric<n> v(v1);
  return v+=v2;
}

template <unsigned n>
VectorGeneric<n> operator-(const VectorGeneric<n>&v1,const VectorGeneric<n>&v2) {
  VectorGeneric<n> v(v1);
  return v-=v2;
}

template <unsigned n>
VectorGeneric<n> operator*(double s,const VectorGeneric<n>&v) {
  VectorGeneric<n> vv(v);
  return vv*=s;
}

template <unsigned n>
VectorGeneric<n> operator*(const VectorGeneric<n>&v,double s) {
  return s*v;
}

template <unsigned n>
VectorGeneric<n> operator/(const VectorGeneric<n>&v,double s) {
  return v*(1.0/s);
}

template <unsigned n>
VectorGeneric<n> delta(const VectorGeneric<n>&v1,const VectorGeneric<n>&v2) {
  return v2-v1;
}

template <unsigned n>
double VectorGeneric<n>::modulo2()const {
  return LoopUnroller<n>::_sum2(d.data());
}

template <unsigned n>
double dotProduct(const VectorGeneric<n>& v1,const VectorGeneric<n>& v2) {
  return LoopUnroller<n>::_dot(v1.d.data(),v2.d.data());
}

inline
VectorGeneric<3> crossProduct(const VectorGeneric<3>& v1,const VectorGeneric<3>& v2) {
  return VectorGeneric<3>(
           v1[1]*v2[2]-v1[2]*v2[1],
           v1[2]*v2[0]-v1[0]*v2[2],
           v1[0]*v2[1]-v1[1]*v2[0]);
}

template<unsigned n>
double VectorGeneric<n>::modulo()const {
  return sqrt(modulo2());
}

template<unsigned n>
double modulo2(const VectorGeneric<n>&v) {
  return v.modulo2();
}

template<unsigned n>
double modulo(const VectorGeneric<n>&v) {
  return v.modulo();
}

template<unsigned n>
std::ostream & operator<<(std::ostream &os, const VectorGeneric<n>& v) {
  for(unsigned i=0; i<n-1; i++) os<<v(i)<<" ";
  os<<v(n-1);
  return os;
}


/// \ingroup TOOLBOX
/// Alias for one dimensional vectors
typedef VectorGeneric<1> Vector1d;
/// \ingroup TOOLBOX
/// Alias for two dimensional vectors
typedef VectorGeneric<2> Vector2d;
/// \ingroup TOOLBOX
/// Alias for three dimensional vectors
typedef VectorGeneric<3> Vector3d;
/// \ingroup TOOLBOX
/// Alias for four dimensional vectors
typedef VectorGeneric<4> Vector4d;
/// \ingroup TOOLBOX
/// Alias for five dimensional vectors
typedef VectorGeneric<5> Vector5d;
/// \ingroup TOOLBOX
/// Alias for three dimensional vectors
typedef Vector3d Vector;

}

#endif

