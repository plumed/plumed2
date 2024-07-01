/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
  VectorGeneric<3,double> v1;
  v1[0]=3.0;
// use equivalently () and [] syntax:
  v1(1)=5.0;
// initialize with components
  VectorGeneric<3,double> v2=VectorGeneric<3,double>(1.0,2.0,3.0);
  VectorGeneric<3,double> v3=crossProduct(v1,v2);
  double d=dotProduct(v1,v2);
  v3+=v1;
  v2=v1+2.0*v3;
}
\endverbatim

*/
template<typename T, unsigned n> class VectorGeneric;
template<typename T, unsigned n>
std::ostream & operator<<(std::ostream &os, const VectorGeneric<T, n>& v);

template<typename T, unsigned n>
VectorGeneric<T, n> delta(const VectorGeneric<T, n>&,const VectorGeneric<T, n>&);

template<typename T, unsigned n>
class VectorGeneric {
  std::array<T,n> d;
/// Auxiliary private function for constructor
  void auxiliaryConstructor();
/// Auxiliary private function for constructor
  template<typename... Args>
  void auxiliaryConstructor(T first,Args... arg);
public:
/// Constructor accepting n T parameters.
/// Can be used as Vector<3>(1.0,2.0,3.0) or Vector<2>(2.0,3.0).
/// In case a wrong number of parameters is given, a static assertion will fail.
  template<typename... Args>
  VectorGeneric(T first,Args... arg);
/// create it null
  VectorGeneric();
/// get the underline pointer to the data
  T* data();
/// get the underline pointer to the data
  const T* data() const;
/// set it to zero
  void zero();
/// array-like access [i]
  T & operator[](unsigned i);
/// array-like access [i]
  const T & operator[](unsigned i)const;
/// parenthesis access (i)
  T & operator()(unsigned i);
/// parenthesis access (i)
  const T & operator()(unsigned i)const;
/// increment
  VectorGeneric& operator +=(const VectorGeneric& b);
/// decrement
  VectorGeneric& operator -=(const VectorGeneric& b);
/// multiply
  VectorGeneric& operator *=(T s);
/// divide
  VectorGeneric& operator /=(T s);
/// sign +
  VectorGeneric operator +()const;
/// sign -
  VectorGeneric operator -()const;
/// return v1+v2
  template<typename U, unsigned m>
  friend VectorGeneric<U, m> operator+(const VectorGeneric<U, m>&,const VectorGeneric<U, m>&);
/// return v1-v2
  template<typename U, unsigned m>
  friend VectorGeneric<U, m> operator-(VectorGeneric<U, m>,const VectorGeneric<U, m>&);
/// return s*v
  template<typename U, typename J, unsigned m>
  friend VectorGeneric<U, m> operator*(J,VectorGeneric<U, m>);
/// return v*s
  template<typename U, typename J, unsigned m>
  friend VectorGeneric<U, m> operator*(VectorGeneric<U, m>,J);
/// return v/s
  template<typename U, typename J, unsigned m>
  friend VectorGeneric<U, m> operator/(const VectorGeneric<U, m>&,J);
/// return v2-v1
  friend VectorGeneric delta<>(const VectorGeneric&v1,const VectorGeneric&v2);
/// return v1 .scalar. v2
  template<typename U, unsigned m>
  friend U dotProduct(const VectorGeneric<U, m>&,const VectorGeneric<U, m>&);
  //this bad boy produces a warning (in fact becasue declrare the crossproduc as a friend for ALL thhe possible combinations of n and T)
/// return v1 .vector. v2
/// Only available for size 3
  template<typename U>
  friend VectorGeneric<U, 3> crossProduct(const VectorGeneric<U, 3>&,const VectorGeneric<U, 3>&);
/// compute the squared modulo
  T modulo2()const;
/// Compute the modulo.
/// Shortcut for sqrt(v.modulo2())
  T modulo()const;
/// friend version of modulo2 (to simplify some syntax)
  template<typename U, unsigned m>
  friend U modulo2(const VectorGeneric<U, m>&);
/// friend version of modulo (to simplify some syntax)
  template<typename U, unsigned m>
  friend U modulo(const VectorGeneric<U, m>&);
/// << operator.
/// Allows printing vector `v` with `std::cout<<v;`
  friend std::ostream & operator<< <> (std::ostream &os, const VectorGeneric&);
};

template<typename T, unsigned n>
void VectorGeneric<T, n>::auxiliaryConstructor()
{}

template<typename T, unsigned n>
template<typename... Args>
void VectorGeneric<T, n>::auxiliaryConstructor(T first,Args... arg)
{
  d[n-(sizeof...(Args))-1]=first;
  auxiliaryConstructor(arg...);
}

template<typename T, unsigned n>
template<typename... Args>
VectorGeneric<T, n>::VectorGeneric(T first,Args... arg)
{
  static_assert((sizeof...(Args))+1==n,"you are trying to initialize a Vector with the wrong number of arguments");
  auxiliaryConstructor(first,arg...);
}

template<typename T, unsigned n>
T* VectorGeneric<T, n>::data() {return d.data();}

template<typename T, unsigned n>
const T* VectorGeneric<T, n>::data() const {return d.data();}

template<typename T, unsigned n>
void VectorGeneric<T, n>::zero() {
  LoopUnroller<T, n>::_zero(d.data());
}

template<typename T, unsigned n>
VectorGeneric<T, n>::VectorGeneric() {
  LoopUnroller<T, n>::_zero(d.data());
}

template<typename T, unsigned n>
T & VectorGeneric<T, n>::operator[](unsigned i) {
  return d[i];
}

template<typename T, unsigned n>
const T & VectorGeneric<T, n>::operator[](unsigned i)const {
  return d[i];
}

template<typename T, unsigned n>
T & VectorGeneric<T, n>::operator()(unsigned i) {
  return d[i];
}

template<typename T, unsigned n>
const T & VectorGeneric<T, n>::operator()(unsigned i)const {
  return d[i];
}

template<typename T, unsigned n>
VectorGeneric<T, n>& VectorGeneric<T, n>::operator +=(const VectorGeneric<T, n>& b) {
  LoopUnroller<T, n>::_add(d.data(),b.d.data());
  return *this;
}

template<typename T, unsigned n>
VectorGeneric<T, n>& VectorGeneric<T, n>::operator -=(const VectorGeneric<T, n>& b) {
  LoopUnroller<T, n>::_sub(d.data(),b.d.data());
  return *this;
}

template<typename T, unsigned n>
VectorGeneric<T, n>& VectorGeneric<T, n>::operator *=(T s) {
  LoopUnroller<T, n>::_mul(d.data(),s);
  return *this;
}

template<typename T, unsigned n>
VectorGeneric<T, n>& VectorGeneric<T, n>::operator /=(T s) {
  LoopUnroller<T, n>::_mul(d.data(),1.0/s);
  return *this;
}

template<typename T, unsigned n>
VectorGeneric<T, n>  VectorGeneric<T, n>::operator +()const {
  return *this;
}

template<typename T, unsigned n>
VectorGeneric<T, n> VectorGeneric<T, n>::operator -()const {
  VectorGeneric<T, n> r;
  LoopUnroller<T, n>::_neg(r.d.data(),d.data());
  return r;
}

template<typename T, unsigned n>
VectorGeneric<T, n> operator+(const VectorGeneric<T, n>&v1,const VectorGeneric<T, n>&v2) {
  VectorGeneric<T, n> v(v1);
  return v+=v2;
}

template<typename T, unsigned n>
VectorGeneric<T, n> operator-(VectorGeneric<T, n>v1,const VectorGeneric<T, n>&v2) {
  return v1-=v2;
}

template<typename T, typename J, unsigned n>
VectorGeneric<T, n> operator*(J s,VectorGeneric<T, n>v) {
  return v*=s;
}

template<typename T, typename J, unsigned n>
VectorGeneric<T, n> operator*(VectorGeneric<T, n> v,J s) {
  return v*=s;
}

template<typename T, typename J, unsigned n>
VectorGeneric<T, n> operator/(const VectorGeneric<T, n>&v,J s) {
  return v*(T(1.0)/s);
}

template<typename T, unsigned n>
VectorGeneric<T, n> delta(const VectorGeneric<T, n>&v1,const VectorGeneric<T, n>&v2) {
  return v2-v1;
}

template<typename T, unsigned n>
T VectorGeneric<T, n>::modulo2()const {
  return LoopUnroller<T, n>::_sum2(d.data());
}

template<typename T, unsigned n>
T dotProduct(const VectorGeneric<T, n>& v1,const VectorGeneric<T, n>& v2) {
  return LoopUnroller<T, n>::_dot(v1.d.data(),v2.d.data());
}

template<typename T>
inline
VectorGeneric<T, 3> crossProduct(const VectorGeneric<T, 3>& v1,const VectorGeneric<T, 3>& v2) {
  return VectorGeneric<T, 3>(
           v1[1]*v2[2]-v1[2]*v2[1],
           v1[2]*v2[0]-v1[0]*v2[2],
           v1[0]*v2[1]-v1[1]*v2[0]);
}

template<typename T, unsigned n>
T VectorGeneric<T, n>::modulo()const {
  return sqrt(modulo2());
}

template<typename T, unsigned n>
T modulo2(const VectorGeneric<T, n>&v) {
  return v.modulo2();
}

template<typename T, unsigned n>
T modulo(const VectorGeneric<T, n>&v) {
  return v.modulo();
}

template<typename T, unsigned n>
std::ostream & operator<<(std::ostream &os, const VectorGeneric<T, n>& v) {
  for(unsigned i=0; i<n-1; i++) os<<v(i)<<" ";
  os<<v(n-1);
  return os;
}


/// \ingroup TOOLBOX
/// Alias for one dimensional vectors
typedef VectorGeneric<double, 1> Vector1d;
/// \ingroup TOOLBOX
/// Alias for two dimensional vectors
typedef VectorGeneric<double, 2> Vector2d;
/// \ingroup TOOLBOX
/// Alias for three dimensional vectors
// typedef VectorGeneric<double, 3> Vector3d;
using Vector3d=VectorGeneric<double, 3>;
/// \ingroup TOOLBOX
/// Alias for four dimensional vectors
typedef VectorGeneric<double, 4> Vector4d;
/// \ingroup TOOLBOX
/// Alias for five dimensional vectors
typedef VectorGeneric<double, 5> Vector5d;
/// \ingroup TOOLBOX
/// Alias for three dimensional vectors
using Vector=Vector3d;
//using the using keyword seems to be more trasparent

static_assert(sizeof(VectorGeneric<double, 2>)==2*sizeof(double), "code cannot work if this is not satisfied");
static_assert(sizeof(VectorGeneric<double, 3>)==3*sizeof(double), "code cannot work if this is not satisfied");
static_assert(sizeof(VectorGeneric<double, 4>)==4*sizeof(double), "code cannot work if this is not satisfied");

} //PLMD

#endif //__PLUMED_tools_Vector_h
