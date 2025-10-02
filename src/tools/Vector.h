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
#include <initializer_list>
#include <iosfwd>
#include <array>
#include <type_traits>
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
template<typename T, unsigned n> class VectorTyped;

template<typename T, unsigned n>
constexpr VectorTyped<T, n> delta(const VectorTyped<T, n>&,const VectorTyped<T, n>&);

template<typename T, unsigned n>
constexpr T dotProduct(const VectorTyped<T, n>&,const VectorTyped<T, n>&);


template<typename T, unsigned n>
std::ostream & operator<<(std::ostream &os, const VectorTyped<T, n>& v);
template<typename T, unsigned n>
class VectorTyped {
  std::array<T,n> d;
/// Auxiliary private function for constructor
  void auxiliaryConstructor();
/// Auxiliary private function for constructor
  template<typename... Args>
  void auxiliaryConstructor(T first,Args... arg);
public:
  template<unsigned m=n,typename=std::enable_if_t<m==3,void>>
  constexpr VectorTyped(T first, T second, T third):d{first,second,third} {}
/// Constructor accepting n T parameters.
/// Can be used as Vector<3>(1.0,2.0,3.0) or Vector<2>(2.0,3.0).
/// In case a wrong number of parameters is given, a static assertion will fail.
  template<typename... Args,unsigned m=n,typename=std::enable_if_t<m!=3,void>>
  VectorTyped(T first,Args... arg);
/// create it null
  constexpr VectorTyped();
/// Returns a pointer to the underlying array serving as element storage.
  constexpr T* data() noexcept;
/// Returns a pointer to the underlying array serving as element storage.
  constexpr const T* data() const noexcept;
/// set it to zero
  constexpr void zero();
/// array-like access [i]
  constexpr T & operator[](unsigned i);
/// array-like access [i]
  constexpr const T & operator[](unsigned i)const;
/// parenthesis access (i)
  constexpr T & operator()(unsigned i);
/// parenthesis access (i)
  constexpr const T & operator()(unsigned i)const;
/// increment
  constexpr VectorTyped& operator +=(const VectorTyped& b);
/// decrement
  constexpr VectorTyped& operator -=(const VectorTyped& b);
/// multiply
  constexpr VectorTyped& operator *=(T s);
/// divide
  constexpr VectorTyped& operator /=(T s);
/// sign +
  constexpr VectorTyped operator +()const;
/// sign -
  constexpr VectorTyped operator -()const;
/// return v1+v2
  template<typename U, unsigned m>
  friend constexpr VectorTyped<U, m> operator+(const VectorTyped<U, m>&,const VectorTyped<U, m>&);
/// return v1-v2
  template<typename U, unsigned m>
  friend constexpr VectorTyped<U, m> operator-(VectorTyped<U, m>,const VectorTyped<U, m>&);
/// return s*v
  template<typename U, typename J, unsigned m>
  friend constexpr VectorTyped<U, m> operator*(J,VectorTyped<U, m>);
/// return v*s
  template<typename U, typename J, unsigned m>
  friend constexpr VectorTyped<U, m> operator*(VectorTyped<U, m>,J);
/// return v/s
  template<typename U, typename J, unsigned m>
  friend constexpr VectorTyped<U, m> operator/(const VectorTyped<U, m>&,J);
/// return v2-v1
  friend constexpr VectorTyped delta<>(const VectorTyped&v1,const VectorTyped&v2);
/// return v1 .scalar. v2
  friend constexpr T dotProduct<>(const VectorTyped&,const VectorTyped&);
  //this bad boy produces a warning (in fact becasue declrare the crossproduc as a friend for ALL thhe possible combinations of n and T)
/// return v1 .vector. v2
/// Only available for size 3
  template<typename U>
  friend constexpr VectorTyped<U, 3> crossProduct(const VectorTyped<U, 3>&,const VectorTyped<U, 3>&);
/// compute the squared modulo
  constexpr T modulo2()const;
/// Compute the modulo.
/// Shortcut for sqrt(v.modulo2())
  constexpr T modulo()const;
/// friend version of modulo2 (to simplify some syntax)
  template<typename U, unsigned m>
  friend constexpr U modulo2(const VectorTyped<U, m>&);
/// friend version of modulo (to simplify some syntax)
  template<typename U, unsigned m>
  friend constexpr U modulo(const VectorTyped<U, m>&);
///conversion to another type
  template<typename TT>
  constexpr VectorTyped<TT,n> convert()const;
/// << operator.
/// Allows printing vector `v` with `std::cout<<v;`
  friend std::ostream & operator<< <>(std::ostream &os, const VectorTyped&);
};

template<typename T, unsigned n>
void VectorTyped<T, n>::auxiliaryConstructor()
{}

template<typename T, unsigned n>
template<typename... Args>
void VectorTyped<T, n>::auxiliaryConstructor(T first,Args... arg) {
  d[n-(sizeof...(Args))-1]=first;
  auxiliaryConstructor(arg...);
}

template<typename T, unsigned n>
template<typename... Args, unsigned,typename >
VectorTyped<T, n>::VectorTyped(T first,Args... arg) {
  static_assert((sizeof...(Args))+1==n,"you are trying to initialize a Vector with the wrong number of arguments");
  auxiliaryConstructor(first,arg...);
}

template<typename T, unsigned n>
constexpr T* VectorTyped<T, n>::data() noexcept {
  return d.data();
}

template<typename T, unsigned n>
constexpr const T* VectorTyped<T, n>::data() const noexcept {
  return d.data();
}

template<typename T, unsigned n>
constexpr void VectorTyped<T, n>::zero() {
  LoopUnroller<n>::_zero(d.data());
}

template<typename T, unsigned n>
constexpr VectorTyped<T, n>::VectorTyped() {
  LoopUnroller<n>::_zero(d.data());
}

template<typename T, unsigned n>
constexpr T & VectorTyped<T, n>::operator[](unsigned i) {
  return d[i];
}

template<typename T, unsigned n>
constexpr const T & VectorTyped<T, n>::operator[](unsigned i)const {
  return d[i];
}

template<typename T, unsigned n>
constexpr T & VectorTyped<T, n>::operator()(unsigned i) {
  return d[i];
}

template<typename T, unsigned n>
constexpr const T & VectorTyped<T, n>::operator()(unsigned i)const {
  return d[i];
}

template<typename T, unsigned n>
constexpr VectorTyped<T, n>& VectorTyped<T, n>::operator +=(const VectorTyped<T, n>& b) {
  LoopUnroller<n>::_add(d.data(),b.d.data());
  return *this;
}

template<typename T, unsigned n>
constexpr VectorTyped<T, n>& VectorTyped<T, n>::operator -=(const VectorTyped<T, n>& b) {
  LoopUnroller<n>::_sub(d.data(),b.d.data());
  return *this;
}

template<typename T, unsigned n>
constexpr VectorTyped<T, n>& VectorTyped<T, n>::operator *=(T s) {
  LoopUnroller<n>::_mul(d.data(),s);
  return *this;
}

template<typename T, unsigned n>
constexpr VectorTyped<T, n>& VectorTyped<T, n>::operator /=(T s) {
  LoopUnroller<n>::_mul(d.data(),1.0/s);
  return *this;
}

template<typename T, unsigned n>
constexpr VectorTyped<T, n>  VectorTyped<T, n>::operator +()const {
  return *this;
}

template<typename T, unsigned n>
constexpr VectorTyped<T, n> VectorTyped<T, n>::operator -()const {
  VectorTyped<T, n> r;
  LoopUnroller<n>::_neg(r.d.data(),d.data());
  return r;
}

template<typename T, unsigned n>
constexpr VectorTyped<T, n> operator+(const VectorTyped<T, n>&v1,const VectorTyped<T, n>&v2) {
  VectorTyped<T, n> v(v1);
  return v+=v2;
}

template<typename T, typename  TT, unsigned n>
constexpr VectorTyped<T, n> sumT(const VectorTyped<T, n>&v1,const VectorTyped<TT, n>&v2) {
  VectorTyped<T, n> v(v1);
  LoopUnroller<n>::_add(v.data(),v2.data());
  return v;
}

template<typename T, unsigned n>
constexpr VectorTyped<T, n> operator-(VectorTyped<T, n>v1,const VectorTyped<T, n>&v2) {
  return v1-=v2;
}

template<typename T, typename J, unsigned n>
constexpr VectorTyped<T, n> operator*(J s,VectorTyped<T, n>v) {
  return v*=s;
}

template<typename T, typename J, unsigned n>
constexpr VectorTyped<T, n> operator*(VectorTyped<T, n> v,J s) {
  return v*=s;
}

template<typename T, typename J, unsigned n>
constexpr VectorTyped<T, n> operator/(const VectorTyped<T, n>&v,J s) {
  return v*(T(1.0)/s);
}

template<typename T, unsigned n>
constexpr VectorTyped<T, n> delta(const VectorTyped<T, n>&v1,const VectorTyped<T, n>&v2) {
  return v2-v1;
}

template<typename T, unsigned n>
constexpr T VectorTyped<T, n>::modulo2()const {
  return LoopUnroller<n>::_sum2(d.data());
}

template<typename T, unsigned n>
constexpr T dotProduct(const VectorTyped<T, n>& v1,const VectorTyped<T, n>& v2) {
  return LoopUnroller<n>::_dot(v1.d.data(),v2.d.data());
}

template<typename T>
constexpr inline
VectorTyped<T, 3> crossProduct(const VectorTyped<T, 3>& v1,const VectorTyped<T, 3>& v2) {
  return VectorTyped<T, 3>(
           v1[1]*v2[2]-v1[2]*v2[1],
           v1[2]*v2[0]-v1[0]*v2[2],
           v1[0]*v2[1]-v1[1]*v2[0]);
}

template<typename T, unsigned n>
constexpr T VectorTyped<T, n>::modulo()const {
  return sqrt(modulo2());
}

template<typename T, unsigned n>
constexpr T modulo2(const VectorTyped<T, n>&v) {
  return v.modulo2();
}

template<typename T, unsigned n>
constexpr T modulo(const VectorTyped<T, n>&v) {
  return v.modulo();
}

template<typename T, unsigned n>
template<typename TT>
constexpr VectorTyped<TT,n> VectorTyped<T,n>::convert()const {
  VectorTyped<TT,n> target;
  LoopUnroller<n>::_copy(target.data(),d.data());
  return target;
}

template<typename T, unsigned n>
std::ostream & operator<<(std::ostream &os, const VectorTyped<T, n>& v) {
  for(unsigned i=0; i<n-1; i++) {
    os<<v(i)<<" ";
  }
  os<<v(n-1);
  return os;
}

template<unsigned n>
using VectorGeneric=VectorTyped<double, n>;

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

/// \ingroup TOOLBOX
/// Alias for one dimensional vectors
typedef VectorTyped<float,1> Vector1f;
/// \ingroup TOOLBOX
/// Alias for two dimensional vectors
typedef VectorTyped<float,2> Vector2f;
/// \ingroup TOOLBOX
/// Alias for three dimensional vectors
typedef VectorTyped<float,3> Vector3f;
/// \ingroup TOOLBOX
/// Alias for four dimensional vectors
typedef VectorTyped<float,4> Vector4f;
/// \ingroup TOOLBOX
/// Alias for five dimensional vectors
typedef VectorTyped<float,5> Vector5f;

/// \ingroup TOOLBOX
/// Alias for three dimensional vectors of any type:
template <typename T>
using VectorT=VectorTyped<T,3>;

static_assert(sizeof(VectorGeneric<2>)==2*sizeof(double), "code cannot work if this is not satisfied");
static_assert(sizeof(VectorGeneric<3>)==3*sizeof(double), "code cannot work if this is not satisfied");
static_assert(sizeof(VectorGeneric<4>)==4*sizeof(double), "code cannot work if this is not satisfied");
namespace Versors {
template <typename T>
constexpr auto xp= PLMD::VectorTyped<T,3>( {
  T( 1),T( 0),T( 0)
});
template <typename T>
constexpr auto xm= PLMD::VectorTyped<T,3>( {
  T(-1),T( 0),T( 0)
});
template <typename T>
constexpr auto yp= PLMD::VectorTyped<T,3>( {
  T( 0),T( 1),T( 0)
});
template <typename T>
constexpr auto ym= PLMD::VectorTyped<T,3>( {
  T( 0),T(-1),T( 0)
});
template <typename T>
constexpr auto zp= PLMD::VectorTyped<T,3>( {
  T( 0),T( 0),T( 1)
});
template <typename T>
constexpr auto zm= PLMD::VectorTyped<T,3>( {
  T( 0),T( 0),T(-1)
});
} // namespace Versors
} //PLMD

#endif //__PLUMED_tools_Vector_h
