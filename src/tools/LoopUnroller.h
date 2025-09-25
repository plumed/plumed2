/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2023 The plumed team
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
#ifndef __PLUMED_tools_LoopUnroller_h
#define __PLUMED_tools_LoopUnroller_h

namespace PLMD {

/**
\ingroup TOOLBOX
Utiliy class for loop unrolling.

Many c++ compilers do not unroll small loops such as those
used in the PLMD::Vector and PLMD::Tensor classes.
This class provides methods to perform basic vector
operations with unrolled loops. The methods work on double*
so that they can be used in principles in other places of the code,
but they are designed to be used in PLMD::Vector and PLMD::Tensor .

In case in the future we see that some compiler better optimize explicit loops,
it should be easy to replace the methods here with loops. Alternatively,
we could provide two paths using a cpp macro (e.g. __PLUMED_UNROLL_LOOPS or so).

All the methods for class LoopUnroller<n> act on n elements.
Implementation is made using template metaprogramming, that is:
- LoopUnroller<1>::xxx acts on the element [0] of the array.
- LoopUnroller<n>::xxx calls LoopUnroller<n-1>::xxx then acts on element [n-1] of the array.

Here xxx is any of the methods of the class.

*/
//the typename is the second parameter argument so that it can be deduced
template<unsigned n>
class LoopUnroller {
public:
/// Set to zero.
/// Same as `for(unsigned i=0;i<n;i++) d[i]=0.0;`
  template<typename T>
  static constexpr void _zero(T*d) noexcept;
/// copy the content o v into v
/// Same as `for(unsigned i=0;i<n;i++) d[i]=v[i];`
  template<typename T,typename TT>
  static constexpr void _copy(T*d,TT const*v) noexcept;
/// Add v to d.
/// Same as `for(unsigned i=0;i<n;i++) d[i]+=v[i];`
  template<typename T,typename TT>
  static constexpr void _add(T*d,TT const*v) noexcept;
/// Subtract v from d.
/// Same as `for(unsigned i=0;i<n;i++) d[i]-=v[i];`
  template<typename T,typename TT>
  static constexpr void _sub(T*d,TT const*v) noexcept;
/// Multiply d by s.
/// Same as `for(unsigned i=0;i<n;i++) d[i]*=s;`
  template<typename T,typename TT>
  static constexpr void _mul(T*d,TT const s) noexcept;
/// Set d to -v.
/// Same as `for(unsigned i=0;i<n;i++) d[i]=-v[i];`
  template<typename T,typename TT>
  static constexpr void _neg(T*d,TT const*v) noexcept;
/// Squared modulo of d;
/// Same as `r=0.0; for(unsigned i=0;i<n;i++) r+=d[i]*d[i]; return r;`
  template<typename T>
  static constexpr T _sum2(T const*d) noexcept;
/// Dot product of d and v
/// Same as `r=0.0; for(unsigned i=0;i<n;i++) r+=d[i]*v[i]; return r;`
  template<typename T,typename TT>
  static constexpr T _dot(T const*d,TT const*v) noexcept;
};

template<unsigned n>
template<typename T>
constexpr void LoopUnroller<n>::_zero(T*d) noexcept {
  if constexpr (n>1) {
    LoopUnroller<n-1>::_zero(d);
  }
  d[n-1]=T(0);
}

template<unsigned n>
template<typename T,typename TT>
constexpr void LoopUnroller<n>::_copy(T*d,TT const*v) noexcept {
  if constexpr (n==1) {
    d[0]=v[0];
  } else {
    LoopUnroller<n-1>::_copy(d,v);
    d[n-1]=v[n-1];
  }
}

template<unsigned n>
template<typename T, typename TT>
constexpr void LoopUnroller<n>::_add(T*d,TT const*a) noexcept {
  if constexpr (n>1) {
    LoopUnroller<n-1>::_add(d,a);
  }
  d[n-1]+=a[n-1];
}

template<unsigned n>
template<typename T, typename TT>
constexpr void LoopUnroller<n>::_sub(T*d,TT const*a) noexcept {
  if constexpr (n>1) {
    LoopUnroller<n-1>::_sub(d,a);
  }
  d[n-1]-=a[n-1];
}

template<unsigned n>
template<typename T, typename TT>
constexpr void LoopUnroller<n>::_mul(T*d,TT const s) noexcept {
  if constexpr (n>1) {
    LoopUnroller<n-1>::_mul(d,s);
  }
  d[n-1]*=s;
}

template<unsigned n>
template<typename T, typename TT>
constexpr void LoopUnroller<n>::_neg(T*d,TT const*a ) noexcept {
  if constexpr (n>1) {
    LoopUnroller<n-1>::_neg(d,a);
  }
  d[n-1]=-a[n-1];
}

template<unsigned n>
template<typename T>
constexpr T LoopUnroller<n>::_sum2(T const*d) noexcept {
  if constexpr (n>1) {
    return LoopUnroller<n-1>::_sum2(d)+d[n-1]*d[n-1];
  } else {
    return d[0]*d[0];
  }
}
template<unsigned n>
template<typename T, typename TT>
constexpr T LoopUnroller<n>::_dot(T const*d,TT const*v) noexcept {
  if constexpr (n>1) {
    return LoopUnroller<n-1>::_dot(d,v)+d[n-1]*v[n-1];
  } else {
    return d[0]*v[0];
  }
}
} //PLMD
#endif //__PLUMED_tools_LoopUnroller_h
