/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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
template<unsigned n>
class LoopUnroller {
public:
/// Set to zero.
/// Same as `for(unsigned i=0;i<n;i++) d[i]=0.0;`
  static void _zero(double*d);
/// Add v to d.
/// Same as `for(unsigned i=0;i<n;i++) d[i]+=v[i];`
  static void _add(double*d,const double*v);
/// Subtract v from d.
/// Same as `for(unsigned i=0;i<n;i++) d[i]-=v[i];`
  static void _sub(double*d,const double*v);
/// Multiply d by s.
/// Same as `for(unsigned i=0;i<n;i++) d[i]*=s;`
  static void _mul(double*d,const double s);
/// Set d to -v.
/// Same as `for(unsigned i=0;i<n;i++) d[i]=-v[i];`
  static void _neg(double*d,const double*v);
/// Squared modulo of d;
/// Same as `r=0.0; for(unsigned i=0;i<n;i++) r+=d[i]*d[i]; return r;`
  static double _sum2(const double*d);
/// Dot product of d and v
/// Same as `r=0.0; for(unsigned i=0;i<n;i++) r+=d[i]*v[i]; return r;`
  static double _dot(const double*d,const double*v);
};

template<unsigned n>
void LoopUnroller<n>::_zero(double*d) {
  LoopUnroller<n-1>::_zero(d);
  d[n-1]=0.0;
}

template<>
inline
void LoopUnroller<1>::_zero(double*d) {
  d[0]=0.0;
}

template<unsigned n>
void LoopUnroller<n>::_add(double*d,const double*a) {
  LoopUnroller<n-1>::_add(d,a);
  d[n-1]+=a[n-1];
}

template<>
inline
void LoopUnroller<1>::_add(double*d,const double*a) {
  d[0]+=a[0];
}

template<unsigned n>
void LoopUnroller<n>::_sub(double*d,const double*a) {
  LoopUnroller<n-1>::_sub(d,a);
  d[n-1]-=a[n-1];
}

template<>
inline
void LoopUnroller<1>::_sub(double*d,const double*a) {
  d[0]-=a[0];
}

template<unsigned n>
void LoopUnroller<n>::_mul(double*d,const double s) {
  LoopUnroller<n-1>::_mul(d,s);
  d[n-1]*=s;
}

template<>
inline
void LoopUnroller<1>::_mul(double*d,const double s) {
  d[0]*=s;
}

template<unsigned n>
void LoopUnroller<n>::_neg(double*d,const double*a ) {
  LoopUnroller<n-1>::_neg(d,a);
  d[n-1]=-a[n-1];
}

template<>
inline
void LoopUnroller<1>::_neg(double*d,const double*a) {
  d[0]=-a[0];
}

template<unsigned n>
double LoopUnroller<n>::_sum2(const double*d) {
  return LoopUnroller<n-1>::_sum2(d)+d[n-1]*d[n-1];
}

template<>
inline
double LoopUnroller<1>::_sum2(const double*d) {
  return d[0]*d[0];
}

template<unsigned n>
double LoopUnroller<n>::_dot(const double*d,const double*v) {
  return LoopUnroller<n-1>::_dot(d,v)+d[n-1]*v[n-1];
}

template<>
inline
double LoopUnroller<1>::_dot(const double*d,const double*v) {
  return d[0]*v[0];
}

}

#endif
