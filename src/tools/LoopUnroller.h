/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
template<unsigned n>
class LoopUnroller{
public:
  static void _zero(double*);
  static void _add(double*,const double*);
  static void _sub(double*,const double*);
  static void _mul(double*,const double);
  static void _neg(double*,const double*);
  static double _sum2(const double*);
  static double _dot(const double*,const double*);
};

template<unsigned n>
void LoopUnroller<n>::_zero(double*d){
  LoopUnroller<n-1>::_zero(d);
  d[n-1]=0.0;
}

template<>
inline
void LoopUnroller<1>::_zero(double*d){
  d[0]=0.0;
}

template<unsigned n>
void LoopUnroller<n>::_add(double*d,const double*a){
  LoopUnroller<n-1>::_add(d,a);
  d[n-1]+=a[n-1];
}

template<>
inline
void LoopUnroller<1>::_add(double*d,const double*a){
  d[0]+=a[0];
}

template<unsigned n>
void LoopUnroller<n>::_sub(double*d,const double*a){
  LoopUnroller<n-1>::_sub(d,a);
  d[n-1]-=a[n-1];
}

template<>
inline
void LoopUnroller<1>::_sub(double*d,const double*a){
  d[0]-=a[0];
}

template<unsigned n>
void LoopUnroller<n>::_mul(double*d,const double s){
  LoopUnroller<n-1>::_mul(d,s);
  d[n-1]*=s;
}

template<>
inline
void LoopUnroller<1>::_mul(double*d,const double s){
  d[0]*=s;
}

template<unsigned n>
void LoopUnroller<n>::_neg(double*d,const double*a ){
  LoopUnroller<n-1>::_neg(d,a);
  d[n-1]=-a[n-1];
}

template<>
inline
void LoopUnroller<1>::_neg(double*d,const double*a){
  d[0]=-a[0];
}

template<unsigned n>
double LoopUnroller<n>::_sum2(const double*d){
  return LoopUnroller<n-1>::_sum2(d)+d[n-1]*d[n-1];
}

template<>
inline
double LoopUnroller<1>::_sum2(const double*d){
  return d[0]*d[0];
}

template<unsigned n>
double LoopUnroller<n>::_dot(const double*d,const double*v){
  return LoopUnroller<n-1>::_dot(d,v)+d[n-1]*v[n-1];
}

template<>
inline
double LoopUnroller<1>::_dot(const double*d,const double*v){
  return d[0]*v[0];
}

#endif
