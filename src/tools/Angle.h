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
#ifndef __PLUMED_tools_Angle_h
#define __PLUMED_tools_Angle_h

#include "Vector.h"
#include "Tools.h"
#include <cmath>

namespace PLMD {

/// \ingroup TOOLBOX
/// Class to compute angles.
/// I define it as a class even if it does not contain anything. The reason
/// is that in the future I would like to extend it to contain options about
/// how the calculation should be done. So, for now use it as
/// Angle a;
/// double angle=a.compute(v1,v2);
/// I know it is a bit misleading. If we really do not need to store "options"
/// inside the Angle class, we can remove it later and write compute as
/// a static function.
class Angle {
// still empty, but may accommodate some options in the future
public:
/// Compute the angle between vectors v1 and v2
  template <typename T>
  static T compute(const VectorTyped<T,3> v1,const VectorTyped<T,3> v2) {
    return std::acos(dotProduct(v1,v2)/(v1.modulo()*v2.modulo()));
  }
/// Compute the angle between vectors v1 and v2 and its derivatives wrt v1 and v2
  template <typename T>
  static T compute(const VectorTyped<T,3> v1,
                   const VectorTyped<T,3> v2,
                   VectorTyped<T,3>& d1,
                   VectorTyped<T,3>& d2) {
    using V3=VectorTyped<T,3>;
    const T dp(dotProduct(v1,v2));
    const V3& dp_dv1(v2);
    const V3& dp_dv2(v1);
    const T sv1(v1.modulo2());
    const T sv2(v2.modulo2());
    const V3 dsv1_dv1(2*v1);
    const V3 dsv2_dv2(2*v2);
    const T nn(1.0/std::sqrt(sv1*sv2));
    const V3 dnn_dv1(-0.5*nn/sv1*dsv1_dv1);
    const V3 dnn_dv2(-0.5*nn/sv2*dsv2_dv2);

    const T dpnn(dp*nn);
    if(dpnn>=1.0-epsilon) {
      d1=V3(0.0,0.0,0.0);
      d2=V3(0.0,0.0,0.0);
      return 0.0;
    }
    if(dpnn<=-1.0+epsilon) {
      d1=V3(0.0,0.0,0.0);
      d2=V3(0.0,0.0,0.0);
      return pi;
    }
    const V3 ddpnn_dv1(dp*dnn_dv1+dp_dv1*nn);
    const V3 ddpnn_dv2(dp*dnn_dv2+dp_dv2*nn);

    const T x(-1.0/std::sqrt(1-dpnn*dpnn));

    d1=x*ddpnn_dv1;
    d2=x*ddpnn_dv2;


    return std::acos(dpnn);
  }

};

}

#endif
