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
#ifndef __PLUMED_tools_Torsion_h
#define __PLUMED_tools_Torsion_h

#include "Tensor.h"

namespace PLMD {

/// \ingroup TOOLBOX
/// Class to compute torsional angles.
/// I define it as a class even if it does not contain anything. The reason
/// is that in the future I would like to extend it to contain options about
/// how the calculation should be done. So, for now use it as
/// Torsion t;
/// T angle=t.compute(v1,v2,v3);
/// I know it is a bit misleading. If we really do not need to store "options"
/// inside the Torsion class, we can remove it later and write compute as
/// a static function.
class Torsion {
// still empty, but may accommodate some options in the future
public:
/// Compute the angle between the projections of v1 and v3 on the plane orthogonal
/// to v2. To have a "normal" definition (= plumed1), use it as
/// compute(r01,r12,r23);
/// See ColvarTorsion for a practical usage...
  template <typename T>
  static T compute(const VectorTyped<T,3>& v1,
                   const VectorTyped<T,3>& v2,
                   const VectorTyped<T,3>& v3) {
    using V3=PLMD::VectorTyped<T,3>;

    const V3 nv2(v2*(1.0/v2.modulo()));
    const V3 a(crossProduct(nv2,v1));
    const V3 b(crossProduct(v3,nv2));
    const T cosangle=dotProduct(a,b);
    const T sinangle=dotProduct(crossProduct(a,b),nv2);
    return std::atan2(-sinangle,cosangle);
  }

/// This is the version which also computes the derivatives wrt the arguments.
  template <typename T>
  static T compute(const VectorTyped<T,3>& v1,
                   const VectorTyped<T,3>& v2,
                   const VectorTyped<T,3>& v3,
                   VectorTyped<T,3>& d1,
                   VectorTyped<T,3>& d2,
                   VectorTyped<T,3>& d3) {
    using V3=PLMD::VectorTyped<T,3>;
    using T33=PLMD::TensorTyped<T,3,3>;

    const T modv2(T(1.0)/v2.modulo());
    const V3 nv2(v2*modv2);
    const T33 dnv2_v2((T33::identity()-extProduct(nv2,nv2))*modv2);

    const V3 a(crossProduct(v2,v1));
    const T33 da_dv2(dcrossDv1(v2,v1));
    const T33 da_dv1(dcrossDv2(v2,v1));
    const V3 b(crossProduct(v3,v2));
    const T33 db_dv3(dcrossDv1(v3,v2));
    const T33 db_dv2(dcrossDv2(v3,v2));
    const T cosangle=dotProduct(a,b);
    const V3 dcosangle_dv1=matmul(b,da_dv1);
    const V3 dcosangle_dv2=matmul(b,da_dv2) + matmul(a,db_dv2);
    const V3 dcosangle_dv3=matmul(a,db_dv3);

    const V3 cab(crossProduct(a,b));
    const T33 dcab_dv1(matmul(dcrossDv1(a,b),da_dv1));
    const T33 dcab_dv2(matmul(dcrossDv1(a,b),da_dv2) + matmul(dcrossDv2(a,b),db_dv2));
    const T33 dcab_dv3(matmul(dcrossDv2(a,b),db_dv3));

    const T sinangle=dotProduct(cab,nv2);
    const V3 dsinangle_dv1=matmul(nv2,dcab_dv1);
    const V3 dsinangle_dv2=matmul(nv2,dcab_dv2)+matmul(cab,dnv2_v2);
    const V3 dsinangle_dv3=matmul(nv2,dcab_dv3);

    const T torsion=std::atan2(-sinangle,cosangle);
// this is required since v1 and v3 are not normalized:
    const T invR2=T(1.0)/(cosangle*cosangle+sinangle*sinangle);

    d1= ( -dsinangle_dv1*cosangle + sinangle * dcosangle_dv1 ) *invR2;
    d2= ( -dsinangle_dv2*cosangle + sinangle * dcosangle_dv2 ) *invR2;
    d3= ( -dsinangle_dv3*cosangle + sinangle * dcosangle_dv3 ) *invR2;

    return torsion;
  }

};

}

#endif
