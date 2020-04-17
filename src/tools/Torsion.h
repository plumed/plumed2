/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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

#include "Vector.h"

namespace PLMD {

/// \ingroup TOOLBOX
/// Class to compute torsional angles.
/// I define it as a class even if it does not contain anything. The reason
/// is that in the future I would like to extend it to contain options about
/// how the calculation should be done. So, for now use it as
/// Torsion t;
/// double angle=t.compute(v1,v2,v3);
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
  double compute(const Vector& v1,const Vector& v2,const Vector& v3)const;
/// This is the version which also computes the derivatives wrt the arguments.
  double compute(const Vector& v1,const Vector& v2,const Vector& v3,Vector& d1,Vector& d2,Vector& d3)const;
};

}

#endif
