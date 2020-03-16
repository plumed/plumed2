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
#ifndef __PLUMED_tools_Angle_h
#define __PLUMED_tools_Angle_h

#include "Vector.h"

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
  double compute(const Vector& v1,const Vector& v2)const;
/// Compute the angle between vectors v1 and v2 and its derivatives wrt v1 and v2
  double compute(const Vector& v1,const Vector& v2,Vector& d1,Vector& d2)const;
};

}

#endif
