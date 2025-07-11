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
#ifndef __PLUMED_gridtools_Interpolator_h
#define __PLUMED_gridtools_Interpolator_h

#include "core/Value.h"
#include "GridCoordinatesObject.h"

namespace PLMD {
namespace gridtools {

class Interpolator {
private:
  Value* values;
  const GridCoordinatesObject& gridobject;

public:
  Interpolator( Value* myval, const GridCoordinatesObject& mygrid ) : values(myval), gridobject(mygrid) {}
  /// Interpolate the function using splines
  double splineInterpolation( const std::vector<double>& x, std::vector<double>& der ) const ;
  double splineInterpolation( View<const double> x, std::vector<double>& der ) const ;
};

}
}
#endif
