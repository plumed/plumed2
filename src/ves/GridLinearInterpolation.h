/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The VES code team
   (see the PEOPLE-VES file at the root of this folder for a list of names)

   See http://www.ves-code.org for more information.

   This file is part of VES code module.

   The VES code module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The VES code module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the VES code module.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_ves_GridLinearInterpolation_h
#define __PLUMED_ves_GridLinearInterpolation_h

#include <vector>


namespace PLMD {


class GridBase;


namespace ves {

class GridLinearInterpolation {
private:
  static double getGridValueWithLinearInterpolation_1D(GridBase* grid_pntr, const std::vector<double>& arg);
  static double getGridValueWithLinearInterpolation_2D(GridBase* grid_pntr, const std::vector<double>& arg);
  static double linearInterpolation(const double x, const double x0, const double x1, const double y0, const double y1);
public:
  static double getGridValueWithLinearInterpolation(GridBase* grid_pntr, const std::vector<double>& arg);
};


inline
double GridLinearInterpolation::linearInterpolation(const double x, const double x0, const double x1, const double y0, const double y1) {
  // https://en.wikipedia.org/wiki/Linear_interpolation
  if(x1!=x0) {
    return y0 + (x-x0) * ((y1-y0)/(x1-x0));
  }
  else {
    return y0;
  }
}


}
}

#endif
