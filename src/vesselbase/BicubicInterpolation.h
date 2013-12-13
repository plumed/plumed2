/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
#ifndef __PLUMED_tools_BicubicInterpolation_h
#define __PLUMED_tools_BicubicInterpolation_h

#include "InterpolationBase.h"

namespace PLMD {
namespace vesselbase {

class BicubicInterpolation : public InterpolationBase {
private:
  Matrix<double> wt, dcross;
  std::vector<double> t1, t2, clist;
/// Build the elements of the interpolation thing
  void IBicCoeff( const std::vector<double>& y, const std::vector<double>& dy1, const std::vector<double>& dy2,
                  const std::vector<double>& d2y12, const double& d1, const double& d2, Matrix<double>& c );
public:
  BicubicInterpolation( GridVesselBase* gg, const unsigned dstart=0 );
/// Setup interpolation tables
  void setInterpolationTables();
/// Get the value of the function at point pos
  double interpolateFunction( const unsigned& mybox, const std::vector<double>& dd );
};

}
}

#endif
