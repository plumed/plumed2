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

#include "GridLinearInterpolation.h"

#include "tools/Grid.h"
#include "tools/Exception.h"
#include "tools/Tools.h"


namespace PLMD {
namespace ves {


double GridLinearInterpolation::getGridValueWithLinearInterpolation_1D(GridBase* grid_pntr, const std::vector<double>& arg) {

  plumed_massert(grid_pntr->getDimension()==1,"The grid is of the wrong dimension, should be one-dimensional");
  plumed_massert(arg.size()==1,"input value is of the wrong size");

  double x = arg[0];
  double grid_dx = grid_pntr->getDx()[0];
  double grid_min; Tools::convert( grid_pntr->getMin()[0], grid_min);
  std::vector<unsigned int> i0(1); i0[0] = unsigned( std::floor( (x-grid_min)/grid_dx ) );
  std::vector<unsigned int> i1(1); i1[0] = unsigned( std::ceil(  (x-grid_min)/grid_dx ) );
  //
  double x0 = grid_pntr->getPoint(i0)[0];
  double x1 = grid_pntr->getPoint(i1)[0];
  double y0 = grid_pntr->getValue(i0);
  double y1 = grid_pntr->getValue(i1);
  //
  return linearInterpolation(x,x0,x1,y0,y1);
}


double GridLinearInterpolation::getGridValueWithLinearInterpolation_2D(GridBase* grid_pntr, const std::vector<double>& arg) {
  plumed_massert(grid_pntr->getDimension()==2,"The grid is of the wrong dimension, should be two-dimensional");
  plumed_massert(arg.size()==2,"input value is of the wrong size");

  std::vector<double> grid_delta = grid_pntr->getDx();
  std::vector<double> grid_min(2);
  Tools::convert( grid_pntr->getMin()[0], grid_min[0]);
  Tools::convert( grid_pntr->getMin()[1], grid_min[1]);

  std::vector<unsigned int> i00(2);
  std::vector<unsigned int> i01(2);
  std::vector<unsigned int> i10(2);
  std::vector<unsigned int> i11(2);

  i00[0] = i01[0] = unsigned( std::floor( (arg[0]-grid_min[0])/grid_delta[0] ) );
  i10[0] = i11[0] = unsigned( std::ceil(  (arg[0]-grid_min[0])/grid_delta[0] ) );

  i00[1] = i10[1] = unsigned( std::floor( (arg[1]-grid_min[1])/grid_delta[1] ) );
  i01[1] = i11[1] = unsigned( std::ceil(  (arg[1]-grid_min[1])/grid_delta[1] ) );

  // https://en.wikipedia.org/wiki/Bilinear_interpolation
  double x = arg[0];
  double y = arg[1];

  double x0 = grid_pntr->getPoint(i00)[0];
  double x1 = grid_pntr->getPoint(i10)[0];
  double y0 = grid_pntr->getPoint(i00)[1];
  double y1 = grid_pntr->getPoint(i11)[1];

  double f00 = grid_pntr->getValue(i00);
  double f01 = grid_pntr->getValue(i01);
  double f10 = grid_pntr->getValue(i10);
  double f11 = grid_pntr->getValue(i11);

  // linear interpolation in x-direction
  double fx0 = linearInterpolation(x,x0,x1,f00,f10);
  double fx1 = linearInterpolation(x,x0,x1,f01,f11);
  // linear interpolation in y-direction
  double fxy = linearInterpolation(y,y0,y1,fx0,fx1);
  //
  return fxy;
}


double GridLinearInterpolation::getGridValueWithLinearInterpolation(GridBase* grid_pntr, const std::vector<double>& arg) {
  unsigned int dim = grid_pntr->getDimension();
  if(dim==1) {
    return getGridValueWithLinearInterpolation_1D(grid_pntr,arg);
  }
  else if(dim==2) {
    return getGridValueWithLinearInterpolation_2D(grid_pntr,arg);
  }
  else {
    plumed_merror("GridLinearInterpolation::getGridValueWithLinearInterpolation only defined for 1D or 2D grids");
  }
}



}
}
