/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2017 The plumed team
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
#include "FunctionOfGrid.h"
#include "core/ActionRegister.h"
#include "function/Sum.h"
#include "function/Custom.h"

//+PLUMEDOC GRIDCALC SUM_GRID
/*
Sum the values of all the function at the points on a grid

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC GRIDCALC INTEGRATE_GRID
/*
Calculate the numerical integral of the function stored on the grid

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC GRIDCALC CUSTOM_GRID
/*
Calculate a function of the grid or grids that are input and return a new grid

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC GRIDCALC MATHEVAL_GRID
/*
Calculate a function of the grid or grids that are input and return a new grid

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

typedef FunctionOfGrid<function::Sum> GridSum;
PLUMED_REGISTER_ACTION(GridSum,"SUM_GRID")
PLUMED_REGISTER_ACTION(GridSum,"INTEGRATE_GRID")
typedef FunctionOfGrid<function::Custom> GridCustom;
PLUMED_REGISTER_ACTION(GridCustom,"CUSTOM_GRID")
PLUMED_REGISTER_ACTION(GridCustom,"MATHEVAL_GRID")

}
}
