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

This action takes a function, $f(x)$, that is stored on a grid of points between $a$ and $b$.  It then computes the following definite integral
numerically:

$$
y = \int_a^b f(x) \textrm{d}x
$$

The approximate value of this integral would be computed using:

$$
y = \frac{b-a}{n} \sum_{i=0}^{n} f\left( a + \frac{i(b-a)}{n} \right)
$$

where $n$ is the number of points at which the grid has been computed.

The following example input demonstrates how this action can be used:

```plumed
c1: COORDINATIONNUMBER SPECIESA=1-10 SPECIESB=1-200 SWITCH={RATIONAL R_0=1.0}
hu: KDE ARG=c1 GRID_BIN=200 GRID_MIN=0 GRID_MAX=10 BANDWIDTH=0.1
h: CUSTOM ARG=hu FUNC=x/10 PERIODIC=NO
iv: INTEGRATE_GRID ARG=h PERIODIC=NO
```

In this case, the [KDE](KDE.md) action is used to compute a histogram that shows the instanenous distribution for 10 coordination numbers. This distribution is
unormalised so that [CUSTOM](CUSTOM.md) command here is used to normalise the distribution.  The resulting integral that is computed here should thsu be one.

Notice also that the action can still be used if the function stored on the grid is a function of 2 or more variables.  For example, the following input performs an integral
that is very similar to the one performed in the previous.  Now, however, the input grid is two dimensional.

```plumed
c1: COORDINATIONNUMBER SPECIESA=1-10 SPECIESB=1-200 SWITCH={RATIONAL R_0=1.0}
q6: Q6 SPECIESA=1-10 SPECIESB=1-200 SWITCH={RATIONAL R_0=1.0}
hu: KDE ARG=c1,q6 GRID_BIN=200,200 GRID_MIN=0,-0.5 GRID_MAX=10,2 BANDWIDTH=0.1,0.1
h: CUSTOM ARG=hu FUNC=x/10 PERIODIC=NO
iv: INTEGRATE_GRID ARG=h PERIODIC=NO
```

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
