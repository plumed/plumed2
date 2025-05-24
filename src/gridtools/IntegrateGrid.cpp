/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The plumed team
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
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"

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

namespace PLMD {
namespace gridtools {

class IntegrateGrid : public ActionShortcut {
public:
  static void registerKeywords(Keywords&);
  explicit IntegrateGrid(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(IntegrateGrid,"INTEGRATE_GRID")

void IntegrateGrid::registerKeywords(Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","grid","the label of the function on a grid that is being integrated");
  keys.add("compulsory","PERIODIC","if the value of the output integral has a periodid domain then you use this keyword to specify the periodicity.  If the output integral is not periodic you must state this using PERIODIC=NO");
  keys.setValueDescription("scalar","the numerical integral of the input function over its whole domain");
  keys.needsAction("GET_VOLUME_ELEMENT");
  keys.needsAction("SUM");
  keys.needsAction("CUSTOM");
}

IntegrateGrid::IntegrateGrid(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string gridname;
  parse("ARG", gridname);
  readInputLine( getShortcutLabel() + "_volelem: GET_VOLUME_ELEMENT ARG=" + gridname );
  readInputLine( getShortcutLabel() + "_sum: SUM ARG=" + gridname + " PERIODIC=NO" );
  std::string period;
  parse("PERIODIC", period);
  readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_volelem," + getShortcutLabel() + "_sum FUNC=x*y PERIODIC=" + period );
}


}
}
