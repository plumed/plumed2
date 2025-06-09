/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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
#include "core/PlumedMain.h"

//+PLUMEDOC GRIDANALYSIS CONVERT_TO_FES
/*
Convert a histogram to a free energy surface.

This action allows you to take a free energy surface that was calculated using the [HISTOGRAM](HISTOGRAM.md)
action and to convert it to a free energy surface.  This transformation performed by doing:

$$
F(x) = -k_B T \ln H(x)
$$

The free energy calculated on a grid is output by this action and can be printed using [DUMPGRID](DUMPGRID.md)

## Examples

This is a typical example showing how CONVERT_TO_FES might be used when post processing a trajectory.
The input below calculates the free energy as a function of the distance between atom 1 and atom 2.
This is done by accumulating a histogram as a function of this distance using kernel density estimation
and the HISTOGRAM action.  All the data within this trajectory is used in the construction of this
HISTOGRAM.  Finally, once all the data has been read in, the histogram is converted to a free energy
using the formula above and the free energy is output to a file called fes.dat

```plumed
x: DISTANCE ATOMS=1,2
hA1: HISTOGRAM ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1
ff: CONVERT_TO_FES ARG=hA1 TEMP=300
DUMPGRID ARG=ff FILE=fes.dat
```

In this example we ensure that all the reported values of the free energy are non-negative by finding the minimum in
the free energy and setting reporting the free energies relative to the estimate of the value of the free energy
at this  minimum.

```plumed
x: DISTANCE ATOMS=1,2
hA1: HISTOGRAM ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1
ff: CONVERT_TO_FES ARG=hA1 TEMP=300 MINTOZERO
DUMPGRID ARG=ff FILE=fes.dat
```

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class ConvertToFES : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit ConvertToFES(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(ConvertToFES,"CONVERT_TO_FES")

void ConvertToFES::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("compulsory","ARG","the histogram that you would like to convert into a free energy surface");
  keys.addDeprecatedKeyword("GRID","ARG");
  keys.add("optional","TEMP","the temperature at which you are operating");
  keys.addFlag("MINTOZERO",false,"set the minimum in the free energy to be equal to zero");
  keys.setValueDescription("grid","the free energy surface");
  keys.needsAction("FIND_GRID_MINIMUM");
  keys.needsAction("CUSTOM");
}

ConvertToFES::ConvertToFES(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  bool minzero=false;
  parseFlag("MINTOZERO",minzero);
  double simtemp=getkBT();
  if( simtemp==0 ) {
    error("TEMP not set - use keyword TEMP");
  }

  std::vector<std::string> argv;
  parseVector("GRID",argv);
  if( argv.size()==0 ) {
    parseVector("ARG",argv);
  }
  if( argv.size()!=1 ) {
    error("should only have one argument");
  }

  std::string str_temp;
  Tools::convert( simtemp, str_temp );
  std::string flab="";
  if( minzero ) {
    flab="_unz";
  }
  readInputLine( getShortcutLabel() + flab + ": CUSTOM ARG=" + argv[0] + " FUNC=-" + str_temp + "*log(x) PERIODIC=NO");
  if( minzero ) {
    readInputLine( getShortcutLabel() + "_min: FIND_GRID_MINIMUM ARG=" + getShortcutLabel() + "_unz" );
    readInputLine( getShortcutLabel() + ": CUSTOM ARG=" + getShortcutLabel() + "_unz," + getShortcutLabel() + "_min.optval FUNC=x-y PERIODIC=NO");
  }
}

}
}
