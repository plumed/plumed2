/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2023 The plumed team
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
#include "ActionWithGrid.h"

//+PLUMEDOC GRIDANALYSIS GET_VOLUME_ELEMENT
/*
Get the volume element from the input grid

This action is used in the shortcut action that is used within PLUMED to calculate integrals as you can see below

```plumed
c1: COORDINATIONNUMBER SPECIESA=1-10 SPECIESB=1-200 SWITCH={RATIONAL R_0=1.0}
hu: KDE ARG=c1 GRID_BIN=200 GRID_MIN=0 GRID_MAX=10 BANDWIDTH=0.1
h: CUSTOM ARG=hu FUNC=x/10 PERIODIC=NO
iv: INTEGRATE_GRID ARG=h PERIODIC=NO
```

As discused in the documentation for [INTEGRATE_GRID](INTEGRATE_GRID.md) we compute (approximate) integrals in PLUMED using:

$$
y = \frac{b-a}{n} \sum_{i=0}^{n} f\left( a + \frac{i(b-a)}{n} \right)
$$

This action gets the the value of $\frac{b-a}{n}$ so that we can multiply by this value by the sum in the expression
above, which is computed by a [SUM](SUM.md) action.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class GetGridVolumeElement :
  public ActionWithValue,
  public ActionWithArguments {
public:
  static void registerKeywords( Keywords& keys );
  explicit GetGridVolumeElement(const ActionOptions&ao);
  unsigned getNumberOfDerivatives() override ;
  void calculate() override ;
  void apply() override {}
};

PLUMED_REGISTER_ACTION(GetGridVolumeElement,"GET_VOLUME_ELEMENT")

void GetGridVolumeElement::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.addInputKeyword("compulsory","ARG","grid","the label for function on the grid that you need the volume element for");
  keys.setValueDescription("scalar","the volume element ");
  keys.remove("NUMERICAL_DERIVATIVES");
}

GetGridVolumeElement::GetGridVolumeElement(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao) {
  if( getNumberOfArguments()!=1 ) {
    error("should only be one argument to this action");
  }
  if( getPntrToArgument(0)->getRank()==0 || !getPntrToArgument(0)->hasDerivatives() ) {
    error("input to this action should be a grid");
  }
  ActionWithGrid* ag = dynamic_cast<ActionWithGrid*>( getPntrToArgument(0)->getPntrToAction() );
  if( !ag ) {
    error("could not find grid details from input value");
  }
  addValue();
  setNotPeriodic();
}

unsigned GetGridVolumeElement::getNumberOfDerivatives() {
  return 0;
}

void GetGridVolumeElement::calculate() {
  double volume;
  ActionWithGrid* ag = dynamic_cast<ActionWithGrid*>( getPntrToArgument(0)->getPntrToAction() );
  plumed_assert( ag );
  const GridCoordinatesObject& mygrid = ag->getGridCoordinatesObject();
  unsigned npoints = getPntrToArgument(0)->getNumberOfValues();
  if( mygrid.getGridType()=="flat" ) {
    std::vector<double> vv( mygrid.getGridSpacing() );
    volume=vv[0];
    for(unsigned i=1; i<vv.size(); ++i) {
      volume *=vv[i];
    }
  } else {
    volume = 4*pi / static_cast<double>( npoints );
  }
  getPntrToComponent(0)->set( volume );
}

}
}
