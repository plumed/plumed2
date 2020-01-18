/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#include "ActionWithIntegral.h"

//+PLUMEDOC GRIDANALYSIS INTEGRATE_GRID
/*
Calculate the total integral of the function on the input grid

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class IntegrateGrid : public ActionWithIntegral {
public:
  static void registerKeywords( Keywords& keys );
  explicit IntegrateGrid(const ActionOptions&ao);
  void compute( const unsigned& current, MultiValue& myvals ) const override;
};

PLUMED_REGISTER_ACTION(IntegrateGrid,"INTEGRATE_GRID")

void IntegrateGrid::registerKeywords( Keywords& keys ) {
  ActionWithIntegral::registerKeywords( keys );
}

IntegrateGrid::IntegrateGrid(const ActionOptions&ao):
  Action(ao),
  ActionWithIntegral(ao)
{
}

void IntegrateGrid::compute( const unsigned& current, MultiValue& myvals ) const {
  myvals.setValue( 0, 1.0 ); myvals.setValue( 1, getVolume()*getFunctionValue( current ) );
  if( !doNotCalculateDerivatives() ) myvals.addDerivative( 1, current, getVolume() );
}

}
}
