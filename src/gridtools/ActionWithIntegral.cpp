/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#include "ActionWithIntegral.h"

namespace PLMD {
namespace gridtools {

void ActionWithIntegral::registerKeywords( Keywords& keys ) {
  ActionWithInputGrid::registerKeywords( keys );
}

ActionWithIntegral::ActionWithIntegral(const ActionOptions&ao):
  Action(ao),
  ActionWithInputGrid(ao),
  volume(1)
{
  // Now create task list - number of tasks is equal to the number of grid points
  done_over_stream=false; std::vector<unsigned> shape;
  addValueWithDerivatives( shape ); setNotPeriodic();
}

void ActionWithIntegral::finishOutputSetup() {
  // Retrieve the volume of the grid (for integration)
  volume = gridobject.getCellVolume();
  // as we have to evaluate the function at each grid points
  for(unsigned i=0; i<getPntrToArgument(0)->getNumberOfValues( getLabel() ); ++i) addTaskToList(i);
  plumed_assert( arg_ends.size()==0 ); arg_ends.push_back(0); arg_ends.push_back(1);
  forcesToApply.resize( getPntrToArgument(0)->getNumberOfValues( getLabel() ) );
  // Turn on derivatives to resize derivatives in values now that we know the size
  if( !doNotCalculateDerivatives() ) turnOnDerivatives();
}

void ActionWithIntegral::apply() {
  if( doNotCalculateDerivatives() ) return;
  std::fill(forcesToApply.begin(),forcesToApply.end(),0); unsigned fstart=0;
  if( getForcesFromValues( forcesToApply ) ) setForcesOnArguments( 0, forcesToApply, fstart );
}

}
}
