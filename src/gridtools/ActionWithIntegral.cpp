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
  keys.remove("KERNEL"); keys.remove("BANDWIDTH");
  keys.remove("CLEAR"); keys.add("compulsory","CLEAR","1","the frequency with which to clear all the accumulated data.");
}

ActionWithIntegral::ActionWithIntegral(const ActionOptions&ao):
  Action(ao),
  ActionWithInputGrid(ao)
{
  plumed_assert( ingrid->getNumberOfComponents()==1 );
  // Retrieve the volume of the grid (for integration)
  volume = ingrid->getCellVolume();
  // Create something that is going to calculate the sum of all the values
  // at the various grid points - this is going to be the integral
  std::string fake_input; addVessel( "SUM", fake_input, -1 ); readVesselKeywords();
  // Now create task list - number of tasks is equal to the number of grid points
  // as we have to evaluate the function at each grid points
  for(unsigned i=0; i<ingrid->getNumberOfPoints(); ++i) addTaskToList(i);
  // And activate all tasks
  deactivateAllTasks();
  for(unsigned i=0; i<ingrid->getNumberOfPoints(); ++i) taskFlags[i]=1;
  lockContributors();
}

void ActionWithIntegral::turnOnDerivatives() {
  ActionWithGrid::turnOnDerivatives();
  forcesToApply.resize( ingrid->getNumberOfPoints() );
}

void ActionWithIntegral::apply() {
  if( getForcesFromVessels( forcesToApply ) ) ingrid->setForce( forcesToApply );
}

}
}
