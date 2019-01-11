/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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
#include "ManyRestraintsBase.h"
#include "vesselbase/Vessel.h"

namespace PLMD {
namespace manyrestraints {

void ManyRestraintsBase::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithVessel::registerKeywords( keys );
  ActionWithInputVessel::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  keys.add("hidden","STRIDE","the frequency with which the forces due to the bias should be calculated.  This can be used to correctly set up multistep algorithms");
  keys.remove("TOL");
  keys.addOutputComponent("bias","default","the instantaneous value of the bias potentials");
}

ManyRestraintsBase::ManyRestraintsBase(const ActionOptions& ao):
  Action(ao),
  ActionWithValue(ao),
  ActionPilot(ao),
  ActionWithVessel(ao),
  ActionWithInputVessel(ao)
{
  // Read in the vessel we are action on
  readArgument("bridge");
  aves=dynamic_cast<ActionWithVessel*>( getDependencies()[0] );

  plumed_assert( getDependencies().size()==1 && aves );
  log.printf("  adding restraints on variables calculated by %s action with label %s\n",
             aves->getName().c_str(),aves->getLabel().c_str());

  // Add a task list in order to avoid problems
  for(unsigned i=0; i<aves->getFullNumberOfTasks(); ++i) addTaskToList( aves->getTaskCode(i) );
  // And turn on the derivatives (note problems here because of ActionWithValue)
  turnOnDerivatives(); needsDerivatives();

  // Now create the vessel
  std::string fake_input="LABEL=bias";
  addVessel( "SUM", fake_input, 0 );
  readVesselKeywords();
}

void ManyRestraintsBase::doJobsRequiredBeforeTaskList() {
  ActionWithVessel::doJobsRequiredBeforeTaskList();
  ActionWithValue::clearDerivatives();
}

void ManyRestraintsBase::transformBridgedDerivatives( const unsigned& current, MultiValue& invals, MultiValue& outvals ) const {
  outvals.setValue( 0, invals.get(0) );

  // Get the potential
  double dval=0, val=calcPotential( invals.get(1), dval );

  outvals.setValue( 1, val );
  for(unsigned i=0; i<invals.getNumberActive(); ++i) {
    unsigned jder=invals.getActiveIndex(i);
    outvals.addDerivative( 1, jder, dval*invals.getDerivative( 1, jder ) );
  }

  // Now update the outvals derivatives lists
  outvals.emptyActiveMembers();
  for(unsigned j=0; j<invals.getNumberActive(); ++j) outvals.updateIndex( invals.getActiveIndex(j) );
  outvals.completeUpdate();
  return;
}

void ManyRestraintsBase::apply() {
  plumed_dbg_assert( getNumberOfComponents()==1 );
  getPntrToComponent(0)->addForce( -1.0*getStride() );
}

}
}

