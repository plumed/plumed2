/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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

void ManyRestraintsBase::registerKeywords( Keywords& keys ){
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

  // And turn on the derivatives (note problems here because of ActionWithValue)
  turnOnDerivatives(); needsDerivatives();

  // Now create the vessel
  std::string fake_input="LABEL=bias";
  addVessel( "SUM", fake_input, 0 ); 
  readVesselKeywords();
}

void ManyRestraintsBase::doJobsRequiredBeforeTaskList(){
  ActionWithVessel::doJobsRequiredBeforeTaskList();
  ActionWithValue::clearDerivatives();
}

void ManyRestraintsBase::applyChainRuleForDerivatives( const double& df ){
   // Value (this could be optimized more -- GAT)
   for(unsigned i=0;i<aves->getNumberOfDerivatives();++i){
       setElementDerivative( i, df*aves->getElementDerivative(i) );   
   }
   // And weights
   unsigned nder=aves->getNumberOfDerivatives();
   for(unsigned i=0;i<aves->getNumberOfDerivatives();++i){
       setElementDerivative( nder+i, aves->getElementDerivative(nder+i) );
   }
}

void ManyRestraintsBase::apply(){
  plumed_dbg_assert( getNumberOfComponents()==1 );
  getPntrToComponent(0)->addForce( -1.0*getStride() );
}

}
}

