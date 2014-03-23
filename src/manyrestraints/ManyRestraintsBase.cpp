/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
  ActionAtomistic::registerKeywords( keys );
  ActionWithVessel::registerKeywords( keys );
  ActionPilot::registerKeywords( keys );
  keys.add("hidden","STRIDE","the frequency with which the forces due to the bias should be calculated.  This can be used to correctly set up multistep algorithms");
  keys.remove("TOL");
}

ManyRestraintsBase::ManyRestraintsBase(const ActionOptions& ao):
Action(ao),
ActionAtomistic(ao),
ActionWithValue(ao),
ActionPilot(ao),
ActionWithVessel(ao)
{
}

void ManyRestraintsBase::turnOnDerivatives(){
  error("restraints cannot be used as collective variables");
}

void ManyRestraintsBase::createRestraints( const unsigned& nrestraints ){
  std::string fake_input; 
  for(unsigned i=0;i<nrestraints;++i) addTaskToList(i);
  addVessel( "SUM", fake_input, 0, "bias" );
  readVesselKeywords();
  forcesToApply.resize( getNumberOfDerivatives() );
}

void ManyRestraintsBase::apply(){
  plumed_dbg_assert( getNumberOfComponents()==1 );
  getPntrToComponent(0)->addForce(-1.0);
  bool wasforced=getForcesFromVessels( forcesToApply );
  plumed_assert( wasforced ); setForcesOnAtoms( forcesToApply );
}

}
}

