/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2019 The plumed team
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
#include "BridgedMultiColvarFunction.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "CatomPack.h"

namespace PLMD {
namespace multicolvar {

void BridgedMultiColvarFunction::registerKeywords( Keywords& keys ) {
  MultiColvarBase::registerKeywords( keys );
  keys.add("compulsory","DATA","The multicolvar that calculates the set of base quantities that we are interested in");
}

BridgedMultiColvarFunction::BridgedMultiColvarFunction(const ActionOptions&ao):
  Action(ao),
  MultiColvarBase(ao)
{
  std::string mlab; parse("DATA",mlab);
  mycolv = plumed.getActionSet().selectWithLabel<MultiColvarBase*>(mlab);
  if(!mycolv) error("action labeled " + mlab + " does not exist or is not a multicolvar");

  // When using numerical derivatives here we must use numerical derivatives
  // in base multicolvar
  if( checkNumericalDerivatives() ) mycolv->useNumericalDerivatives();

  myBridgeVessel = mycolv->addBridgingVessel( this ); addDependency(mycolv);
  weightHasDerivatives=true; usespecies=mycolv->usespecies;
  // Number of tasks is the same as the number in the underlying MultiColvar
  for(unsigned i=0; i<mycolv->getFullNumberOfTasks(); ++i) addTaskToList( mycolv->getTaskCode(i) );
}

void BridgedMultiColvarFunction::turnOnDerivatives() {
  BridgedMultiColvarFunction* check = dynamic_cast<BridgedMultiColvarFunction*>( mycolv );
  if( check ) {
    if( check->getNumberOfAtoms()>0 ) error("cannot calculate required derivatives of this quantity");
  }
  MultiColvarBase::turnOnDerivatives();
}

void BridgedMultiColvarFunction::transformBridgedDerivatives( const unsigned& current, MultiValue& invals, MultiValue& outvals ) const {
  completeTask( current, invals, outvals );

  // Now update the outvals derivatives lists
  if( derivativesAreRequired() ) {
    outvals.emptyActiveMembers();
    if( mycolv->isDensity() ) {
      for(unsigned j=0; j<3; ++j) outvals.putIndexInActiveArray( 3*current+j );
      for(unsigned j=invals.getNumberOfDerivatives()-9; j<invals.getNumberOfDerivatives(); ++j) outvals.putIndexInActiveArray(j);
    } else {
      for(unsigned j=0; j<invals.getNumberActive(); ++j) outvals.putIndexInActiveArray( invals.getActiveIndex(j) );
    }
    for(unsigned j=invals.getNumberOfDerivatives(); j<outvals.getNumberOfDerivatives(); ++j) outvals.putIndexInActiveArray( j );
    outvals.completeUpdate();
  }
}

void BridgedMultiColvarFunction::performTask( const unsigned& taskIndex, const unsigned& current, MultiValue& myvals ) const {
  // This allows us to speed up the code as we don't need to reallocate memory on every call of perform task
  MultiValue& invals=myBridgeVessel->getTemporyMultiValue();
  if( invals.getNumberOfValues()!=mycolv->getNumberOfQuantities() ||
      invals.getNumberOfDerivatives()!=mycolv->getNumberOfDerivatives() ) {
    invals.resize( mycolv->getNumberOfQuantities(), mycolv->getNumberOfDerivatives() );
  }
  invals.clearAll(); mycolv->performTask( taskIndex, current, invals );
  transformBridgedDerivatives( taskIndex, invals, myvals );
}

void BridgedMultiColvarFunction::calculateNumericalDerivatives( ActionWithValue* a ) {
  if(!a) {
    a=dynamic_cast<ActionWithValue*>(this);
    plumed_massert(a,"cannot compute numerical derivatives for an action without values");
  }
  if( myBridgeVessel ) {
    myBridgeVessel->completeNumericalDerivatives();
  } else {
    error("numerical derivatives are not implemented");
  }
}

void BridgedMultiColvarFunction::applyBridgeForces( const std::vector<double>& bb ) {
  if( getNumberOfAtoms()==0 ) return ;

  std::vector<Vector>& f( modifyForces() );
  for(unsigned i=0; i<getNumberOfAtoms(); ++i) {
    f[i][0]+=bb[3*i+0]; f[i][1]+=bb[3*i+1]; f[i][2]+=bb[3*i+2];
  }
  applyForces();
}

bool BridgedMultiColvarFunction::isPeriodic() {
  return mycolv->isPeriodic();
}

void BridgedMultiColvarFunction::deactivate_task( const unsigned& taskno ) {
  plumed_merror("This should never be called");
}

void BridgedMultiColvarFunction::getCentralAtomPack( const unsigned& basn, const unsigned& curr, CatomPack& mypack ) {
  return mycolv->getCentralAtomPack( basn, curr, mypack );
}

}
}
