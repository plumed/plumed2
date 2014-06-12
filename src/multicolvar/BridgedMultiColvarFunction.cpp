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
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "BridgedMultiColvarFunction.h"

namespace PLMD {
namespace multicolvar {

void BridgedMultiColvarFunction::registerKeywords( Keywords& keys ){
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
  myBridgeVessel = mycolv->addBridgingVessel( this ); addDependency(mycolv);
  // Checks for neighbor list
  if( mycolv->updateFreq>0 && updateFreq>0 ){
      if( updateFreq%mycolv->updateFreq!=0 ) error("update frequency must be a multiple of the update frequency in the base colvar");
  }
  weightHasDerivatives=true;
}

void BridgedMultiColvarFunction::finishTaskListUpdate(){
  plumed_assert( all_atoms.fullSize()==0 );
  resizeLocalArrays();
}

void BridgedMultiColvarFunction::mergeDerivatives( const unsigned& ider, const double& df ){
  unsigned vstart=getNumberOfDerivatives()*ider;
  // Merge atom derivatives
  for(unsigned i=0;i<atoms_with_derivatives.getNumberActive();++i){
     unsigned iatom=3*atoms_with_derivatives[i];
     accumulateDerivative( iatom, df*getElementDerivative(vstart+iatom) ); iatom++;
     accumulateDerivative( iatom, df*getElementDerivative(vstart+iatom) ); iatom++;
     accumulateDerivative( iatom, df*getElementDerivative(vstart+iatom) );
  }
  // Merge virial derivatives
  unsigned nvir=3*mycolv->getNumberOfAtoms();
  for(unsigned j=0;j<9;++j){
     accumulateDerivative( nvir, df*getElementDerivative(vstart+nvir) ); nvir++;
  }
  // Merge local atom derivatives
  for(unsigned j=0;j<getNumberOfAtoms();++j){
     accumulateDerivative( nvir, df*getElementDerivative(vstart+nvir) ); nvir++;
     accumulateDerivative( nvir, df*getElementDerivative(vstart+nvir) ); nvir++;
     accumulateDerivative( nvir, df*getElementDerivative(vstart+nvir) ); nvir++;
  }
  plumed_dbg_assert( nvir==getNumberOfDerivatives() );
}

void BridgedMultiColvarFunction::clearDerivativesAfterTask( const unsigned& ider ){
  unsigned vstart=getNumberOfDerivatives()*ider;
  if( derivativesAreRequired() ){
     // Clear atom derivatives
     for(unsigned i=0;i<atoms_with_derivatives.getNumberActive();++i){
        unsigned iatom=vstart+3*atoms_with_derivatives[i];
        setElementDerivative( iatom, 0.0 ); iatom++;
        setElementDerivative( iatom, 0.0 ); iatom++;
        setElementDerivative( iatom, 0.0 );
     }
     // Clear virial contribution
     unsigned nvir=vstart+3*mycolv->getNumberOfAtoms();
     for(unsigned j=0;j<9;++j){
        setElementDerivative( nvir, 0.0 ); nvir++;
     }
     // Clear derivatives of local atoms
     for(unsigned j=0;j<getNumberOfAtoms();++j){
        setElementDerivative( nvir, 0.0 ); nvir++;
        setElementDerivative( nvir, 0.0 ); nvir++;
        setElementDerivative( nvir, 0.0 ); nvir++;
     }
     plumed_dbg_assert( (nvir-vstart)==getNumberOfDerivatives() );
  }
  // Clear values
  thisval_wasset[ider]=false; setElementValue( ider, 0.0 ); thisval_wasset[ider]=false;
}

void BridgedMultiColvarFunction::calculateNumericalDerivatives( ActionWithValue* a ){
  if(!a){
    a=dynamic_cast<ActionWithValue*>(this);
    plumed_massert(a,"cannot compute numerical derivatives for an action without values");
  }
  if( myBridgeVessel ){
     myBridgeVessel->completeNumericalDerivatives();
  } else {
     error("numerical derivatives are not implemented");
  }
}

bool BridgedMultiColvarFunction::isPeriodic(){
  return mycolv->isPeriodic();
}

void BridgedMultiColvarFunction::deactivate_task(){
  plumed_merror("This should never be called");
}

}
}
