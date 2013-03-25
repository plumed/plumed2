/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "ActionVolume.h"

namespace PLMD {
namespace multicolvar {

void ActionVolume::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithVessel::registerKeywords( keys );
  keys.use("MEAN"); keys.use("LESS_THAN"); keys.use("MORE_THAN");
  keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MIN");
  keys.add("compulsory","ARG","the label of the action that calculates the multicolvar we are interested in"); 
  keys.add("compulsory","SIGMA","the width of the function to be used for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian","the type of kernel function to be used");
  keys.addFlag("OUTSIDE",false,"calculate quantities for colvars that are on atoms outside the region of interest");
}

ActionVolume::ActionVolume(const ActionOptions&ao):
Action(ao),
ActionAtomistic(ao),
ActionWithValue(ao),
ActionWithVessel(ao)
{
  std::string mlab; parse("ARG",mlab);
  mycolv = plumed.getActionSet().selectWithLabel<multicolvar::MultiColvar*>(mlab);
  if(!mycolv) error("action labeled " + mlab + " does not exist or is not a multicolvar");
  addDependency(mycolv); mycolv->addBridgingVessel( this, myBridgeVessel );
  std::string functype=mycolv->getName();
  std::transform( functype.begin(), functype.end(), functype.begin(), tolower );
  log.printf("  calculating %s inside region of insterest\n",functype.c_str() ); 
 
  if( checkNumericalDerivatives() ){
      // If we use numerical derivatives we have to force the base
      // multicolvar to also use numerical derivatives
      ActionWithValue* vv=dynamic_cast<ActionWithValue*>( mycolv );
      plumed_assert( vv ); vv->useNumericalDerivatives();
  }

  parseFlag("OUTSIDE",not_in); parse("SIGMA",sigma); 
  bead.isPeriodic( 0, 1.0 ); 
  std::string kerneltype; parse("KERNEL",kerneltype); 
  bead.setKernelType( kerneltype );

  if( mycolv->isDensity() ){
     std::string input;
     addVessel( "SUM", input, -1 );  // -1 here means that this value will be named getLabel()
     resizeFunctions();
  } else {
     weightHasDerivatives=true;
     readVesselKeywords();
  }
}

void ActionVolume::doJobsRequiredBeforeTaskList(){
  ActionWithValue::clearDerivatives();
  retrieveAtoms(); setupRegion();
  ActionWithVessel::doJobsRequiredBeforeTaskList();
}

bool ActionVolume::performTask( const unsigned& j ){
  Vector catom_pos=mycolv->retrieveCentralAtomPos( derivativesOfFractionalCoordinates() );

  double weight; Vector wdf; 
  weight=calculateNumberInside( catom_pos, bead, wdf ); 
  if( not_in ){ weight = 1.0 - weight; wdf *= -1.; }  

  if( mycolv->isDensity() ){
     setElementValue( 0, weight ); setElementValue( 1, 1.0 ); 
     for(unsigned i=0;i<mycolv->getNAtoms();++i){
        unsigned nx=3*linkIndex( i, mycolv->colvar_atoms[mycolv->current], mycolv->all_atoms );
        addElementDerivative( nx+0, mycolv->getCentralAtomDerivative(i, 0, wdf ) );
        addElementDerivative( nx+1, mycolv->getCentralAtomDerivative(i, 1, wdf ) );
        addElementDerivative( nx+2, mycolv->getCentralAtomDerivative(i, 2, wdf ) );
     }
  } else {
     // Copy derivatives of the colvar and the value of the colvar
     double colv=mycolv->getElementValue(0); setElementValue( 0, colv );
     for(unsigned i=mycolv->getFirstDerivativeToMerge();i<mycolv->getNumberOfDerivatives();i=mycolv->getNextDerivativeToMerge(i)){
        addElementDerivative( i, mycolv->getElementDerivative(i) ); 
     }

     double ww=mycolv->getElementValue(1);
     setElementValue( 1, ww*weight ); 
     unsigned nder=getNumberOfDerivatives();

     // Add derivatives of weight if we have a weight
     if( mycolv->weightHasDerivatives ){
        for(unsigned i=mycolv->getFirstDerivativeToMerge();i<mycolv->getNumberOfDerivatives();i=mycolv->getNextDerivativeToMerge(i)){
           addElementDerivative( nder+i, weight*mycolv->getElementDerivative(nder+i) );
        } 
     }

     // Add derivatives of central atoms
     for(unsigned i=0;i<mycolv->atomsWithCatomDer.getNumberActive();++i){
         unsigned n=mycolv->atomsWithCatomDer[i];
         unsigned nx=nder+3*linkIndex( n, mycolv->colvar_atoms[mycolv->current], mycolv->all_atoms);
         addElementDerivative( nx+0, ww*mycolv->getCentralAtomDerivative(i, 0, wdf ) );
         addElementDerivative( nx+1, ww*mycolv->getCentralAtomDerivative(i, 1, wdf ) );
         addElementDerivative( nx+2, ww*mycolv->getCentralAtomDerivative(i, 2, wdf ) );
     }
  }

  // Only continue if the weight is greater than tolerance
  return ( weight>=getTolerance() );
}

void ActionVolume::calculateNumericalDerivatives(){
  myBridgeVessel->completeNumericalDerivatives();
}

bool ActionVolume::isPeriodic(){
  return mycolv->isPeriodic();
}

void ActionVolume::deactivate_task(){
  plumed_merror("This should never be called");
}

}
}
