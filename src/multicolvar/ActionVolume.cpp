/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
  keys.use("BETWEEN"); keys.use("HISTOGRAM"); 
  keys.add("compulsory","ARG","the label of the action that calculates the multicolvar we are interested in"); 
  keys.add("compulsory","SIGMA","the width of the function to be used for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian","the type of kernel function to be used");
  keys.addFlag("OUTSIDE",false,"calculate quantities for colvars that are on atoms outside the region of interest");
  keys.use("NL_TOL");
  keys.add("hidden","NL_STRIDE","the frequency with which the neighbor list should be updated. Between neighbour list update steps all quantities "
                                "that contributed less than TOL at the previous neighbor list update step are ignored.");
}

ActionVolume::ActionVolume(const ActionOptions&ao):
Action(ao),
ActionAtomistic(ao),
ActionWithValue(ao),
ActionWithVessel(ao),
updateFreq(0),
lastUpdate(0)
{
  std::string mlab; parse("ARG",mlab);
  mycolv = plumed.getActionSet().selectWithLabel<multicolvar::MultiColvarBase*>(mlab);
  if(!mycolv) error("action labeled " + mlab + " does not exist or is not a multicolvar");
  std::string functype=mycolv->getName();
  std::transform( functype.begin(), functype.end(), functype.begin(), tolower );
  log.printf("  calculating %s inside region of insterest\n",functype.c_str() ); 
 
  if( checkNumericalDerivatives() ){
      // If we use numerical derivatives we have to force the base
      // multicolvar to also use numerical derivatives
      ActionWithValue* vv=dynamic_cast<ActionWithValue*>( mycolv );
      plumed_assert( vv ); vv->useNumericalDerivatives();
  }

  // Neighbor list readin
  parse("NL_STRIDE",updateFreq);
  if(updateFreq>0){
     if( !mycolv->isDensity() && updateFreq%mycolv->updateFreq!=0 ){ 
        error("update frequency must be a multiple of update frequency for base multicolvar");  
     }
     log.printf("  Updating contributors every %d steps.\n",updateFreq);
  } else {
     log.printf("  Updating contributors every step.\n");
  }

  parseFlag("OUTSIDE",not_in); parse("SIGMA",sigma); 
  bead.isPeriodic( 0, 1.0 ); 
  std::string kerneltype; parse("KERNEL",kerneltype); 
  bead.setKernelType( kerneltype );
  weightHasDerivatives=true;
  
  if( mycolv->isDensity() ){
     std::string input;
     addVessel( "SUM", input, -1 );  // -1 here means that this value will be named getLabel()
     resizeFunctions();
  } else {
     readVesselKeywords();
  }

  // Now set up the bridging vessel (has to be done this way for internal arrays to be resized properly)
  addDependency(mycolv); myBridgeVessel = mycolv->addBridgingVessel( this );
}

void ActionVolume::doJobsRequiredBeforeTaskList(){
  ActionWithValue::clearDerivatives();
  retrieveAtoms(); setupRegion();
  ActionWithVessel::doJobsRequiredBeforeTaskList();
}

void ActionVolume::prepare(){
  if( contributorsAreUnlocked ) lockContributors();
  if( updateFreq>0 && (getStep()-lastUpdate)>=updateFreq ){
      if( !mycolv->isDensity() ){
          mycolv->taskList.activateAll();
          for(unsigned i=0;i<mycolv->taskList.getNumberActive();++i) mycolv->colvar_atoms[i].activateAll();
          mycolv->unlockContributors(); mycolv->resizeDynamicArrays();
          plumed_dbg_assert( mycolv->getNumberOfVessels()==0 );
      } else {
          plumed_massert( mycolv->contributorsAreUnlocked, "contributors are not unlocked in base multicolvar" );
      }
      unlockContributors(); 
  }
}

void ActionVolume::performTask( const unsigned& j ){
  Vector catom_pos=mycolv->retrieveCentralAtomPos( derivativesOfFractionalCoordinates() );

  double weight; Vector wdf; 
  weight=calculateNumberInside( catom_pos, bead, wdf ); 
  if( not_in ){ weight = 1.0 - weight; wdf *= -1.; }  

  if( mycolv->isDensity() ){
     unsigned nder=getNumberOfDerivatives(); 
     setElementValue( 1, weight ); setElementValue( 0, 1.0 ); 
     for(unsigned i=0;i<mycolv->atomsWithCatomDer.getNumberActive();++i){
        unsigned n=mycolv->atomsWithCatomDer[i], nx=nder + 3*n;
        addElementDerivative( nx+0, mycolv->getCentralAtomDerivative(n, 0, wdf ) );
        addElementDerivative( nx+1, mycolv->getCentralAtomDerivative(n, 1, wdf ) );
        addElementDerivative( nx+2, mycolv->getCentralAtomDerivative(n, 2, wdf ) );
     }
     // Virial contribution
     if( !derivativesOfFractionalCoordinates() ){
        Tensor vir( catom_pos, wdf ); nder+=3*mycolv->getNumberOfAtoms();
        addElementDerivative( nder+0, -vir(0,0) ); 
        addElementDerivative( nder+1, -vir(0,1) );
        addElementDerivative( nder+2, -vir(0,2) );
        addElementDerivative( nder+3, -vir(1,0) );
        addElementDerivative( nder+4, -vir(1,1) );
        addElementDerivative( nder+5, -vir(1,2) );
        addElementDerivative( nder+6, -vir(2,0) );
        addElementDerivative( nder+7, -vir(2,1) );
        addElementDerivative( nder+8, -vir(2,2) );
     }
  } else {
     // Copy derivatives of the colvar and the value of the colvar
     double colv=mycolv->getElementValue(0); setElementValue( 0, colv );
     for(unsigned i=0;i<mycolv->atoms_with_derivatives.getNumberActive();++i){
        unsigned n=3*mycolv->atoms_with_derivatives(i);
        addElementDerivative( n, mycolv->getElementDerivative(n) ); n++;
        addElementDerivative( n, mycolv->getElementDerivative(n) ); n++;
        addElementDerivative( n, mycolv->getElementDerivative(n) ); 
     }
     unsigned nvir=3*mycolv->getNumberOfAtoms();
     for(unsigned i=0;i<9;++i){ 
       addElementDerivative( nvir, mycolv->getElementDerivative(nvir) ); nvir++;
     }

     double ww=mycolv->getElementValue(1);
     setElementValue( 1, ww*weight ); 
     unsigned nder=getNumberOfDerivatives();

     // Add derivatives of weight if we have a weight
     if( mycolv->weightHasDerivatives ){
        for(unsigned i=0;i<mycolv->atoms_with_derivatives.getNumberActive();++i){
           unsigned n=3*mycolv->atoms_with_derivatives(i);
           addElementDerivative( nder+n, weight*mycolv->getElementDerivative(nder+n) ); n++;
           addElementDerivative( nder+n, weight*mycolv->getElementDerivative(nder+n) ); n++;
           addElementDerivative( nder+n, weight*mycolv->getElementDerivative(nder+n) );
        } 
        unsigned nwvir=3*mycolv->getNumberOfAtoms();
        for(unsigned i=0;i<9;++i){
           addElementDerivative( nwvir, mycolv->getElementDerivative(nwvir) ); nwvir++; 
        }
     }

     // Add derivatives of central atoms
     for(unsigned i=0;i<mycolv->atomsWithCatomDer.getNumberActive();++i){
         unsigned n=mycolv->atomsWithCatomDer[i], nx=nder+3*n; 
         addElementDerivative( nx+0, ww*mycolv->getCentralAtomDerivative(n, 0, wdf ) );
         addElementDerivative( nx+1, ww*mycolv->getCentralAtomDerivative(n, 1, wdf ) );
         addElementDerivative( nx+2, ww*mycolv->getCentralAtomDerivative(n, 2, wdf ) );
     }
     // Virial contribution for central atoms
     if( !derivativesOfFractionalCoordinates() ){
        Tensor vir( catom_pos, wdf ); nder+=3*mycolv->getNumberOfAtoms();
        addElementDerivative( nder+0, -ww*vir(0,0) );
        addElementDerivative( nder+1, -ww*vir(0,1) );
        addElementDerivative( nder+2, -ww*vir(0,2) );
        addElementDerivative( nder+3, -ww*vir(1,0) );
        addElementDerivative( nder+4, -ww*vir(1,1) );
        addElementDerivative( nder+5, -ww*vir(1,2) );
        addElementDerivative( nder+6, -ww*vir(2,0) );
        addElementDerivative( nder+7, -ww*vir(2,1) );
        addElementDerivative( nder+8, -ww*vir(2,2) );
     }
  }
}

void ActionVolume::calculateNumericalDerivatives( ActionWithValue* a ){
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
