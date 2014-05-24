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
#include "VolumeGradientBase.h"

namespace PLMD {
namespace multicolvar {

void VolumeGradientBase::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithVessel::registerKeywords( keys );
  ActionWithInputVessel::registerKeywords( keys );
  keys.setComponentsIntroduction("This Action can be used to calculate the following quantities by employing the keywords listed below. "
                                 "You must select which from amongst these quantities you wish to calculate - this command cannot be run "
                                 "unless one of the quantities below is being calculated."
                                 "These quantities can then be referenced elsewhere in the input file by using this Action's label "
                                 "followed by a dot and the name of the quantity. Some amongst them can be calculated multiple times "
                                 "with different parameters.  In this case the quantities calculated can be referenced elsewhere in the "
                                 "input by using the name of the quantity followed by a numerical identifier "
                                 "e.g. <em>label</em>.less_than-1, <em>label</em>.less_than-2 etc.");
  keys.use("NL_TOL");
  keys.add("hidden","NL_STRIDE","the frequency with which the neighbor list should be updated. Between neighbour list update steps all quantities "
                                "that contributed less than TOL at the previous neighbor list update step are ignored.");
}

VolumeGradientBase::VolumeGradientBase(const ActionOptions&ao):
Action(ao),
ActionAtomistic(ao),
ActionWithValue(ao),
ActionWithVessel(ao),
ActionWithInputVessel(ao),
updateFreq(0)
{
  readArgument("bridge");
  mycolv = dynamic_cast<MultiColvarBase*>( getDependencies()[0] ); 
  plumed_assert( getDependencies().size()==1 );
  if(!mycolv) error("action labeled " + mycolv->getLabel() + " is not a multicolvar"); 

  // Neighbor list readin
  parse("NL_STRIDE",updateFreq);
  if(updateFreq>0){
     if( !mycolv->isDensity() && updateFreq%mycolv->updateFreq!=0 ){ 
        error("update frequency must be a multiple of update frequency for base multicolvar");  
     }
     firsttime=true;
     log.printf("  Updating contributors every %d steps.\n",updateFreq);
  } else {
     firsttime=false; contributorsAreUnlocked=true; // This will lock during first prepare step
     log.printf("  Updating contributors every step.\n");
  }
  weightHasDerivatives=true;
}

void VolumeGradientBase::turnOnDerivatives(){
  ActionWithValue::turnOnDerivatives();
  needsDerivatives();
} 

void VolumeGradientBase::requestAtoms( const std::vector<AtomNumber>& atoms ){
  ActionAtomistic::requestAtoms(atoms); bridgeVariable=3*atoms.size();
  addDependency( mycolv ); 
  tmpforces.resize( 3*atoms.size()+9 );
}

void VolumeGradientBase::doJobsRequiredBeforeTaskList(){
  ActionWithValue::clearDerivatives();
  retrieveAtoms(); setupRegions();
  ActionWithVessel::doJobsRequiredBeforeTaskList();
}

void VolumeGradientBase::prepare(){
  bool updatetime=false;
  if( contributorsAreUnlocked ){ updatetime=true; lockContributors(); }
  if( updateFreq>0 ){
      if( firsttime || getStep()%updateFreq==0 ){
         if( !mycolv->isDensity() ){
             mycolv->unlockContributors(); 
         } else {
             plumed_massert( mycolv->contributorsAreUnlocked, "contributors are not unlocked in base multicolvar" );
         }
         unlockContributors(); 
         firsttime=false;
      }
  }
}

void VolumeGradientBase::finishTaskListUpdate(){
  activeAtoms.clear();
  for(unsigned i=0;i<mycolv->getNumberOfAtoms();++i) activeAtoms.addIndexToList(i);
  activeAtoms.deactivateAll();
}

void VolumeGradientBase::performTask(){
  activeAtoms.deactivateAll(); // Currently no atoms are active so deactivate them all
//  Vector catom_pos=mycolv->retrieveCentralAtomPos();

//  double weight; Vector wdf; 
//  weight=calculateNumberInside( catom_pos, bead, wdf ); 
//  if( not_in ){ weight = 1.0 - weight; wdf *= -1.; }  

  if( mycolv->isDensity() ){ 
     setElementValue( 0, 1.0 );
     // unsigned nder=getNumberOfDerivatives();
     // setElementValue( 1, weight ); setElementValue( 0, 1.0 );
     // for(unsigned i=0;i<mycolv->atomsWithCatomDer.getNumberActive();++i){
     //    unsigned n=mycolv->atomsWithCatomDer[i], nx=nder + 3*n;
     //    activeAtoms.activate(n);
     //    addElementDerivative( nx+0, mycolv->getCentralAtomDerivative(n, 0, wdf ) );
     //    addElementDerivative( nx+1, mycolv->getCentralAtomDerivative(n, 1, wdf ) );
     //    addElementDerivative( nx+2, mycolv->getCentralAtomDerivative(n, 2, wdf ) );
     // }
  } else {
     // Copy derivatives of the colvar and the value of the colvar
     mycolv->copyElementsToBridgedColvar( this );
     //  // Add derivatives of central atoms
     //  double ww=mycolv->getElementValue(1);
     //  setElementValue( 1, ww*weight );
     //  unsigned nder=getNumberOfDerivatives(); 
     //  for(unsigned i=0;i<mycolv->atomsWithCatomDer.getNumberActive();++i){
     //      unsigned n=mycolv->atomsWithCatomDer[i], nx=nder+3*n;
     //      activeAtoms.activate(n); 
     //      addElementDerivative( nx+0, ww*mycolv->getCentralAtomDerivative(n, 0, wdf ) );
     //      addElementDerivative( nx+1, ww*mycolv->getCentralAtomDerivative(n, 1, wdf ) );
     //      addElementDerivative( nx+2, ww*mycolv->getCentralAtomDerivative(n, 2, wdf ) );
     //  }
  }
  calculateAllVolumes(); 
  activeAtoms.updateActiveMembers();
}

void VolumeGradientBase::setNumberInVolume( const unsigned& ivol, const double& weight, const Vector& wdf ){
  if( !mycolv->weightHasDerivatives ){
      unsigned nstart=ivol*getNumberOfDerivatives();
      setElementValue( ivol, weight );
      for(unsigned i=0;i<mycolv->atomsWithCatomDer.getNumberActive();++i){
         unsigned n=mycolv->atomsWithCatomDer[i], nx=nstart + 3*n;
         activeAtoms.activate(n);
         addElementDerivative( nx+0, mycolv->getCentralAtomDerivative(n, 0, wdf ) );
         addElementDerivative( nx+1, mycolv->getCentralAtomDerivative(n, 1, wdf ) );
         addElementDerivative( nx+2, mycolv->getCentralAtomDerivative(n, 2, wdf ) );
      } 
  } else {
      unsigned nstart=ivol*getNumberOfDerivatives();
      double ww=mycolv->getElementValue(1); setElementValue( ivol, ww*weight );
      for(unsigned i=0;i<mycolv->atomsWithCatomDer.getNumberActive();++i){
          unsigned n=mycolv->atomsWithCatomDer[i], nx=nstart + 3*n;
          activeAtoms.activate(n); 
          addElementDerivative( nx+0, ww*mycolv->getCentralAtomDerivative(n, 0, wdf ) );
          addElementDerivative( nx+1, ww*mycolv->getCentralAtomDerivative(n, 1, wdf ) );
          addElementDerivative( nx+2, ww*mycolv->getCentralAtomDerivative(n, 2, wdf ) );
     }
     if( mycolv->weightHasDerivatives ){
         unsigned nder=mycolv->getNumberOfDerivatives(); 
         for(unsigned i=0;i<mycolv->atoms_with_derivatives.getNumberActive();++i){
            unsigned n=mycolv->atoms_with_derivatives[i], nx=nder + 3*n, ny=nstart + 3*n;
            activeAtoms.activate(n);
            addElementDerivative( ny+0, weight*mycolv->getElementDerivative(nx+0) );
            addElementDerivative( ny+1, weight*mycolv->getElementDerivative(nx+1) );
            addElementDerivative( ny+2, weight*mycolv->getElementDerivative(nx+2) );
         }
         unsigned nwvir=mycolv->getNumberOfDerivatives()-9, nwvir2=nstart + 3*mycolv->getNumberOfAtoms();
         for(unsigned i=0;i<9;++i){
            addElementDerivative( nwvir2, weight*mycolv->getElementDerivative(nwvir) ); nwvir++; nwvir2++;
         }
     }
  }
}

void VolumeGradientBase::mergeDerivatives( const unsigned& ider, const double& df ){
  unsigned vstart=getNumberOfDerivatives()*ider;
  // Merge atom derivatives
  for(unsigned i=0;i<activeAtoms.getNumberActive();++i){
     unsigned iatom=3*activeAtoms[i];
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

void VolumeGradientBase::clearDerivativesAfterTask( const unsigned& ider ){
  unsigned vstart=getNumberOfDerivatives()*ider;
  if( derivativesAreRequired() ){
     // Clear atom derivatives
     for(unsigned i=0;i<activeAtoms.getNumberActive();++i){
        unsigned iatom=vstart+3*activeAtoms[i];
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

void VolumeGradientBase::calculateNumericalDerivatives( ActionWithValue* a ){
  ActionWithInputVessel::calculateNumericalDerivatives(a);
}

bool VolumeGradientBase::isPeriodic(){
  return mycolv->isPeriodic();
}

void VolumeGradientBase::deactivate_task(){
  plumed_merror("This should never be called");
}

void VolumeGradientBase::addBridgeForces( const std::vector<double>& bb ){ 
  plumed_dbg_assert( bb.size()==tmpforces.size()-9 );
  // Forces on local atoms
  for(unsigned i=0;i<bb.size();++i) tmpforces[i]=bb[i];
  // Virial contribution is zero
  for(unsigned i=bb.size();i<bb.size()+9;++i) tmpforces[i]=0.0;
  setForcesOnAtoms( tmpforces, 0 );
}

}
}
