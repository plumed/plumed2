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
#include "MultiColvarBase.h"
#include "MultiColvarFunction.h"
#include "vesselbase/Vessel.h"
#include "tools/Pbc.h"
#include <vector>
#include <string>

using namespace std;

namespace PLMD{
namespace multicolvar{

void MultiColvarBase::registerKeywords( Keywords& keys ){
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  ActionWithVessel::registerKeywords( keys );
  keys.use("NL_TOL");
  keys.add("hidden","NL_STRIDE","the frequency with which the neighbor list should be updated. Between neighbour list update steps all quantities "
                                "that contributed less than TOL at the previous neighbor list update step are ignored.");
  keys.setComponentsIntroduction("This Action can be used to calculate the following quantities by employing the keywords listed below. "
                                 "You must select which from amongst these quantities you wish to calculate - this command cannot be run "
                                 "unless one of the quantities below is being calculated."
                                 "These quantities can then be referenced elsewhere in the input file by using this Action's label "
                                 "followed by a dot and the name of the quantity. Some amongst them can be calculated multiple times "
                                 "with different parameters.  In this case the quantities calculated can be referenced elsewhere in the "
                                 "input by using the name of the quantity followed by a numerical identifier "
                                 "e.g. <em>label</em>.less_than-1, <em>label</em>.less_than-2 etc.");
} 

MultiColvarBase::MultiColvarBase(const ActionOptions&ao):
Action(ao),
ActionAtomistic(ao),
ActionWithValue(ao),
ActionWithVessel(ao),
usepbc(false),
updateFreq(0),
mycatoms(NULL),        // This will be destroyed by ActionWithVesel
myvalues(NULL),        // This will be destroyed by ActionWithVesel 
usespecies(false)
{
  if( keywords.exists("NOPBC") ){ 
    bool nopbc=!usepbc; parseFlag("NOPBC",nopbc);
    usepbc=!nopbc;
  } 
  if( keywords.exists("SPECIES") ) usespecies=true;
  if( keywords.exists("NL_STRIDE") ) parse("NL_STRIDE",updateFreq);
  if(updateFreq>0){
     firsttime=true;
     log.printf("  Updating contributors every %d steps.\n",updateFreq);
  } else {
     firsttime=false; contributorsAreUnlocked=true;   // This will lock during first prepare step methinks
     log.printf("  Updating contributors every step.\n");
  }
}

void MultiColvarBase::copyAtomListToFunction( MultiColvarBase* myfunction ){
  for(unsigned i=0;i<all_atoms.fullSize();++i) myfunction->all_atoms.addIndexToList( all_atoms(i) );
}

void MultiColvarBase::copyActiveAtomsToFunction( MultiColvarBase* myfunction, const unsigned& start ){
  for(unsigned i=0;i<all_atoms.getNumberActive();++i){
      unsigned iatom=all_atoms.linkIndex( i );
      myfunction->all_atoms.activate( start + iatom );
  }
}

void MultiColvarBase::setupMultiColvarBase(){
  // Setup decoder array
  if( !usespecies && ablocks.size()<4 ){
     decoder.resize( ablocks.size() ); unsigned code=1;
     for(unsigned i=0;i<ablocks.size();++i){ decoder[ablocks.size()-1-i]=code; code *= nblock; } 
  } else if( ablocks.size()==1 ) {
     // Setup coordination sphere
     csphere_atoms.resize( getFullNumberOfTasks() ); unsigned nflags=0;
     for(unsigned i=0;i<getFullNumberOfTasks();++i){
        for(unsigned j=0;j<ablocks[0].size();++j){
           if( !same_index( getActiveTask(i), ablocks[0][j] ) ){
               csphere_atoms[i].addIndexToList( ablocks[0][j] ); nflags++;
           }
        }
        csphere_atoms[i].activateAll();
     } 
     csphere_flags.resize( nflags, 0 );
  } 
  // Do an initial task list update
  finishTaskListUpdate();
  // Setup underlying ActionWithVessel
  readVesselKeywords();
}

void MultiColvarBase::turnOnDerivatives(){
  ActionWithValue::turnOnDerivatives();
  needsDerivatives();
} 

void MultiColvarBase::prepare(){
  if( contributorsAreUnlocked ) lockContributors();
  if( updateFreq>0 ){
     if( firsttime || getStep()%updateFreq==0 ){ firsttime=false; unlockContributors(); }
  }
}

void MultiColvarBase::updateCSphereArrays(){
  if( !usespecies || isDensity() ) return ;

  if( contributorsAreUnlocked ){
     if( !serialCalculation() ) comm.Sum( csphere_flags );
     unsigned istart=0;
     for(unsigned i=0;i<getCurrentNumberOfActiveTasks();++i){
         csphere_atoms[i].deactivateAll();
         for(unsigned j=0;j<csphere_atoms[i].fullSize();++j){
             if( csphere_flags[istart+j]==0 ) csphere_atoms[i].activate(j);
         }
         csphere_atoms[i].updateActiveMembers();
         istart += csphere_atoms[i].fullSize();
     }
     plumed_assert( istart==csphere_flags.size() );
  } else {
     for(unsigned i=0;i<csphere_flags.size();++i) csphere_flags[i]=0;
  } 
}

void MultiColvarBase::resizeLocalArrays(){
  atoms_with_derivatives.clear(); 
  for(unsigned i=0;i<getNumberOfAtoms();++i) atoms_with_derivatives.addIndexToList( i );
  atoms_with_derivatives.deactivateAll();
  // Set up stuff for central atoms
  atomsWithCatomDer.clear();
  for(unsigned i=0;i<getNumberOfAtoms();++i) atomsWithCatomDer.addIndexToList( i );
  atomsWithCatomDer.deactivateAll();
  // Resize tempory forces array
  if( !doNotCalculateDerivatives() ) forcesToApply.resize( getNumberOfDerivatives() );
  else forcesToApply.resize( 0 );
}

bool MultiColvarBase::setupCurrentAtomList( const unsigned& taskCode ){
  if( usespecies ){
     natomsper=1;
     current_atoms[0]=getBaseQuantityIndex( taskCode );
     if( contributorsAreUnlocked ){
        csphere_start=0; for(unsigned i=0;i<taskCode;++i) csphere_start+=csphere_atoms[i].fullSize();
     }
     for(unsigned j=0;j<ablocks.size();++j){
        for(unsigned i=0;i<csphere_atoms[taskCode].getNumberActive();++i){
           current_atoms[natomsper]=getBaseQuantityIndex( csphere_atoms[taskCode][i] );
           natomsper++; 
        }
     }
     if( natomsper==1 ) return isDensity();
  } else if( current_atoms.size()<4 ){
     natomsper=current_atoms.size();
     unsigned scode = taskCode;
     for(unsigned i=0;i<ablocks.size();++i){
        unsigned ind=std::floor( scode / decoder[i] );
        current_atoms[i]=getBaseQuantityIndex( ablocks[i][ind] );
        scode -= ind*decoder[i]; 
     }
  } else {
     natomsper=current_atoms.size(); 
     for(unsigned i=0;i<ablocks.size();++i) current_atoms[i]=getBaseQuantityIndex( ablocks[i][taskCode] );
  } 
  return true;
}

void MultiColvarBase::performTask(){
  // Currently no atoms have derivatives so deactivate those that are active
  atoms_with_derivatives.deactivateAll();
  // Currently no central atoms have derivatives so deactive them all
  atomsWithCatomDer.deactivateAll();
  // Retrieve the atom list
  if( !setupCurrentAtomList( getCurrentTask() ) ) return;

  // Do a quick check on the size of this contribution  
  calculateWeight();
  if( getElementValue(1)<getTolerance() ){
     updateActiveAtoms();
     return;   
  }

  // Compute everything
  double vv=doCalculation();
  // Set the value of this element in ActionWithVessel
  setElementValue( 0, vv );
  return;
}

double MultiColvarBase::doCalculation(){
  double val=compute(); updateActiveAtoms();
  return val;
}

Vector MultiColvarBase::retrieveCentralAtomPos(){
  ibox=getPbc().getInvBox().transpose();

  if( atomsWithCatomDer.getNumberActive()==0 ){
      Vector cvec = calculateCentralAtomPosition();
      for(unsigned i=0;i<3;++i) setElementValue( 2+i, cvec[i] );
      return cvec;
  }
  Vector cvec; 
  for(unsigned i=0;i<3;++i) cvec[i]=getElementValue( 2+i );
  return cvec;
}

void MultiColvarBase::addCentralAtomDerivatives( const unsigned& iatom, const Tensor& der ){
  plumed_dbg_assert( iatom<getNumberOfAtoms() );
  atomsWithCatomDer.activate(iatom);
  unsigned nder = 3*getNumberOfAtoms() + 9;
  for(unsigned i=0;i<3;++i){ 
    for(unsigned j=0;j<3;++j){
        addElementDerivative( (2+j)*nder + 3*iatom + i, der(j,i) );         
     }
  }
}

double MultiColvarBase::getCentralAtomDerivative( const unsigned& iatom, const unsigned& jcomp, const Vector& df ){
  plumed_dbg_assert( atomsWithCatomDer.isActive(iatom) && jcomp<3 );
  unsigned nder = 3*getNumberOfAtoms() + 9;
  return df[0]*getElementDerivative( 2*nder + 3*iatom + jcomp ) + 
         df[1]*getElementDerivative( 3*nder + 3*iatom + jcomp ) + 
         df[2]*getElementDerivative( 4*nder + 3*iatom + jcomp ); 
}

Vector MultiColvarBase::getSeparation( const Vector& vec1, const Vector& vec2 ) const {
  if(usepbc){ return pbcDistance( vec1, vec2 ); }
  else{ return delta( vec1, vec2 ); }
}

unsigned MultiColvarBase::getInternalIndex( const AtomNumber& iatom ) const {
  plumed_massert( usespecies && ablocks.size()==1, "This should only be used to interogate atom centered multicolvars");
  unsigned katom=0; bool found=false;
  for(unsigned i=0;i<ablocks[0].size();++i){
      if( all_atoms[ all_atoms.linkIndex(ablocks[0][i]) ]==iatom ){
         katom=i; found=true;
      }
  }
  if(!true) error("could not find required atom in any of the quantities calculated by the base multicolvar");
  return katom;
}

void MultiColvarBase::getIndexList( const unsigned& ntotal, const unsigned& jstore, const unsigned& maxder, std::vector<unsigned>& indices ){
  plumed_dbg_assert( !doNotCalculateDerivatives() );
  indices[jstore]=3*atoms_with_derivatives.getNumberActive() + 9;
  if( indices[jstore]>maxder ) error("too many derivatives to store. Run with LOWMEM");

  unsigned kder = ntotal + jstore*maxder;
  for(unsigned jder=0;jder<atoms_with_derivatives.getNumberActive();++jder){
     unsigned iatom = 3*atoms_with_derivatives[jder];
     for(unsigned icomp=0;icomp<3;++icomp){ indices[ kder ] = iatom+icomp; kder++; }
  }
  unsigned nbase = 3*getNumberOfAtoms(); 
  for(unsigned icomp=0;icomp<9;++icomp){ indices[ kder ] = nbase + icomp; kder++; }   
}   

void MultiColvarBase::getCentralAtomIndexList( const unsigned& ntotal, const unsigned& jstore, const unsigned& maxder, std::vector<unsigned>& indices ) const {
  plumed_dbg_assert( !doNotCalculateDerivatives() );

  indices[jstore]=3*atomsWithCatomDer.getNumberActive();
  if( indices[jstore]>maxder ) error("too many derivatives to store. Run with LOWMEM");

  unsigned kder  = ntotal + jstore*maxder;
  for(unsigned jder=0;jder<atomsWithCatomDer.getNumberActive();++jder){
     unsigned iatom = 3*atomsWithCatomDer[jder];
     for(unsigned icomp=0;icomp<3;++icomp){ indices[ kder ] = iatom+icomp; kder++; }
  }
}

void MultiColvarBase::activateIndexes( const unsigned& istart, const unsigned& number, const std::vector<unsigned>& indexes ){
  plumed_assert( number>0 );
  for(unsigned i=0;i<number-9;i+=3){
      plumed_dbg_assert( indexes[istart+i]%3==0 ); unsigned iatom=indexes[istart+i]/3; 
      atoms_with_derivatives.activate( iatom ); 
  }
}

void MultiColvarBase::quotientRule( const unsigned& uder, const unsigned& vder, const unsigned& iout ){
  unsigned ustart=uder*getNumberOfDerivatives();
  unsigned vstart=vder*getNumberOfDerivatives();
  unsigned istart=iout*getNumberOfDerivatives();
  double weight = getElementValue( vder ), pref = getElementValue( uder ) / (weight*weight);
  if( !doNotCalculateDerivatives() ){
      for(unsigned i=0;i<atoms_with_derivatives.getNumberActive();++i){
          unsigned n=3*atoms_with_derivatives[i], nx=n, ny=n+1, nz=n+2;
          setElementDerivative( istart + nx, getElementDerivative(ustart+nx) / weight - pref*getElementDerivative(vstart+nx) );
          setElementDerivative( istart + ny, getElementDerivative(ustart+ny) / weight - pref*getElementDerivative(vstart+ny) );
          setElementDerivative( istart + nz, getElementDerivative(ustart+nz) / weight - pref*getElementDerivative(vstart+nz) );
      }
      unsigned vbase=3*getNumberOfAtoms();
      for(unsigned i=0;i<9;++i){ 
          setElementDerivative( istart + vbase + i, getElementDerivative(ustart+vbase+i) / weight - pref*getElementDerivative(vstart+vbase+i) );
      }
  }
  thisval_wasset[iout]=false; setElementValue( iout, getElementValue(uder) / weight );
}

void MultiColvarBase::mergeDerivatives( const unsigned& ider, const double& df ){
  unsigned vstart=getNumberOfDerivatives()*ider;
  for(unsigned i=0;i<atoms_with_derivatives.getNumberActive();++i){
     unsigned iatom=3*atoms_with_derivatives[i];
     accumulateDerivative( iatom, df*getElementDerivative(vstart+iatom) ); iatom++;
     accumulateDerivative( iatom, df*getElementDerivative(vstart+iatom) ); iatom++;
     accumulateDerivative( iatom, df*getElementDerivative(vstart+iatom) );
  }
  unsigned nvir=3*getNumberOfAtoms();
  for(unsigned j=0;j<9;++j){
     accumulateDerivative( nvir, df*getElementDerivative(vstart+nvir) ); nvir++;
  }
}

void MultiColvarBase::clearDerivativesAfterTask( const unsigned& ider ){
  unsigned vstart=getNumberOfDerivatives()*ider;
  thisval_wasset[ider]=false; setElementValue( ider, 0.0 );
  thisval_wasset[ider]=false;
  if( ider>1 && ider<5 && derivativesAreRequired() ){
     for(unsigned i=0;i<atomsWithCatomDer.getNumberActive();++i){
        unsigned iatom=vstart+3*atomsWithCatomDer[i];
        setElementDerivative( iatom, 0.0 ); iatom++;
        setElementDerivative( iatom, 0.0 ); iatom++;
        setElementDerivative( iatom, 0.0 );
     }  
  } else if( derivativesAreRequired() ) {
     for(unsigned i=0;i<atoms_with_derivatives.getNumberActive();++i){
        unsigned iatom=vstart+3*atoms_with_derivatives[i];
        setElementDerivative( iatom, 0.0 ); iatom++;
        setElementDerivative( iatom, 0.0 ); iatom++;
        setElementDerivative( iatom, 0.0 );
     }   
     unsigned nvir=vstart+3*getNumberOfAtoms();
     for(unsigned j=0;j<9;++j){
        setElementDerivative( nvir, 0.0 ); nvir++;
     }
  }
}

void MultiColvarBase::apply(){
  if( getForcesFromVessels( forcesToApply ) ) setForcesOnAtoms( forcesToApply );
}

vesselbase::StoreDataVessel* MultiColvarBase::buildDataStashes(){
  // Check if vessels have already been setup
  for(unsigned i=0;i<getNumberOfVessels();++i){
     StoreColvarVessel* ssc=dynamic_cast<StoreColvarVessel*>( getPntrToVessel(i) );
     if(ssc) return ssc;
  }
 
  // Setup central atoms
  vesselbase::VesselOptions da("","",0,"",this);
  mycatoms=new StoreCentralAtomsVessel(da);
  addVessel(mycatoms);

  // Setup store values vessel
  vesselbase::VesselOptions ta("","",0,"",this);
  myvalues=new StoreColvarVessel(ta);   // Currently ignoring weights - good thing?
  addVessel(myvalues);

  // Make sure resizing of vessels is done
  resizeFunctions();
  return myvalues;
}


Vector MultiColvarBase::getCentralAtomPosition( const unsigned& iatom ) const {
  plumed_dbg_assert( mycatoms );
  return mycatoms->getPosition( iatom );
}

void MultiColvarBase::addCentralAtomDerivativeToFunction( const unsigned& iatom, const unsigned& jout, const unsigned& base_cv_no, const Vector& der, MultiColvarFunction* func ){
  plumed_dbg_assert( mycatoms ); 
  if( usingLowMem() ){
      mycatoms->recompute( iatom, 0 ); ;
      mycatoms->addAtomsDerivatives( 0, jout, base_cv_no, der, func );
  } else{
      mycatoms->addAtomsDerivatives( iatom, jout, base_cv_no, der, func ); 
  }
}

void MultiColvarBase::getValueForTask( const unsigned& iatom, std::vector<double>& vals ){
  plumed_dbg_assert( myvalues && vals.size()==1 );
  vals[0]=myvalues->getValue( iatom );
}

void MultiColvarBase::addWeightedValueDerivatives( const unsigned& iatom, const unsigned& base_cv_no, const double& weight, MultiColvarFunction* func ){
  plumed_dbg_assert( myvalues );
  if( usingLowMem() ){
     myvalues->recompute( iatom, 0 );
     myvalues->chainRuleForComponent( 0, 0, base_cv_no, weight, func );
  } else {
     myvalues->chainRuleForComponent( iatom, 0, base_cv_no, weight, func );
  }
}

void MultiColvarBase::finishWeightedAverageCalculation( MultiColvarFunction* func ){
  func->quotientRule( 0, 1, 0 );
}

void MultiColvarBase::addOrientationDerivativesToBase( const unsigned& iatom, const unsigned& jstore, const unsigned& base_cv_no, 
                                                       const std::vector<double>& weight, MultiColvarFunction* func ) {
  plumed_merror("This should not be called - invalid use of multicolvar in function");
} 
     
}
}
