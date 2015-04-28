/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
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
#include "BridgedMultiColvarFunction.h"
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
  keys.setComponentsIntroduction("When the label of this action is used as the input for a second you are not referring to a scalar quantity as you are in "
                                 "regular collective variables.  The label is used to reference the full set of quantities calculated by "
                                 "the action.  This is usual when using \\ref multicolvarfunction. Generally when doing this the previously calculated "
                                 "multicolvar will be referenced using the DATA keyword rather than ARG.\n\n"
                                 "This Action can be used to calculate the following scalar quantities directly.  These quantities are calculated by "
                                 "employing the keywords listed below. "
                                 "These quantities can then be referenced elsewhere in the input file by using this Action's label "
                                 "followed by a dot and the name of the quantity. Some amongst them can be calculated multiple times "
                                 "with different parameters.  In this case the quantities calculated can be referenced elsewhere in the "
                                 "input by using the name of the quantity followed by a numerical identifier "
                                 "e.g. <em>label</em>.lessthan-1, <em>label</em>.lessthan-2 etc.  When doing this and, for clarity we have "
                                 "made the label of the components customizable. As such by using the LABEL keyword in the description of the keyword "
                                 "input you can customize the component name");
} 

MultiColvarBase::MultiColvarBase(const ActionOptions&ao):
Action(ao),
ActionAtomistic(ao),
ActionWithValue(ao),
ActionWithVessel(ao),
usepbc(false),
linkcells(comm),
mycatoms(NULL),        // This will be destroyed by ActionWithVesel
myvalues(NULL),        // This will be destroyed by ActionWithVesel 
usespecies(false)
{
  if( keywords.exists("NOPBC") ){ 
    bool nopbc=!usepbc; parseFlag("NOPBC",nopbc);
    usepbc=!nopbc;
  } 
  if( keywords.exists("SPECIES") ) usespecies=true;
}

void MultiColvarBase::addTaskToList( const unsigned& taskCode ){
  plumed_assert( getNumberOfVessels()==0 );
  ActionWithVessel::addTaskToList( taskCode );
}

void MultiColvarBase::resizeBookeepingArray( const unsigned& num1, const unsigned& num2 ){
  bookeeping.resize( num1, num2 );
  for(unsigned i=0;i<num1;++i){
      for(unsigned j=0;j<num2;++j){ bookeeping(i,j).first=0; bookeeping(i,j).second=0; }
  }
}

void MultiColvarBase::setupMultiColvarBase(){
  // Setup decoder array
  if( !usespecies && ablocks.size()<4 ){
     decoder.resize( ablocks.size() ); unsigned code=1;
     for(unsigned i=0;i<ablocks.size();++i){ decoder[ablocks.size()-1-i]=code; code *= nblock; } 
  } 
  // Resize stuff here
  resizeLocalArrays();

  // Setup underlying ActionWithVessel
  readVesselKeywords();
}

void MultiColvarBase::turnOnDerivatives(){
  ActionWithValue::turnOnDerivatives();
  needsDerivatives(); 
  forcesToApply.resize( getNumberOfDerivatives() );
} 

void MultiColvarBase::setLinkCellCutoff( const double& lcut ){
  plumed_assert( usespecies || current_atoms.size()<4 );
  linkcells.setCutoff( lcut );
}

void MultiColvarBase::setupLinkCells(){
  if( !linkcells.enabled() ) return ;

  unsigned iblock, jblock;
  if( usespecies ){
      iblock=0; 
  } else if( current_atoms.size()<4 ){ 
      iblock=1;  
  } else {
      plumed_error();
  }
 
  // Count number of currently active atoms
  unsigned nactive_atoms=0;
  for(unsigned i=0;i<ablocks[iblock].size();++i){
      if( isCurrentlyActive( ablocks[iblock][i] ) ) nactive_atoms++;
  }

  std::vector<Vector> ltmp_pos( nactive_atoms ); 
  std::vector<unsigned> ltmp_ind( nactive_atoms );

  nactive_atoms=0;
  if( usespecies ){
     ltmp_pos.resize( ablocks[0].size() ); ltmp_ind.resize( ablocks[0].size() );
     for(unsigned i=0;i<ablocks[0].size();++i){
        if( !isCurrentlyActive( ablocks[0][i] ) ) continue; 
        ltmp_ind[nactive_atoms]=ablocks[0][i];
        ltmp_pos[nactive_atoms]=getPositionOfAtomForLinkCells( ltmp_ind[nactive_atoms] );
        nactive_atoms++;
     }
  } else {
     ltmp_pos.resize( ablocks[1].size() ); ltmp_ind.resize( ablocks[1].size() ); 
     for(unsigned i=0;i<ablocks[1].size();++i){
        if( !isCurrentlyActive( ablocks[1][i] ) ) continue;
        ltmp_ind[nactive_atoms]=i; 
        ltmp_pos[nactive_atoms]=getPositionOfAtomForLinkCells( ablocks[1][i] );
        nactive_atoms++; 
     }
  }

  // Build the lists for the link cells
  linkcells.buildCellLists( ltmp_pos, ltmp_ind, getPbc() );

  if( !usespecies ){
     // Get some parallel info
     unsigned stride=comm.Get_size();
     unsigned rank=comm.Get_rank(); 
     if( serialCalculation() ){ stride=1; rank=0; }

     // Ensure we only do tasks where atoms are in appropriate link cells
     std::vector<unsigned> linked_atoms( 1+ablocks[1].size() );
     std::vector<unsigned>  active_tasks( getFullNumberOfTasks(), 0 );
     for(unsigned i=rank;i<ablocks[0].size();i+=stride){
         if( !isCurrentlyActive( ablocks[0][i] ) ) continue;
         natomsper=1; linked_atoms[0]=ltmp_ind[0];  // Note we always check atom 0 because it is simpler than changing LinkCells.cpp
         linkcells.retrieveNeighboringAtoms( getPositionOfAtomForLinkCells( ablocks[0][i] ), natomsper, linked_atoms );
         for(unsigned j=0;j<natomsper;++j){
             for(unsigned k=bookeeping(i,linked_atoms[j]).first;k<bookeeping(i,linked_atoms[j]).second;++k) active_tasks[k]=1;
         }
     }
     if( !serialCalculation() ) comm.Sum( active_tasks ); 

     deactivateAllTasks(); 
     activateTheseTasks( active_tasks );
     contributorsAreUnlocked=false;
  }
}

void MultiColvarBase::resizeLocalArrays(){
  atoms_with_derivatives.clear(); 
  for(unsigned i=0;i<getSizeOfAtomsWithDerivatives();++i) atoms_with_derivatives.addIndexToList( i );
  atoms_with_derivatives.deactivateAll();
  // Set up stuff for central atoms
  atomsWithCatomDer.clear();
  for(unsigned i=0;i<getSizeOfAtomsWithDerivatives();++i) atomsWithCatomDer.addIndexToList( i );
  atomsWithCatomDer.deactivateAll();
}

bool MultiColvarBase::setupCurrentAtomList( const unsigned& taskCode ){
  if( usespecies ){
     natomsper=1;
     if( isDensity() ) return true;
     current_atoms[0]=taskCode;
     linkcells.retrieveNeighboringAtoms( getPositionOfAtomForLinkCells(current_atoms[0]), natomsper, current_atoms );
     return natomsper>1;
  } else if( current_atoms.size()<4 ){
     natomsper=current_atoms.size();
     unsigned scode = taskCode;
     for(unsigned i=0;i<ablocks.size();++i){
        unsigned ind=( scode / decoder[i] );
        current_atoms[i]=ablocks[i][ind];
        scode -= ind*decoder[i]; 
     }
  } else {
     natomsper=current_atoms.size(); 
     for(unsigned i=0;i<ablocks.size();++i) current_atoms[i]=ablocks[i][taskCode];
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
  calculateWeight(); // printf("HELLO WEIGHT %f \n",getElementValue(1) );
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
  if( atomsWithCatomDer.getNumberActive()==0 ){
      Vector cvec = calculateCentralAtomPosition();
      for(unsigned i=0;i<3;++i) setElementValue( getCentralAtomElementIndex()+i, cvec[i] );
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
        addElementDerivative( (getCentralAtomElementIndex()+j)*nder + 3*iatom + i, der(j,i) );
     }
  }
}

double MultiColvarBase::getCentralAtomDerivative( const unsigned& iatom, const unsigned& jcomp, const Vector& df ){
  plumed_dbg_assert( atomsWithCatomDer.isActive(iatom) && jcomp<3 );
  unsigned nder = 3*getNumberOfAtoms() + 9;
  return df[0]*getElementDerivative( (getCentralAtomElementIndex()+0)*nder + 3*iatom + jcomp ) +
         df[1]*getElementDerivative( (getCentralAtomElementIndex()+1)*nder + 3*iatom + jcomp ) +
         df[2]*getElementDerivative( (getCentralAtomElementIndex()+2)*nder + 3*iatom + jcomp ); 
}

Vector MultiColvarBase::getSeparation( const Vector& vec1, const Vector& vec2 ) const {
  if(usepbc){ return pbcDistance( vec1, vec2 ); }
  else{ return delta( vec1, vec2 ); }
}

void MultiColvarBase::applyPbc(std::vector<Vector>& dlist, unsigned int max_index) const {
   if (usepbc) pbcApply(dlist, max_index);
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

vesselbase::StoreDataVessel* MultiColvarBase::buildDataStashes( const bool& allow_wcutoff, const double& wtol ){
  // Check if vessels have already been setup
  for(unsigned i=0;i<getNumberOfVessels();++i){
     StoreColvarVessel* ssc=dynamic_cast<StoreColvarVessel*>( getPntrToVessel(i) );
     if(ssc){
        if( allow_wcutoff && !ssc->weightCutoffIsOn() ) error("Cannot have more than one data stash with different properties");
        if( !allow_wcutoff && ssc->weightCutoffIsOn() ) error("Cannot have more than one data stash with different properties");
        return ssc;
     }
  }
 
  // Setup central atoms
  vesselbase::VesselOptions da("","",0,"",this);
  mycatoms=new StoreCentralAtomsVessel(da);
  if( allow_wcutoff ) mycatoms->setHardCutoffOnWeight( wtol );
  addVessel(mycatoms); 

  // Setup store values vessel
  vesselbase::VesselOptions ta("","",0,"",this);
  myvalues=new StoreColvarVessel(ta);   
  if( allow_wcutoff ) myvalues->setHardCutoffOnWeight( wtol );
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
  // plumed_dbg_assert( myvalues && vals.size()==1 );  // Problem if we want to do MultiColvar + Bridge + MultiColvarFunction
  vals[0]=myvalues->getValue( iatom );
}

void MultiColvarBase::copyElementsToBridgedColvar( BridgedMultiColvarFunction* func ){
  func->setElementValue( 0, getElementValue(0) ); 
  for(unsigned i=0;i<atoms_with_derivatives.getNumberActive();++i){
     unsigned n=atoms_with_derivatives[i], nx=3*n;
     func->atoms_with_derivatives.activate(n);
     func->addElementDerivative( nx+0, getElementDerivative(nx+0) );
     func->addElementDerivative( nx+1, getElementDerivative(nx+1) );
     func->addElementDerivative( nx+2, getElementDerivative(nx+2) ); 
  }
  unsigned nvir=3*getNumberOfAtoms();
  for(unsigned i=0;i<9;++i){ 
     func->addElementDerivative( nvir, getElementDerivative(nvir) ); nvir++;
  }
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

bool MultiColvarBase::storedValueIsActive( const unsigned& iatom ){
  return myvalues->storedValueIsActive( iatom );
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
