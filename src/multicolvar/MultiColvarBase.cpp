/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2016 The plumed team
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
#include "MultiColvarBase.h"
#include "BridgedMultiColvarFunction.h"
#include "ActionVolume.h"
#include "MultiColvarFilter.h"
#include "vesselbase/Vessel.h"
#include "vesselbase/BridgeVessel.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "tools/Pbc.h"
#include "AtomValuePack.h"
#include "CatomPack.h"
#include "CatomPack.h"
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
allthirdblockintasks(false),
uselinkforthree(false),
linkcells(comm),
threecells(comm),
setup_completed(false),
atomsWereRetrieved(false),
usespecies(false)
{
  if( keywords.exists("NOPBC") ){ 
    bool nopbc=!usepbc; parseFlag("NOPBC",nopbc);
    usepbc=!nopbc;
  } 
  if( keywords.exists("SPECIES") ) usespecies=true;
}

bool MultiColvarBase::interpretInputMultiColvars( const std::vector<std::string>& mlabs, const double& wtolerance ){
  if( mlabs.size()==0 ) return false;

  std::string mname;
  for(unsigned i=0;i<mlabs.size();++i){
      MultiColvarBase* mycolv = plumed.getActionSet().selectWithLabel<MultiColvarBase*>(mlabs[i]);
      if(!mycolv) return false;
      // Check all base multicolvars are of same type
      if( i==0 ){
          mname = mycolv->getName();
          if( mycolv->isPeriodic() ) error("multicolvar functions don't work with this multicolvar");
      } else {
          if( mname!=mycolv->getName() ) error("All input multicolvars must be of same type");
      }
      // And track which variable stores each colvar
      for(unsigned j=0;j<mycolv->getFullNumberOfTasks();++j) colvar_label.push_back( mybasemulticolvars.size() );
      // And store the multicolvar base
      mybasemulticolvars.push_back( mycolv );
      // And create a basedata stash
      mybasedata.push_back( mybasemulticolvars[mybasemulticolvars.size()-1]->buildDataStashes( true, wtolerance, this ) );
      // Check if weight has derivatives
      if( mybasemulticolvars[ mybasemulticolvars.size()-1 ]->weightHasDerivatives ) weightHasDerivatives=true;    
      plumed_assert( mybasemulticolvars.size()==mybasedata.size() );
  }

  log.printf("  using colvars calculated by actions "); 
  for(unsigned i=0;i<mlabs.size();++i) log.printf("%s ",mlabs[i].c_str() );
  log.printf("\n"); 
  return true;
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

void MultiColvarBase::setupMultiColvarBase( const std::vector<AtomNumber>& atoms ){
  // Setup decoder array
  if( !usespecies && ablocks.size()<4 ){
     use_for_central_atom.resize( ablocks.size(), true );
     numberForCentralAtom = 1.0 / static_cast<double>( ablocks.size() );
     if( ablocks.size()==3 ){
         allthirdblockintasks=uselinkforthree=true;
         for(unsigned i=0;i<bookeeping.nrows();++i){
             for(unsigned j=0;j<bookeeping.ncols();++j){
                 unsigned ntper = bookeeping(i,j).second - bookeeping(i,j).first;
                 if( i==j && ntper==0 ){
                     continue;
                 } else if( ntper == 1 && allthirdblockintasks ){
                     allthirdblockintasks=true;
                 } else if( ntper != ablocks[2].size() ){
                     allthirdblockintasks=uselinkforthree=false;
                 } else {
                     allthirdblockintasks=false;
                 }
             }
         } 
     }
    
     if( allthirdblockintasks ) decoder.resize(2);
     else decoder.resize( ablocks.size() ); 
     unsigned code=1; for(unsigned i=0;i<decoder.size();++i){ decoder[decoder.size()-1-i]=code; code *= nblock; }
  } else if( !usespecies ){
     use_for_central_atom.resize( ablocks.size(), true );
     numberForCentralAtom = 1.0 / static_cast<double>( ablocks.size() );
  }

  // Copy lists of atoms involved from base multicolvars 
  std::vector<AtomNumber> tmp_atoms, all_atoms;
  for(unsigned i=0;i<mybasemulticolvars.size();++i){
      BridgedMultiColvarFunction* mybr=dynamic_cast<BridgedMultiColvarFunction*>( mybasemulticolvars[i] );
      if( mybr ) tmp_atoms=(mybr->getPntrToMultiColvar())->getAbsoluteIndexes();
      else tmp_atoms=mybasemulticolvars[i]->getAbsoluteIndexes();
      for(unsigned j=0;j<tmp_atoms.size();++j) all_atoms.push_back( tmp_atoms[j] );
  } 
  // Get additional atom requests
  for(unsigned i=0;i<atoms.size();++i) all_atoms.push_back( atoms[i] );

  // Now make sure we get all the atom positions 
  ActionAtomistic::requestAtoms( all_atoms );
  // And setup dependencies
  for(unsigned i=0;i<mybasemulticolvars.size();++i) addDependency( mybasemulticolvars[i] );

  // Setup underlying ActionWithVessel
  readVesselKeywords();
}

void MultiColvarBase::setAtomsForCentralAtom( const std::vector<bool>& catom_ind ){
  unsigned nat=0; plumed_assert( catom_ind.size()==ablocks.size() );
  for(unsigned i=0;i<catom_ind.size();++i){
      use_for_central_atom[i]=catom_ind[i]; 
      if( use_for_central_atom[i] ) nat++;
  }
  plumed_dbg_assert( nat>0 ); numberForCentralAtom = 1.0 / static_cast<double>( nat );
}

void MultiColvarBase::turnOnDerivatives(){
  ActionWithValue::turnOnDerivatives();
  needsDerivatives(); 
  forcesToApply.resize( getNumberOfDerivatives() );
} 

void MultiColvarBase::setLinkCellCutoff( const double& lcut, double tcut ){
  plumed_assert( usespecies || ablocks.size()<4 );
  if( tcut<0 ) tcut=lcut;
  linkcells.setCutoff( lcut ); 
  threecells.setCutoff( tcut );
}

void MultiColvarBase::setupLinkCells(){
  if( !linkcells.enabled() ) return ;
  // Retrieve any atoms that haven't already been retrieved
  for(std::vector<MultiColvarBase*>::iterator p=mybasemulticolvars.begin();p!=mybasemulticolvars.end();++p){
     (*p)->retrieveAtoms();
  }
  retrieveAtoms();

  unsigned iblock;
  if( usespecies ){
      iblock=0; 
  } else if( ablocks.size()<4 ){ 
      iblock=1;  
  } else {
      plumed_error();
  }
 
  // Count number of currently active atoms
  unsigned nactive_atoms=0;
  for(unsigned i=0;i<ablocks[iblock].size();++i){
      if( isCurrentlyActive( iblock, ablocks[iblock][i] ) ) nactive_atoms++;
  }

  std::vector<Vector> ltmp_pos( nactive_atoms ); 
  std::vector<unsigned> ltmp_ind( nactive_atoms );

  nactive_atoms=0;
  if( usespecies ){
     for(unsigned i=0;i<ablocks[0].size();++i){
        if( !isCurrentlyActive( 0, ablocks[0][i] ) ) continue; 
        ltmp_ind[nactive_atoms]=ablocks[0][i];
        ltmp_pos[nactive_atoms]=getPositionOfAtomForLinkCells( ltmp_ind[nactive_atoms] );
        nactive_atoms++;
     }
  } else {
     for(unsigned i=0;i<ablocks[1].size();++i){
        if( !isCurrentlyActive( 1, ablocks[1][i] ) ) continue;
        ltmp_ind[nactive_atoms]=i; 
        ltmp_pos[nactive_atoms]=getPositionOfAtomForLinkCells( ablocks[1][i] );
        nactive_atoms++; 
     }
  }

  // Build the lists for the link cells
  linkcells.buildCellLists( ltmp_pos, ltmp_ind, getPbc() );
}

void MultiColvarBase::setupNonUseSpeciesLinkCells( const unsigned& my_always_active ){
  plumed_assert( !usespecies );
  if( !linkcells.enabled() ) return ;

  if( !uselinkforthree ){
     // Get some parallel info
     unsigned stride=comm.Get_size();
     unsigned rank=comm.Get_rank(); 
     if( serialCalculation() ){ stride=1; rank=0; }

     // Ensure we only do tasks where atoms are in appropriate link cells
     std::vector<unsigned> linked_atoms( 1+ablocks[1].size() ); deactivateAllTasks();
     for(unsigned i=rank;i<ablocks[0].size();i+=stride){
         if( !isCurrentlyActive( 0, ablocks[0][i] ) ) continue;
         unsigned natomsper=1; linked_atoms[0]=my_always_active;  // Note we always check atom 0 because it is simpler than changing LinkCells.cpp
         linkcells.retrieveNeighboringAtoms( getPositionOfAtomForLinkCells( ablocks[0][i] ), natomsper, linked_atoms );
         for(unsigned j=0;j<natomsper;++j){
             for(unsigned k=bookeeping(i,linked_atoms[j]).first;k<bookeeping(i,linked_atoms[j]).second;++k) taskFlags[k]=1;
         }
     }
     if( !serialCalculation() ) comm.Sum( taskFlags ); 
     lockContributors();
  } else { 
     // Get some parallel info
     unsigned stride=comm.Get_size();
     unsigned rank=comm.Get_rank();
     if( serialCalculation() ){ stride=1; rank=0; }

     unsigned nactive_three=0;
     for(unsigned i=0;i<ablocks[2].size();++i){
         if( isCurrentlyActive( 2, ablocks[2][i] ) ) nactive_three++;
     }

     std::vector<Vector> lttmp_pos( nactive_three );
     std::vector<unsigned> lttmp_ind( nactive_three );

     nactive_three=0;
     if( allthirdblockintasks ){
         for(unsigned i=0;i<ablocks[2].size();++i){
             if( !isCurrentlyActive( 2, ablocks[2][i] ) ) continue;
             lttmp_ind[nactive_three]=ablocks[2][i];
             lttmp_pos[nactive_three]=getPositionOfAtomForLinkCells( ablocks[2][i] );
             nactive_three++;
         }
     } else {
         for(unsigned i=0;i<ablocks[2].size();++i){
             if( !isCurrentlyActive( 2, ablocks[2][i] ) ) continue;
             lttmp_ind[nactive_three]=i;
             lttmp_pos[nactive_three]=getPositionOfAtomForLinkCells( ablocks[2][i] );
             nactive_three++;
         }
     }
     // Build the list of the link cells
     threecells.buildCellLists( lttmp_pos, lttmp_ind, getPbc() );

     // Ensure we only do tasks where atoms are in appropriate link cells
     deactivateAllTasks();
     std::vector<unsigned> linked_atoms( 1+ablocks[1].size() );
     std::vector<unsigned> tlinked_atoms( 1+ablocks[2].size() );
     for(unsigned i=rank;i<ablocks[0].size();i+=stride){
         if( !isCurrentlyActive( 0, ablocks[0][i] ) ) continue;
         unsigned natomsper=1; linked_atoms[0]=my_always_active;  // Note we always check atom 0 because it is simpler than changing LinkCells.cpp
         linkcells.retrieveNeighboringAtoms( getPositionOfAtomForLinkCells( ablocks[0][i] ), natomsper, linked_atoms );
         if( allthirdblockintasks ) {
             for(unsigned j=0;j<natomsper;++j){
                 for(unsigned k=bookeeping(i,linked_atoms[j]).first;k<bookeeping(i,linked_atoms[j]).second;++k) taskFlags[k]=1;
             }
         } else {
             unsigned ntatomsper=1; tlinked_atoms[0]=lttmp_ind[0];
             threecells.retrieveNeighboringAtoms( getPositionOfAtomForLinkCells( ablocks[0][i] ), ntatomsper, tlinked_atoms );
             for(unsigned j=0;j<natomsper;++j){
                 for(unsigned k=0;k<ntatomsper;++k) taskFlags[bookeeping(i,linked_atoms[j]).first+tlinked_atoms[k]]=1;
             }
         }
     }
     if( !serialCalculation() ) comm.Sum( taskFlags );
     lockContributors();
  } 
}

void MultiColvarBase::decodeIndexToAtoms( const unsigned& taskCode, std::vector<unsigned>& atoms ) const {
  plumed_dbg_assert( !usespecies && ablocks.size()<4 );
  if( atoms.size()!=decoder.size() ) atoms.resize( decoder.size() );

  unsigned scode = taskCode;
  for(unsigned i=0;i<decoder.size();++i){
      unsigned ind=( scode / decoder[i] );
      atoms[i] = ablocks[i][ind];
      scode -= ind*decoder[i];
  }
}

bool MultiColvarBase::setupCurrentAtomList( const unsigned& taskCode, AtomValuePack& myatoms ) const {
  if( usespecies ){
     if( isDensity() ) return true;
     std::vector<unsigned> task_atoms(1); task_atoms[0]=taskCode;
     unsigned natomsper=myatoms.setupAtomsFromLinkCells( task_atoms, getLinkCellPosition(task_atoms), linkcells );
     return natomsper>1;
  } else if( allthirdblockintasks ){ 
     plumed_dbg_assert( ablocks.size()==3 ); std::vector<unsigned> atoms(2); decodeIndexToAtoms( taskCode, atoms );
     unsigned natomsper=myatoms.setupAtomsFromLinkCells( atoms, getLinkCellPosition(atoms), threecells );
  } else if( ablocks.size()<4 ){
     std::vector<unsigned> atoms( ablocks.size() );
     decodeIndexToAtoms( taskCode, atoms ); myatoms.setNumberOfAtoms( ablocks.size() );
     for(unsigned i=0;i<ablocks.size();++i) myatoms.setAtom( i, atoms[i] ); 
  } else {
     myatoms.setNumberOfAtoms( ablocks.size() );
     for(unsigned i=0;i<ablocks.size();++i) myatoms.setAtom( i, ablocks[i][taskCode] ); 
  } 
  return true;
}

void MultiColvarBase::setupActiveTaskSet( std::vector<unsigned>& active_tasks, const std::string& input_label ){
  if( !setup_completed ){ 
      bool justVolumes=false;
      if( usespecies ){
          justVolumes=true;
          for(unsigned i=0;i<getNumberOfVessels();++i){
              vesselbase::StoreDataVessel* mys=dynamic_cast<vesselbase::StoreDataVessel*>( getPntrToVessel(i) );
              if( mys ) continue;
              vesselbase::BridgeVessel* myb=dynamic_cast<vesselbase::BridgeVessel*>( getPntrToVessel(i) );
              if( !myb ){ justVolumes=false; break; }
              ActionVolume* myv=dynamic_cast<ActionVolume*>( myb->getOutputAction() );
              if( !myv ){ justVolumes=false; break; }
          }
      }
      deactivateAllTasks();
      if( justVolumes && mydata ){
          if( mydata->getNumberOfDataUsers()==0 ) justVolumes=false;

          for(unsigned i=0;i<mydata->getNumberOfDataUsers();++i){
              MultiColvarBase* myu=dynamic_cast<MultiColvarBase*>( mydata->getDataUser(i) );
              if( myu ){
                  myu->setupActiveTaskSet( taskFlags, getLabel() );
              } else {
                  for(unsigned i=0;i<getFullNumberOfTasks();++i) taskFlags[i]=1;
              }
          }
      }
      if( justVolumes ){
          for(unsigned j=0;j<getNumberOfVessels();++j){
              vesselbase::BridgeVessel* myb=dynamic_cast<vesselbase::BridgeVessel*>( getPntrToVessel(j) );
              if( !myb ) continue ;
              ActionVolume* myv=dynamic_cast<ActionVolume*>( myb->getOutputAction() );
              if( !myv ) continue ;
              myv->retrieveAtoms(); myv->setupRegions();
              
              for(unsigned i=0;i<getFullNumberOfTasks();++i){
                 if( myv->inVolumeOfInterest(i) ) taskFlags[i]=1;
              }
          }
      } else { 
          for(unsigned i=0;i<getFullNumberOfTasks();++i) taskFlags[i]=1;
      } 

      // Now activate all this class
      lockContributors();
      // Setup the link cells
      setupLinkCells();  
      // Ensures that setup is not performed multiple times during one cycle
      setup_completed=true;
  }

  // And activate the tasks in input action
  if( getLabel()!=input_label ){
      int input_code=-1;
      for(unsigned i=0;i<mybasemulticolvars.size();++i){
          if( mybasemulticolvars[i]->getLabel()==input_label ){ input_code=i; break; }
      }

      MultiValue my_tvals( getNumberOfQuantities(), getNumberOfDerivatives() ); 
      AtomValuePack mytmp_atoms( my_tvals, this );   
      for(unsigned i=0;i<getFullNumberOfTasks();++i){
          if( !taskIsCurrentlyActive(i) ) continue;
          setupCurrentAtomList( getTaskCode(i), mytmp_atoms );
          for(unsigned j=0;j<mytmp_atoms.getNumberOfAtoms();++j){
              unsigned itask=mytmp_atoms.getIndex(j);
              if( colvar_label[itask]==input_code ) active_tasks[ convertToLocalIndex( itask, input_code ) ]=1;   
          }
      }
  }
}

bool MultiColvarBase::filtersUsedAsInput(){
  bool inputAreFilters=false;
  for(unsigned i=0;i<mybasemulticolvars.size();++i){
      MultiColvarFilter* myfilt=dynamic_cast<MultiColvarFilter*>( mybasemulticolvars[i] );
      if( myfilt || mybasemulticolvars[i]->filtersUsedAsInput() ) inputAreFilters=true;
  }
  return inputAreFilters;
}

void MultiColvarBase::calculate(){ 
  // Recursive function that sets up tasks
  setupActiveTaskSet( taskFlags, getLabel() );

  // Check for filters and rerun setup of link cells if there are any
  if( colvar_label.size()>0 && filtersUsedAsInput() ) setupLinkCells();

  //  Setup the link cells if we are not using species
  if( !usespecies && ablocks.size()>1 ){
     // This loop finds the first active atom, which is always checked because
     // of a peculiarity in linkcells
     unsigned first_active;
     for(unsigned i=0;i<ablocks[0].size();++i){
        if( !isCurrentlyActive( 1, ablocks[1][i] ) ) continue;
        else {
           first_active=i; break;
        }
     }
     setupNonUseSpeciesLinkCells( first_active );
  }
  // And run all tasks
  runAllTasks();
}

void MultiColvarBase::calculateNumericalDerivatives( ActionWithValue* a ){
  if( colvar_label.size()>0 ) plumed_merror("cannot calculate numerical derivatives for this quantity");
  calculateAtomicNumericalDerivatives( this, 0 );
}

void MultiColvarBase::prepare(){
  setup_completed=false; atomsWereRetrieved=false;
}

void MultiColvarBase::retrieveAtoms(){
  if( !atomsWereRetrieved ){ ActionAtomistic::retrieveAtoms(); atomsWereRetrieved=true; }
}

void MultiColvarBase::addAtomDerivatives( const int& ival, const unsigned& iatom, const Vector& der, multicolvar::AtomValuePack& myatoms ) const {
  unsigned jatom=myatoms.getIndex(iatom);

  if( jatom<colvar_label.size() ){
      unsigned mmc=colvar_label[jatom];
      unsigned basen=0; for(unsigned i=0;i<mmc;++i) basen+=mybasemulticolvars[i]->getNumberOfAtoms();
      multicolvar::CatomPack atom0=mybasemulticolvars[mmc]->getCentralAtomPack( basen, convertToLocalIndex(jatom,mmc) );
      myatoms.addComDerivatives( ival, der, atom0 );
  } else {
      if( ival<0 ) myatoms.addTemporyAtomsDerivatives( iatom, der );
      else myatoms.addAtomsDerivatives( ival, iatom, der );
  }
}

void MultiColvarBase::performTask( const unsigned& task_index, const unsigned& current, MultiValue& myvals ) const {

  AtomValuePack myatoms( myvals, this );
  // Retrieve the atom list
  if( !setupCurrentAtomList( current, myatoms ) ) return;

  // Do a quick check on the size of this contribution  
  calculateWeight( current, myatoms ); 
  if( myatoms.getValue(0)<getTolerance() ){
     updateActiveAtoms( myatoms );
     return;   
  }

  // Compute everything
  double vv=doCalculation( task_index, myatoms ); 
  myatoms.setValue( 1, vv );
  return;
}

void MultiColvarBase::calculateWeight( const unsigned& taskCode, AtomValuePack& myatoms ) const {
  if( usespecies && taskCode<colvar_label.size() ){
      unsigned mmc=colvar_label[taskCode]; std::vector<double> old_data( mybasemulticolvars[mmc]->getNumberOfQuantities() );
      plumed_dbg_assert( mybasedata[mmc]->storedValueIsActive( convertToLocalIndex(taskCode,mmc) ) );
      mybasedata[mmc]->retrieveValueWithIndex( convertToLocalIndex(taskCode,mmc), false, old_data );
      myatoms.setValue( 0, old_data[0] );
      if( !doNotCalculateDerivatives() && mybasemulticolvars[mmc]->weightHasDerivatives ){
          MultiValue myder( mybasemulticolvars[mmc]->getNumberOfQuantities(), mybasemulticolvars[mmc]->getNumberOfDerivatives() );
          MultiValue& outder=myatoms.getUnderlyingMultiValue(); mybasedata[mmc]->retrieveDerivatives( convertToLocalIndex(taskCode,mmc), false, myder );
          for(unsigned j=0;j<myder.getNumberActive();++j){ unsigned jder=myder.getActiveIndex(j); outder.addDerivative( 0, jder, myder.getDerivative(0,jder) ); }
      }
  } else {
      myatoms.setValue( 0, 1.0 );
  }
}

double MultiColvarBase::doCalculation( const unsigned& taskIndex, AtomValuePack& myatoms ) const {
  if( colvar_label.size()>0 ) mybasedata[0]->resetTemporyMultiValues();
  double val=compute( taskIndex, myatoms ); updateActiveAtoms( myatoms );
  return val;
}

void MultiColvarBase::updateActiveAtoms( AtomValuePack& myatoms ) const {
  if( colvar_label.size()==0 ) myatoms.updateUsingIndices();
  else myatoms.updateDynamicList();
}

Vector MultiColvarBase::getCentralAtomPos( const unsigned& taskIndex ){
  unsigned curr=getTaskCode( taskIndex );

  if( usespecies || isDensity() ){
     return getPositionOfAtomForLinkCells(curr);
  } else if( ablocks.size()<4 ){
     // double factor=1.0/static_cast<double>( ablocks.size() );
     Vector mypos; mypos.zero(); 
     std::vector<unsigned> atoms( ablocks.size() ); decodeIndexToAtoms( curr, atoms );
     for(unsigned i=0;i<ablocks.size();++i){
         if( use_for_central_atom[i] ) mypos+=numberForCentralAtom*getPositionOfAtomForLinkCells(atoms[i]); 
     }
     return mypos;
  } else {
     Vector mypos; mypos.zero();
     for(unsigned i=0;i<ablocks.size();++i){
         if( use_for_central_atom[i] ) mypos+=numberForCentralAtom*getPositionOfAtomForLinkCells(ablocks[i][curr]);
     }
     return mypos;
  }
}

CatomPack MultiColvarBase::getCentralAtomPack( const unsigned& basn, const unsigned& taskIndex ){
  unsigned curr=getTaskCode( taskIndex );

  CatomPack mypack;
  if(usespecies){
     mypack.resize(1);
     mypack.setIndex( 0, basn + curr );
     mypack.setDerivative( 0, Tensor::identity() );
  } else if( ablocks.size()<4 ){
     mypack.resize(ablocks.size());
     unsigned k=0;
     std::vector<unsigned> atoms( ablocks.size() ); decodeIndexToAtoms( curr, atoms );
     for(unsigned i=0;i<ablocks.size();++i){
         if( use_for_central_atom[i] ){
             mypack.setIndex( k, basn + atoms[i] );
             mypack.setDerivative( k, numberForCentralAtom*Tensor::identity() );
             k++;
         }
     }
  } else {
     unsigned k=0;
     for(unsigned i=0;i<ablocks.size();++i){
         if( use_for_central_atom[i] ){
             mypack.setIndex( k, basn + ablocks[i][curr] );
             mypack.setDerivative( k, numberForCentralAtom*Tensor::identity() );
             k++;
         }
     }
  }
  return mypack;
} 

Vector MultiColvarBase::getSeparation( const Vector& vec1, const Vector& vec2 ) const {
  if(usepbc){ return pbcDistance( vec1, vec2 ); }
  else{ return delta( vec1, vec2 ); }
}

void MultiColvarBase::applyPbc(std::vector<Vector>& dlist, unsigned int max_index) const {
   if (usepbc) pbcApply(dlist, max_index);
}

void MultiColvarBase::apply(){
  if( getForcesFromVessels( forcesToApply ) ) setForcesOnAtoms( forcesToApply );
}
     
}
}
