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
#include "MultiColvarFunction.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "tools/Pbc.h"
#include "MultiColvarBase.h"

namespace PLMD {
namespace multicolvar { 

void MultiColvarFunction::registerKeywords( Keywords& keys ){
  MultiColvarBase::registerKeywords( keys );
  keys.add("compulsory","DATA","the labels of the action that calculates the multicolvars we are interested in");
  keys.reset_style("NUMERICAL_DERIVATIVES","hidden");
}

MultiColvarFunction::MultiColvarFunction(const ActionOptions& ao):
Action(ao),
MultiColvarBase(ao)
{
  // Read in the arguments
  std::string mname; std::vector<std::string> mlabs; parseVector("DATA",mlabs);
  log.printf("  using colvars calculated by actions ");
  for(unsigned i=0;i<mlabs.size();++i){
      log.printf("%s ",mlabs[i].c_str() );
      MultiColvarBase* mycolv = plumed.getActionSet().selectWithLabel<MultiColvarBase*>(mlabs[i]); 
      if(!mycolv) error("action labeled " + mlabs[i] + " does not exist or is not a multicolvar");
      // Check all base multicolvars are of same type
      if( i==0 ){ 
          mname = mycolv->getName();
          tvals.resize( mycolv->getNumberOfQuantities()-4 );
          if( mycolv->isPeriodic() ) error("multicolvar functions don't work with this multicolvar");
      } else {
          if( mname!=mycolv->getName() ) error("All input multicolvars must be of same type"); 
      }
      // Make sure we use low memory option in base colvar
      mycolv->setLowMemOption( usingLowMem() );
      // Checks for neighbor list
      if( mycolv->updateFreq>0 && updateFreq>0 ){
          if( updateFreq%mycolv->updateFreq!=0 ) error("update frequency must be a multiple of the update frequency in the base colvar");
      }
      // Add the dependency
      addDependency(mycolv);
      // And store the multicolvar base
      mybasemulticolvars.push_back( mycolv );
      // And track which variable stores each colvar
      for(unsigned j=0;j<mycolv->getFullNumberOfTasks();++j) colvar_label.push_back( i );
  }
  log.printf("\n");
}

void MultiColvarFunction::buildSymmetryFunctionLists(){
  if( mybasemulticolvars.size()>2 ) error("Found too many multicolvars in DATA specification. You can use either 1 or 2");

  // Make sure information is stored in the required multicolvars
  for(unsigned i=0;i<mybasemulticolvars.size();++i) mybasemulticolvars[i]->buildDataStashes();

  usespecies=true; ablocks.resize( 1 );
  for(unsigned i=0;i<mybasemulticolvars[0]->getFullNumberOfTasks();++i) addTaskToList( i );

  unsigned ntotal=0;
  for(unsigned i=0;i<mybasemulticolvars.size();++i){
      ntotal += mybasemulticolvars[i]->getFullNumberOfTasks();
  }
  unsigned k=0, start=0;
  current_atoms.resize( 1 + ntotal ); ablocks[0].resize( ntotal ); 
  for(unsigned i=0;i<mybasemulticolvars.size();++i){
      for(unsigned j=0;j<mybasemulticolvars[i]->getFullNumberOfTasks();++j){
          ablocks[0][k]=start + j; k++;
      }
      start += mybasemulticolvars[i]->getFullNumberOfTasks();
  }  
  // Copy lists of atoms involved from base multicolvars
  for(unsigned i=0;i<mybasemulticolvars.size();++i) mybasemulticolvars[i]->copyAtomListToFunction( this ); 
  // Do all setup stuff in MultiColvarBase
  setupMultiColvarBase();
}

void MultiColvarFunction::buildAtomListWithPairs( const bool& allow_intra_group ){
  if( !allow_intra_group && mybasemulticolvars.size()>2 ) error("only two input multicolvars allowed with this function"); 

  // Make sure information is stored in the required multicolvars
  for(unsigned i=0;i<mybasemulticolvars.size();++i) mybasemulticolvars[i]->buildDataStashes();
  
  usespecies=false; ablocks.resize(2); current_atoms.resize( 2 );
  if( !allow_intra_group && mybasemulticolvars.size()==2 ){
     nblock = mybasemulticolvars[0]->getFullNumberOfTasks();
     if( mybasemulticolvars[1]->getFullNumberOfTasks()>nblock ) nblock = mybasemulticolvars[1]->getFullNumberOfTasks();
    
     ablocks[0].resize( mybasemulticolvars[0]->getFullNumberOfTasks() );
     for(unsigned i=0;i<mybasemulticolvars[0]->getFullNumberOfTasks();++i) ablocks[0][i] = i;
     ablocks[1].resize( mybasemulticolvars[1]->getFullNumberOfTasks() ); unsigned istart = ablocks[0].size();
     for(unsigned i=0;i<mybasemulticolvars[1]->getFullNumberOfTasks();++i) ablocks[1][i] = istart + i;
     for(unsigned i=0;i<ablocks[0].size();++i){
         for(unsigned j=0;j<ablocks[1].size();++j) addTaskToList( i*nblock + j );
     }
  } else {
     nblock = 0; for(unsigned i=0;i<mybasemulticolvars.size();++i) nblock += mybasemulticolvars[i]->getFullNumberOfTasks();
     ablocks[0].resize( nblock ); ablocks[1].resize( nblock );
     for(unsigned i=0;i<nblock;++i){ ablocks[0][i] = i; ablocks[1][i] = i; }
     for(unsigned i=1;i<nblock;++i){
        for(unsigned j=0;j<i;++j) addTaskToList( i*nblock + j );
     }
  }
  // Copy lists of atoms involved from base multicolvars
  for(unsigned i=0;i<mybasemulticolvars.size();++i) mybasemulticolvars[i]->copyAtomListToFunction( this );
  // Do all setup stuff in MultiColvarBase
  setupMultiColvarBase();
}

void MultiColvarFunction::finishTaskListUpdate(){
  std::vector< std::vector<bool> > additionalTasks( mybasemulticolvars.size() );
  if( !contributorsAreUnlocked ){
      // Make a list of the tasks required 
      for(unsigned i=0;i<mybasemulticolvars.size();++i){
          additionalTasks[i].resize(mybasemulticolvars[i]->getFullNumberOfTasks(), false );
      }

      // Find what we are required to calculate from base multcolvar
      for(unsigned i=0;i<getCurrentNumberOfActiveTasks();++i){
          bool check=setupCurrentAtomList( getActiveTask(i) );
          plumed_assert( check );
          for(unsigned j=0;j<natomsper;++j){
             unsigned mmc = colvar_label[current_atoms[j]];
             unsigned tl = convertToLocalIndex( current_atoms[j], mmc );
             additionalTasks[mmc][tl] = true;
          }
      }
  } else {
      // Make a list of the tasks required 
      for(unsigned i=0;i<mybasemulticolvars.size();++i){
          additionalTasks[i].resize(mybasemulticolvars[i]->getFullNumberOfTasks(), true );
      }
  }

  // Deactivate all atoms
  all_atoms.deactivateAll(); unsigned start=0;
  for(unsigned i=0;i<mybasemulticolvars.size();++i){ 
     // Add these requirements in the base colvar
     mybasemulticolvars[i]->activateTheseTasks( additionalTasks[i] );
     // Redo preparation step for base colvar
     mybasemulticolvars[i]->prepare();
     // Copy the active atoms here
     mybasemulticolvars[i]->copyActiveAtomsToFunction( this, start );
     // Make sure we are updating the correct atoms
     start += mybasemulticolvars[i]->all_atoms.fullSize();
  }
  // Update all atoms array
  all_atoms.updateActiveMembers();
  // Request the atoms
  ActionAtomistic::requestAtoms( all_atoms.retrieveActiveList() );
  // Re-add all the dependencies
  for(unsigned i=0;i<mybasemulticolvars.size();++i) addDependency(mybasemulticolvars[i]);
  // Resize arrays in MultiColvarBase
  resizeLocalArrays();
  // Make numerical derivatives work
  if( checkNumericalDerivatives() ) numder_store.resize( getNumberOfComponents(), getNumberOfDerivatives() );
}

void MultiColvarFunction::calculate(){
  if( checkNumericalDerivatives() ){
     // Clear any derivatives in base colvar that were 
     // accumulated from previous calculations
     for(unsigned i=0;i<mybasemulticolvars.size();++i) mybasemulticolvars[i]->clearDerivatives(); 
     // And recalculate
     for(unsigned i=0;i<mybasemulticolvars.size();++i) mybasemulticolvars[i]->calculate();
     // Copy the box from the base multicolvar here
     unsigned maxb=mybasemulticolvars.size() - 1;
     changeBox( mybasemulticolvars[maxb]->getBox() );
  }
  setupLinkCells(); runAllTasks();
}

void MultiColvarFunction::calculateNumericalDerivatives( ActionWithValue* a ){
  unsigned pstart=0;
  for(unsigned i=0;i<mybasemulticolvars.size();++i){
     mybasemulticolvars[i]->calculateAtomicNumericalDerivatives( this, pstart );
     for(unsigned k=0;k<getNumberOfComponents();++k){
        Value* val=getPntrToComponent(k);
        for(unsigned j=0;j<3*mybasemulticolvars[i]->getNumberOfAtoms();++j){
           numder_store(k,pstart+j) = val->getDerivative(pstart + j);
        }
     }
     pstart += 3*mybasemulticolvars[i]->getNumberOfAtoms(); 
  }

  // Note numerical derivatives only work for virial if mybasemulticolvars.size()==1
  if( mybasemulticolvars.size()==1 ){
      for(unsigned k=0;k<getNumberOfComponents();++k){
         Value* val=getPntrToComponent(k);
         for(unsigned j=0;j<9;++j) numder_store(k,pstart+j) = val->getDerivative(pstart + j);
      }
  }

  // Now transfer to multicolvar
  clearDerivatives();
  for(unsigned j=0;j<getNumberOfComponents();++j){
     Value* val=getPntrToComponent(j);
     for(unsigned i=0;i<getNumberOfDerivatives();++i) val->addDerivative( i, numder_store(j,i) );
  }
}

void MultiColvarFunction::addStoredDerivative( const unsigned& jout, const unsigned& base_cv_no, const unsigned& base_index, const double& der ){
  plumed_dbg_assert( jout<getNumberOfQuantities() && base_cv_no<mybasemulticolvars.size() && base_index<mybasemulticolvars[base_cv_no]->getNumberOfDerivatives() );

  unsigned mbase = 3*mybasemulticolvars[base_cv_no]->getNumberOfAtoms(), tbase = 3*getNumberOfAtoms();
  if( base_index>=mbase ){
      // Add virial element
      unsigned jindex = base_index - mbase + tbase;
      addElementDerivative( jout*getNumberOfDerivatives() + jindex, der ); 
  } else {
      // Add atomic element
      unsigned offset=0; for(unsigned i=0;i<base_cv_no;++i) offset += 3*mybasemulticolvars[i]->getNumberOfAtoms();
      unsigned jindex = offset + base_index; 
      plumed_dbg_assert( jindex<3*getNumberOfAtoms() );
      addElementDerivative( jout*getNumberOfDerivatives() + jindex, der );
      unsigned iatom = std::floor( jindex / 3 );
      atoms_with_derivatives.activate( iatom );
  }
}

void MultiColvarFunction::updateActiveAtoms(){
  if( atoms_with_derivatives.updateComplete() ) return;
  atoms_with_derivatives.emptyActiveMembers();
  for(unsigned i=0;i<getNumberOfAtoms();++i) atoms_with_derivatives.updateIndex(i);
  atoms_with_derivatives.sortActiveList();
}

Vector MultiColvarFunction::calculateCentralAtomPosition(){
  Vector catom=getCentralAtom();
  atomsWithCatomDer.updateActiveMembers();
  return catom;
}

}
}

