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
#include "MultiColvarFunction.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "MultiColvarBase.h"

namespace PLMD {
namespace multicolvar { 

void MultiColvarFunction::registerKeywords( Keywords& keys ){
  MultiColvarBase::registerKeywords( keys );
  keys.add("compulsory","ARG","the label of the action that calculates the vectors we are interested in");
  keys.reserve("atoms","ONLY_ATOMS","only calculate this quantity for these central atoms");
}

MultiColvarFunction::MultiColvarFunction(const ActionOptions& ao):
Action(ao),
MultiColvarBase(ao)
{
  std::string mlab; parse("ARG",mlab);
  log.printf("  using vectors calculated by action %s\n",mlab.c_str() );
  mycolv = plumed.getActionSet().selectWithLabel<multicolvar::MultiColvarBase*>(mlab); 
  if(!mycolv) error("action labeled " + mlab + " does not exist or is not a multicolvar");
  addDependency(mycolv);
  // Make sure we use low memory option in base colvar
  mycolv->setLowMemOption( usingLowMem() );

  // Checks for neighbor list
  if( mycolv->updateFreq>0 && updateFreq>0 ){
      if( updateFreq%mycolv->updateFreq!=0 ) error("update frequency must be a multiple of the update frequency in the base colvar");
  }

  // Retrieve the central atoms
  catoms = mycolv->getCentralAtoms();

  if( keywords.exists("ONLY_ATOMS") ){
      std::vector<AtomNumber> atoms; parseAtomList("ONLY_ATOMS",atoms);
      if( atoms.size()>0 ){ 
         log.printf("  calculating function for atoms ");
         for(unsigned i=0;i<atoms.size();++i){
            log.printf("%d ",atoms[i].serial() );
            taskList.addIndexToList( mycolv->getInternalIndex(atoms[i]) );
         }
         log.printf("\n");
      } else {
         for(unsigned i=0;i<mycolv->taskList.fullSize();++i) taskList.addIndexToList( i );
      }
      usespecies=true; ablocks.resize(1); ablocks[0].resize( getNumberOfBaseFunctions() );
      for(unsigned i=0;i<getNumberOfBaseFunctions();++i) ablocks[0][i]=i; 
      current_atoms.resize( 1 + ablocks[0].size() );
  } else {
      usespecies=false; nblock=getNumberOfBaseFunctions(); ablocks.resize(2);
      for(unsigned i=0;i<2;++i) ablocks[i].resize(nblock);
      for(unsigned i=0;i<getNumberOfBaseFunctions();++i){ ablocks[0][i]=i; ablocks[1][i]=i; }
      for(unsigned i=1;i<getNumberOfBaseFunctions();++i){
         for(unsigned j=0;j<i;++j) taskList.addIndexToList( i*nblock + j );
      }
      current_atoms.resize( 2 );
  }
}

void MultiColvarFunction::completeSetup(){
  // Copy list of atoms involved from base multicolvar
  mycolv->copyAtomListToFunction( this );
  // Do all setup stuff in MultiColvarBase
  setupMultiColvarBase();
}

Vector MultiColvarFunction::getSeparation( const Vector& vec1, const Vector& vec2 ) const {
  return mycolv->getSeparation( vec1, vec2 );
}

void MultiColvarFunction::unlockContributors(){
  plumed_massert( mycolv->contributorsAreUnlocked,"contributors in base colvar are not unlocked"); 
  ActionWithVessel::unlockContributors();
}

void MultiColvarFunction::lockContributors(){
  // Make a list of the tasks required 
  std::vector<bool> additionalTasks( mycolv->getNumberOfTasks(), false );
  for(unsigned i=0;i<taskList.getNumberActive();++i){
      current=taskList[i];
      for(unsigned j=0;j<getNAtoms();++j) additionalTasks[ getAtomIndex(j) ]=true;
  }
  // And add these to the requirements in the base colvar
  mycolv->activateTheseTasks( additionalTasks ); 

  // Redo preparation step 
  mycolv->prepare();

  // And lock
  ActionWithVessel::lockContributors();
}

void MultiColvarFunction::resizeDynamicArrays(){
  mycolv->copyActiveAtomsToFunction( this );
  // Request the atoms
  requestAtoms();
  // Rerequest the dependency
  addDependency(mycolv);
}

void MultiColvarFunction::calculate(){
  if( checkNumericalDerivatives() ){
     // Clear any derivatives in base colvar that were 
     // accumulated from previous calculations
     mycolv->clearDerivatives(); 
     // And recalculate
     mycolv->calculate();
  }
  runAllTasks();
}

void MultiColvarFunction::calculateNumericalDerivatives( ActionWithValue* a ){
  mycolv->calculateNumericalDerivatives( this );
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

