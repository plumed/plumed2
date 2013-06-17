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

  // Checks for neighbor list
  if( mycolv->updateFreq>0 && updateFreq>0 ){
      if( updateFreq%mycolv->updateFreq!=0 ) error("update frequency must be a multiple of the update frequency in the base colvar");
  }

  // Retrieve the central atoms
  catoms = mycolv->getCentralAtoms();
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

double MultiColvarFunction::doCalculation( const unsigned& j ){
  double val=compute(j);
  atoms_with_derivatives.updateActiveMembers();
  return val;
}

void MultiColvarFunction::calculateNumericalDerivatives( ActionWithValue* a ){
  mycolv->calculateNumericalDerivatives( this );
}

unsigned MultiColvarFunction::getNumberOfAtomsInCentralAtomDerivatives(){
  return getNCAtomDerivatives()*mycolv->getNumberOfAtomsInCentralAtomDerivatives();
}

Vector MultiColvarFunction::calculateCentralAtomPosition(){
  Vector catom=getCentralAtom();
  atomsWithCatomDer.updateActiveMembers();
  return catom;
}

}
}

