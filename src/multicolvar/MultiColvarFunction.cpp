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
  for(unsigned i=0;i<mycolv->all_atoms.fullSize();++i){
     all_atoms.addIndexToList( mycolv->all_atoms(i) );
  }
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
  mycolv->unlockContributors();
  mpi_gatherActiveMembers( comm, colvar_atoms );
 
  // Store tasks that are currently active in base colvar
  unsigned nactive=mycolv->taskList.getNumberActive();
  std::vector<unsigned> active_tasks( nactive );
  for(unsigned i=0;i<nactive;++i) active_tasks[i]=mycolv->taskList[i];

  // Now get the tasks that are active from here
  mycolv->taskList.deactivateAll();
  for(unsigned i=0;i<taskList.getNumberActive();++i){
      unsigned n=taskList[i];
      for(unsigned j=0;j<colvar_atoms[n].getNumberActive();++j){
         mycolv->taskList.activate( colvar_atoms[n][j] );
      }
  }

  // Add in tasks that are active in base
  for(unsigned i=0;i<nactive;++i) mycolv->taskList.activate( active_tasks[i] );

  // Redo preparation step 
  mycolv->prepare();

  // And lock
  ActionWithVessel::lockContributors();
}

void MultiColvarFunction::resizeDynamicArrays(){
  // Copy what is active in the base MultiColvar here
  plumed_dbg_assert( all_atoms.fullSize()==mycolv->all_atoms.fullSize() );
  all_atoms.deactivateAll();
  for(unsigned i=0;i<mycolv->all_atoms.getNumberActive();++i){
     unsigned iatom=mycolv->all_atoms.linkIndex( i ); 
     all_atoms.activate( iatom );
  }
  all_atoms.updateActiveMembers();
  // Request the atoms
  ActionAtomistic::requestAtoms( all_atoms.retrieveActiveList() );
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

Vector MultiColvarFunction::calculateCentralAtomPosition(){
  Vector catom=getCentralAtom();
  atomsWithCatomDer.updateActiveMembers();
  return catom;
}

}
}

